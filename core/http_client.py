# -*- coding: utf-8 -*-
"""Resilient HTTP client with exponential backoff, rate limiting, and fallback.

Key features:
  - Token bucket rate limiter (configurable req/s per host)
  - Exponential backoff with jitter on transient failures
  - HTTP status classification: 429→backoff, 404→skip, 5xx→retry, 403→fail
  - Retry-After header support
  - NCBI_API_KEY env var for 10 req/s instead of 1 req/s
  - ENA ↔ NCBI fallback strategy
"""

import os
import time
import random
import logging
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError
from typing import Optional, Callable

log = logging.getLogger(__name__)


# ============================================================
# Token Bucket Rate Limiter
# ============================================================

class TokenBucket:
    """Simple token bucket for rate limiting.

    Example:
        limiter = TokenBucket(rate=3.0)  # 3 requests/second
        for url in urls:
            limiter.acquire()
            fetch(url)
    """

    def __init__(self, rate: float = 3.0, burst: int = 5):
        self.rate = rate          # tokens per second
        self.burst = burst        # max tokens
        self.tokens = burst       # current token count
        self.last_refill = time.monotonic()

    def acquire(self):
        """Block until a token is available."""
        now = time.monotonic()
        elapsed = now - self.last_refill
        self.tokens = min(self.burst, self.tokens + elapsed * self.rate)
        self.last_refill = now

        if self.tokens >= 1.0:
            self.tokens -= 1.0
            return

        # Not enough tokens — sleep until we have one
        wait = (1.0 - self.tokens) / self.rate
        time.sleep(wait)
        self.tokens = 0.0
        self.last_refill = time.monotonic()


# ============================================================
# Retry Configuration
# ============================================================

class RetryConfig:
    """Exponential backoff configuration."""

    def __init__(self,
                 max_retries: int = 3,
                 base_delay: float = 1.0,
                 max_delay: float = 60.0,
                 jitter: bool = True):
        self.max_retries = max_retries
        self.base_delay = base_delay
        self.max_delay = max_delay
        self.jitter = jitter

    def delay_for(self, attempt: int) -> float:
        """Calculate delay for a given retry attempt (0-indexed)."""
        delay = min(self.base_delay * (2 ** attempt), self.max_delay)
        if self.jitter:
            delay = random.uniform(0, delay)
        return delay

    def should_retry(self, attempt: int, status_code: Optional[int]) -> bool:
        """Determine if we should retry after this failure.

        Retry on: 429, 5xx, network errors, timeouts
        Skip on: 404 (Not Found), 403 (Forbidden after rate limit)
        """
        if attempt >= self.max_retries:
            return False
        if status_code is None:
            return True  # Network error — retry
        if status_code == 429:
            return True  # Rate limited — retry
        if 500 <= status_code < 600:
            return True  # Server error — retry
        return False


# ============================================================
# API Provider Configuration
# ============================================================

class ApiProvider:
    """Configuration for a specific API provider."""

    def __init__(self,
                 name: str,
                 host: str,
                 default_rate: float = 3.0,
                 api_key_env: str = "",
                 api_key_param: str = ""):
        self.name = name
        self.host = host
        self.default_rate = default_rate
        self.api_key_env = api_key_env
        self.api_key_param = api_key_param

    def get_api_key(self) -> Optional[str]:
        if self.api_key_env:
            return os.environ.get(self.api_key_env)
        return None

    def get_rate(self) -> float:
        """Get effective rate, boosted if API key is present."""
        if self.get_api_key():
            return self.default_rate * 3.3  # NCBI: 1→3 req/s without key, ~10 with
        return self.default_rate


# Provider instances
ENA_PROVIDER = ApiProvider("ENA", "www.ebi.ac.uk", default_rate=3.0)
NCBI_PROVIDER = ApiProvider(
    "NCBI", "eutils.ncbi.nlm.nih.gov",
    default_rate=1.0,  # 1 req/s without key (NCBI guidelines)
    api_key_env="NCBI_API_KEY",
    api_key_param="api_key",
)


# ============================================================
# Resilient HTTP Client
# ============================================================

class RateLimitedSession:
    """HTTP client with rate limiting and retry logic.

    Usage:
        session = RateLimitedSession(provider=ENA_PROVIDER)
        text = session.get("https://www.ebi.ac.uk/ena/portal/api/filereport?...")
    """

    def __init__(self,
                 provider: Optional[ApiProvider] = None,
                 retry_config: Optional[RetryConfig] = None,
                 timeout: int = 120,
                 user_agent: str = "PhenotypeCuration/1.0"):
        self.provider = provider
        self.rate_limiter = TokenBucket(rate=provider.get_rate() if provider else 3.0)
        self.retry_config = retry_config or RetryConfig()
        self.timeout = timeout
        self.user_agent = user_agent
        self._limiter_map: dict[str, TokenBucket] = {}  # per-host limiters

    def _get_limiter(self, host: str) -> Optional[TokenBucket]:
        """Get or create a rate limiter for a host."""
        if host not in self._limiter_map:
            if self.provider and host == self.provider.host:
                self._limiter_map[host] = self.rate_limiter
            else:
                self._limiter_map[host] = TokenBucket(rate=3.0)
        return self._limiter_map[host]

    def _build_request(self, url: str) -> Request:
        """Build a Request with appropriate headers and API key injection."""
        # Inject NCBI API key if applicable
        if self.provider and self.provider.host in url:
            api_key = self.provider.get_api_key()
            if api_key and '?' in url and self.provider.api_key_param not in url:
                url = f"{url}&{self.provider.api_key_param}={api_key}"

        headers = {'User-Agent': self.user_agent}
        return Request(url, headers=headers)

    def get(self, url: str,
            retry_config: Optional[RetryConfig] = None,
            timeout: Optional[int] = None) -> str:
        """GET with retry, backoff, and rate limiting.

        Returns response body as string.
        Raises RuntimeError if all retries exhausted.
        """
        retry = retry_config or self.retry_config
        t = timeout or self.timeout
        last_error = None

        for attempt in range(retry.max_retries + 1):
            # Rate limit
            # Extract host from URL
            host = url.split('/')[2] if '//' in url else 'unknown'
            limiter = self._get_limiter(host)
            if limiter:
                limiter.acquire()

            try:
                req = self._build_request(url)
                with urlopen(req, timeout=t) as resp:
                    return resp.read().decode('utf-8')

            except HTTPError as e:
                status = e.code

                # 404 — not found, don't retry
                if status == 404:
                    raise RuntimeError(
                        f"Resource not found: {url} (HTTP 404)") from e

                # 403 — forbidden
                if status == 403:
                    raise RuntimeError(
                        f"Access forbidden: {url} (HTTP 403). "
                        f"Try setting {self.provider.api_key_env if self.provider else ''} env var.") from e

                # 429 — rate limited
                if status == 429:
                    retry_after = e.headers.get('Retry-After', '')
                    delay = int(retry_after) if retry_after.isdigit() else retry.delay_for(attempt)
                    log.warning(f"Rate limited (429) on {host}, waiting {delay:.0f}s")
                    time.sleep(delay)
                    last_error = e
                    continue

                # Other HTTP errors (4xx, 5xx)
                if retry.should_retry(attempt, status):
                    delay = retry.delay_for(attempt)
                    log.warning(f"HTTP {status} on {host}, retry {attempt+1}/{retry.max_retries} in {delay:.1f}s")
                    time.sleep(delay)
                    last_error = e
                    continue
                else:
                    raise RuntimeError(f"HTTP {status} on {url}") from e

            except (URLError, TimeoutError, OSError) as e:
                if retry.should_retry(attempt, None):
                    delay = retry.delay_for(attempt)
                    log.warning(f"Network error on {host}, retry {attempt+1}/{retry.max_retries} in {delay:.1f}s: {e}")
                    time.sleep(delay)
                    last_error = e
                    continue
                else:
                    raise RuntimeError(f"Network error after {attempt+1} attempts: {url}") from e

        raise RuntimeError(
            f"All {retry.max_retries+1} attempts failed for {url}: {last_error}")


# ============================================================
# Convenience: ENA + NCBI sessions
# ============================================================

def ena_session(timeout: int = 120) -> RateLimitedSession:
    """Create a rate-limited ENA API session."""
    return RateLimitedSession(provider=ENA_PROVIDER, timeout=timeout)


def ncbi_session() -> RateLimitedSession:
    """Create a rate-limited NCBI E-utilities session.

    Set NCBI_API_KEY env var for 10 req/s access.
    """
    return RateLimitedSession(provider=NCBI_PROVIDER, timeout=60)


# ============================================================
# Fallback Resolver
# ============================================================

class FallbackFetcher:
    """Try primary API, fall back to secondary on failure."""

    def __init__(self,
                 primary: Callable[[], str],
                 fallback: Optional[Callable[[], str]] = None,
                 name: str = ""):
        self.primary = primary
        self.fallback = fallback
        self.name = name

    def fetch(self) -> tuple[str, str]:  # (content, source_label)
        """Returns (content, source) where source is 'primary' or 'fallback'."""
        try:
            content = self.primary()
            return content, "primary"
        except Exception as e:
            log.warning(f"Primary failed for {self.name}: {e}")
            if self.fallback:
                try:
                    content = self.fallback()
                    return content, "fallback"
                except Exception as fe:
                    log.error(f"Fallback also failed for {self.name}: {fe}")
                    raise RuntimeError(
                        f"Both primary and fallback failed for {self.name}") from fe
            raise