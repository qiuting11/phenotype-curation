# -*- coding: utf-8 -*-
"""Phase 3: Paper-verified study design classification (LLM-assisted).

This is the ONLY phase that requires LLM/AI judgment. It:
1. Reads project_registry.csv, filters to pending-review projects
2. Fetches PMC full text via NCBI E-utilities (when available)
3. Generates structured prompts for study design classification
4. Parses LLM output back into InterventionRegistryRow (schema-validated)
5. Supports batch processing with incremental output (rewrite-safe)

Input:  project_registry.csv (Phase 2)
Output: intervention_registry.csv (Phase 3 → drives Phase 4 build)

Architecture:
    The prompt template encodes the full decision framework from the skill doc:
    study_design taxonomy, intervention_type taxonomy, one-time vs ongoing
    distinction, evidence quality rating. The LLM reads paper Methods and
    returns structured JSON validated against intervention_registry_row.schema.json.
"""

import csv
import json
import os
import re
import sys
import time
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET
from dataclasses import asdict
from pathlib import Path
from typing import Optional

_PARENT = Path(__file__).resolve().parent.parent
if str(_PARENT) not in sys.path:
    sys.path.insert(0, str(_PARENT))

from core.types import InterventionRegistryRow, EvidenceQuality
from core.schemas import validate_row


# ============================================================
# PMC Full-Text Fetcher
# ============================================================

NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# NCBI API key from env for higher rate limit (10 req/s vs 1 req/s)
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")


def _eutils_url(endpoint: str, **params) -> str:
    """Build NCBI E-utilities URL with API key if available."""
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    qs = "&".join(f"{k}={urllib.request.quote(str(v))}" for k, v in params.items())
    return f"{NCBI_EUTILS_BASE}/{endpoint}?{qs}"


def fetch_pmc_full_text(pmcid: str, timeout: int = 30) -> Optional[str]:
    """Fetch PMC full-text XML and extract Methods section.

    Args:
        pmcid: PMC ID (e.g. 'PMC7386241').
        timeout: HTTP timeout in seconds.

    Returns:
        Extracted methods text, or None if unavailable.
    """
    pmcid = pmcid.strip()
    if not pmcid.startswith("PMC"):
        pmcid = f"PMC{pmcid}"

    url = _eutils_url("efetch.fcgi", db="pmc", id=pmcid, rettype="xml", retmode="xml")

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PhenotypeCuration/0.1"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            xml_text = resp.read().decode("utf-8")
    except urllib.error.HTTPError as e:
        print(f"    HTTP {e.code} fetching {pmcid}")
        return None
    except Exception as e:
        print(f"    Error fetching {pmcid}: {e}")
        return None

    # Parse XML
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return None

    # Extract body text, focusing on Methods-like sections
    sections = []
    namespaces = {"ns": "http://www.w3.org/1999/xhtml"}

    # Try to find <body> in PMC XML (may be namespaced)
    body = root.find(".//body") or root.find(".//{http://www.w3.org/1999/xhtml}body")
    if body is None:
        # Maybe the whole document is the body
        body = root

    for sec in body.iter():
        tag = sec.tag.split("}")[-1] if "}" in sec.tag else sec.tag
        if tag in ("sec", "section"):
            title_el = sec.find("title") or sec.find("{http://www.w3.org/1999/xhtml}title")
            title = "".join(title_el.itertext()).strip() if title_el is not None else ""
            text = "".join(sec.itertext()).strip()

            # Look for Methods-like sections
            title_lower = title.lower()
            is_methods = any(kw in title_lower for kw in [
                "method", "study design", "study population",
                "sample collection", "intervention", "treatment",
                "participant", "patient", "subject", "experimental design",
                "protocol", "procedure", "trial design",
                "study design and participants", "materials and methods",
                "patients and methods",
            ])

            if is_methods or (not title and "method" in text[:200].lower()):
                sections.append(f"## {title}\n{text}")
            elif title_lower in ("abstract", "introduction", "background", "results",
                                 "discussion", "conclusion", "references", "acknowledgements",
                                 "supplementary", "figure", "table"):
                # Skip non-Methods sections
                pass
            elif title and len(text) > 100:
                # Include uncategorized sections (might be sub-sections of Methods)
                sections.append(f"## {title}\n{text}")

    if not sections:
        return None

    # Truncate to reasonable size (~8000 chars for LLM context)
    full_text = "\n\n".join(sections)
    if len(full_text) > 12000:
        full_text = full_text[:12000] + "\n\n[... truncated for length ...]"

    return full_text


def fetch_pubmed_abstract(pmid: str, timeout: int = 15) -> Optional[str]:
    """Fetch PubMed abstract as fallback when PMC full text unavailable.

    Args:
        pmid: PubMed ID (digits only).
        timeout: HTTP timeout.

    Returns:
        Abstract text, or None.
    """
    pmid = pmid.strip()
    url = _eutils_url("efetch.fcgi", db="pubmed", id=pmid, rettype="abstract", retmode="text")

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PhenotypeCuration/0.1"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            text = resp.read().decode("utf-8")
            return text[:5000] if len(text) > 5000 else text
    except Exception as e:
        print(f"    Error fetching abstract for PMID {pmid}: {e}")
        return None


def resolve_pmcid_from_pmid(pmid: str, timeout: int = 10) -> Optional[str]:
    """Convert PMID to PMCID via NCBI ID converter.

    Args:
        pmid: PubMed ID.
        timeout: HTTP timeout.

    Returns:
        PMCID string (e.g. 'PMC7386241') or None.
    """
    pmid = pmid.strip()
    url = _eutils_url("elink.fcgi", dbfrom="pubmed", db="pmc", id=pmid, linkname="pubmed_pmc")

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PhenotypeCuration/0.1"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            xml_text = resp.read().decode("utf-8")
    except Exception:
        return None

    root = ET.fromstring(xml_text)

    # NCBI elink returns LinkSetDb/Id elements with PMC IDs
    # Sometimes IDs are returned as plain PubMed IDs (need to use idconvert instead)
    for link in root.iter("Id"):
        if link.text:
            # Check if it looks like a PMC ID
            if link.text.startswith("PMC"):
                return link.text
            # Check if it might be a PMCID with just the numeric part
            if link.text.isdigit() and len(link.text) >= 7:
                # Try to confirm via idconvert API
                pmcid = _resolve_via_idconvert(pmid, timeout)
                if pmcid:
                    return pmcid

    # Fallback: try NCBI ID converter API
    return _resolve_via_idconvert(pmid, timeout)


def _resolve_via_idconvert(pmid: str, timeout: int = 10) -> Optional[str]:
    """Try NCBI ID Converter API (pubmed → pmc)."""
    pmid = pmid.strip()
    url = _eutils_url("epost.fcgi", db="pubmed", id=pmid)

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PhenotypeCuration/0.1"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            pass  # We use efetch/pmc directly instead
    except Exception:
        pass

    # Simpler approach: just try to fetch from PMC directly
    direct_url = _eutils_url("esearch.fcgi", db="pmc", term=f"{pmid}[PMID]", retmax=1)
    try:
        req = urllib.request.Request(direct_url, headers={"User-Agent": "PhenotypeCuration/0.1"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            xml_text = resp.read().decode("utf-8")
        root = ET.fromstring(xml_text)
        for id_el in root.iter("Id"):
            if id_el.text:
                return f"PMC{id_el.text}"
    except Exception:
        pass

    return None


# ============================================================
# Classification Prompt Template
# ============================================================

CLASSIFICATION_SYSTEM_PROMPT = """You are a study-design classifier for metagenomics phenotype curation.
Your job: read the paper Methods section and classify the study into structured fields.

## Classification Framework

### study_design (choose exactly one)
- `observational_cohort` — No study-administered intervention; samples collected for monitoring/observation
- `RCT` — Randomized controlled trial with distinct intervention and control arms
- `RCT_crossover` — Participants cross over between intervention and control arms
- `dose_escalation` — All participants receive intervention at escalating dose tiers

### intervention_type (choose exactly one)
- `observational` — No intervention (use with study_design=observational_cohort)
- `immunotherapy` — Anti-PD-1, anti-CTLA-4, checkpoint inhibitors
- `diet` — Mediterranean diet, low-FODMAP, etc. (whole-diet interventions)
- `diet_prebiotic` — FOS, PDX, inulin, GOS, XOS (isolated dietary glycans/fibers)
- `fmt` — Fecal microbiota transplantation
- `probiotic` — Live bacterial strains (BB-12, LGG, etc.)
- `vaginal_seeding` — Maternal vaginal microbiota transfer at birth

### is_one_time (YES/NO — the critical distinction)
**ONE-TIME ("yes"):** Intervention delivered at a single moment. After that moment, the
participant is in POST-intervention state, even if biological effects persist.
Examples: single-dose FMT capsule, vaginal seeding at birth, alloHSCT transplant.

**ONGOING ("no"):** Intervention maintained over a period. Samples during that period
are "during intervention." Examples: daily probiotic for 10 days, weekly FMT for 8
weeks, 8-week Mediterranean diet, dose-escalation feeding phases.

### arms
Semicolon-separated study arm names (e.g., "FMT;Placebo" or "FOS_Cohort_1;FOS_Cohort_2;FOS_Cohort_3").

### control_type
What is the control condition? `placebo`, `habitual_diet`, `run_in_baseline`,
`within_subject`, `N/A` (for observational).

### blinding
`double_blind`, `single_blind`, `open_label`, or empty string.

### sample_timing_schedule
Human-readable sampling schedule (e.g., "run_in_2wk, feeding_4wk, washout_4wk").

### days_since_intervention
Numeric days from intervention start for each timepoint as a JSON array
(e.g., "[0, 14, 28, 56]"). Only for ongoing interventions.

### classification_confidence
- `verified` — Day-level timeline from paper supplement/tables
- `good` — Clear Methods description, unambiguous mapping
- `moderate` — Some ambiguity in timing or design
- `low` — Paper details unclear or cross-sectional mixed-timing

### notes
Any design details that affect interpretation but don't fit the structured fields.

## Important Rules
1. If the study is observational, ALL samples are is_intervention="no"
2. If intervention is one-time, ALL post-baseline samples are is_intervention="no" (they are follow-up)
3. If intervention is ongoing, only samples DURING the active treatment window are is_intervention="yes"
4. A control group receiving placebo in an RCT does NOT make it "observational" or change is_intervention
5. PMID must be digits only (no "PMC" prefix)

Return ONLY valid JSON matching this schema:
{
  "study_design": "RCT",
  "intervention_type": "fmt",
  "is_one_time": "no",
  "arms": "FMT;Placebo",
  "control_type": "placebo",
  "randomization": "yes",
  "blinding": "double_blind",
  "sample_timing": "run_in_2wk, feeding_4wk, washout_4wk",
  "sample_timing_schedule": "baseline_day_0, week_2, week_4, week_8, week_12",
  "days_since_intervention": "[0, 14, 28, 56, 84]",
  "classification_confidence": "good",
  "notes": "Key design details..."
}
"""


def build_classification_prompt(
    bioproject: str,
    pmid: str,
    pmcid: str = "",
    paper_text: str = "",
    metadata_signals: dict = None,
) -> str:
    """Build a classification prompt for a single project.

    Args:
        bioproject: Project accession.
        pmid: PubMed ID.
        pmcid: PMC ID (optional, used for fetching full text).
        paper_text: Pre-fetched paper Methods text. If empty, prompt includes
                    instructions to fetch it.
        metadata_signals: Auto-detected signals from Phase 2 registry
                          (raw_group_values, study_design_candidate, etc.).

    Returns:
        Formatted prompt string.
    """
    metadata_str = ""
    if metadata_signals:
        signals = []
        for k, v in sorted(metadata_signals.items()):
            if v:
                signals.append(f"  {k}: {v}")
        if signals:
            metadata_str = "## Metadata Signals (auto-detected from sample names)\n" + "\n".join(signals)

    paper_section = ""
    if paper_text:
        # Truncate if very long
        text = paper_text[:15000]
        paper_section = f"## Paper Text (Methods)\n```\n{text}\n```"
    else:
        paper_section = (
            "## Paper Text\n"
            "NOT PROVIDED. Use the following to locate the paper:\n"
            f"- PubMed: https://pubmed.ncbi.nlm.nih.gov/{pmid}/\n"
            + (f"- PMC Full Text: https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/\n" if pmcid else "")
            + "Read the Methods section (especially Study Design, Sample Collection, Intervention) "
            "and return the classification."
        )

    return f"""## Project: {bioproject}
## PMID: {pmid}
{f"## PMCID: {pmcid}" if pmcid else ""}

{metadata_str}

{paper_section}

---
Classify this study using the framework above. Return ONLY valid JSON."""


# ============================================================
# Output Parser
# ============================================================

def parse_classification_output(
    json_text: str,
    bioproject: str,
    pmid: str = "",
    pmcid: str = "",
) -> tuple[Optional[InterventionRegistryRow], list[str]]:
    """Parse LLM JSON output into a validated InterventionRegistryRow.

    Args:
        json_text: Raw LLM output (should contain JSON).
        bioproject: Project accession for the row.
        pmid: PubMed ID.
        pmcid: PMC ID.

    Returns:
        (InterventionRegistryRow or None, list of validation errors)
    """
    errors = []

    # Extract JSON from possibly messy LLM output
    json_text = json_text.strip()

    # Try to find JSON block
    m = re.search(r'\{[\s\S]*\}', json_text)
    if not m:
        return None, ["No JSON object found in output"]

    json_str = m.group(0)

    try:
        data = json.loads(json_str)
    except json.JSONDecodeError as e:
        return None, [f"JSON parse error: {e}"]

    # Build registry row
    row = InterventionRegistryRow(
        bioproject=bioproject,
        study_design=data.get("study_design", ""),
        intervention_type=data.get("intervention_type", ""),
        arms=data.get("arms", ""),
        control_type=data.get("control_type", ""),
        randomization=data.get("randomization", "no"),
        blinding=data.get("blinding", ""),
        sample_timing=data.get("sample_timing", ""),
        primary_pmid=pmid,
        primary_pmc=pmcid,
        notes=data.get("notes", ""),
        is_one_time=data.get("is_one_time", ""),
        sample_timing_schedule=data.get("sample_timing_schedule", ""),
        days_since_intervention=data.get("days_since_intervention", ""),
        classification_confidence=data.get("classification_confidence", "good"),
    )

    # Validate against schema
    row_dict = asdict(row)
    schema_errors = validate_row(row_dict, "intervention_registry_row")

    if schema_errors:
        for e in schema_errors:
            errors.append(str(e))
        # Don't return None — partial results are still useful for human review
        # But mark confidence down
        if row.classification_confidence == "good":
            row.classification_confidence = "moderate"

    return row, errors


# ============================================================
# Batch Classifier
# ============================================================

def classify_projects(
    registry_path: str | Path,
    out_path: Optional[str | Path] = None,
    fetch_papers: bool = True,
    paper_dir: Optional[str | Path] = None,
) -> list[dict]:
    """Generate classification prompts for all pending projects.

    This function does NOT call an LLM — it prepares prompts and optionally
    fetches paper text. The actual LLM call is done by a human or an AI agent
    following the protocol.

    Args:
        registry_path: Path to project_registry.csv.
        out_path: Output path for intervention_registry.csv (incremental).
        fetch_papers: Try to fetch PMC full text via NCBI E-utilities.
        paper_dir: Directory to cache downloaded paper text files.

    Returns:
        List of prompt dicts (one per pending project), each containing:
        {bioproject, pmid, pmcid, prompt, paper_text (if fetched), metadata_signals}
    """
    registry_path = Path(registry_path)
    paper_dir = Path(paper_dir) if paper_dir else Path("paper_texts")
    paper_dir.mkdir(parents=True, exist_ok=True)

    # Load registry
    with open(registry_path, "r", encoding="utf-8-sig") as f:
        registry_rows = list(csv.DictReader(f))

    # Filter to pending
    pending = [
        r for r in registry_rows
        if r.get("manual_review_status", "").startswith("pending")
    ]

    if not pending:
        print("No pending projects found. All projects already verified.")
        return []

    print(f"Found {len(pending)} pending projects out of {len(registry_rows)} total")

    prompts = []
    existing_rows: dict[str, InterventionRegistryRow] = {}

    # Load existing intervention registry if it exists (for incremental update)
    if out_path and Path(out_path).exists():
        with open(out_path, "r", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                bp = row.get("bioproject", "").strip()
                if bp:
                    existing_rows[bp] = InterventionRegistryRow(
                        bioproject=bp,
                        study_design=row.get("study_design", ""),
                        intervention_type=row.get("intervention_type", ""),
                        arms=row.get("arms", ""),
                        control_type=row.get("control_type", ""),
                        randomization=row.get("randomization", "no"),
                        blinding=row.get("blinding", ""),
                        sample_timing=row.get("sample_timing", ""),
                        primary_pmid=row.get("primary_pmid", ""),
                        primary_pmc=row.get("primary_pmc", ""),
                        notes=row.get("notes", ""),
                        is_one_time=row.get("is_one_time", ""),
                        sample_timing_schedule=row.get("sample_timing_schedule", ""),
                        days_since_intervention=row.get("days_since_intervention", ""),
                        classification_confidence=row.get("classification_confidence", "good"),
                    )
        print(f"  Loaded {len(existing_rows)} existing classifications (will skip if present)")

    for r in pending:
        bp = r.get("bioproject_accession", r.get("bioproject", "")).strip()
        pmid = r.get("primary_pmid", "").strip()

        if not bp or not pmid:
            print(f"  SKIP: missing bioproject or PMID")
            continue

        # Skip if already classified
        if bp in existing_rows:
            print(f"  SKIP {bp}: already classified")
            continue

        print(f"\n  Processing {bp} (PMID: {pmid})...")

        # Resolve PMCID from PMID
        pmcid = r.get("primary_pmc", "").strip()
        if not pmcid and fetch_papers:
            pmcid = resolve_pmcid_from_pmid(pmid) or ""
            if pmcid:
                print(f"    Resolved PMC: {pmcid}")
                time.sleep(0.35)  # NCBI rate limit

        # Fetch paper text
        paper_text = ""

        # Check cache first
        cache_path = paper_dir / f"{bp}_paper.txt"
        if cache_path.exists():
            paper_text = cache_path.read_text(encoding="utf-8")
            print(f"    Using cached paper text ({len(paper_text)} chars)")

        if not paper_text and fetch_papers and pmcid:
            print(f"    Fetching PMC full text for {pmcid}...")
            paper_text = fetch_pmc_full_text(pmcid) or ""
            if paper_text:
                cache_path.write_text(paper_text, encoding="utf-8")
                print(f"    Extracted {len(paper_text)} chars, cached to {cache_path}")
            else:
                # Fallback to abstract
                print(f"    PMC full text unavailable, trying PubMed abstract...")
                paper_text = fetch_pubmed_abstract(pmid) or ""
                if paper_text:
                    print(f"    Using abstract ({len(paper_text)} chars)")

        # Build metadata signals
        metadata_signals = {
            "study_design_candidate": r.get("study_design_candidate", ""),
            "raw_group_values": r.get("raw_group_values", ""),
            "has_intervention_metadata": r.get("has_intervention_metadata", ""),
            "has_timepoint_metadata": r.get("has_timepoint_metadata", ""),
        }

        # Build prompt
        prompt = build_classification_prompt(
            bioproject=bp,
            pmid=pmid,
            pmcid=pmcid,
            paper_text=paper_text,
            metadata_signals=metadata_signals,
        )

        prompts.append({
            "bioproject": bp,
            "pmid": pmid,
            "pmcid": pmcid,
            "prompt": prompt,
            "paper_text": paper_text,
            "metadata_signals": metadata_signals,
        })

        time.sleep(0.35)  # NCBI rate limit

    # Write prompts JSON for batch processing
    prompts_path = Path(out_path).with_suffix(".prompts.json") if out_path else Path("classify_prompts.json")
    with open(prompts_path, "w", encoding="utf-8") as f:
        json.dump([
            {k: v for k, v in p.items() if k != "prompt"}
            | {"prompt": p["prompt"][:500] + "..." if len(p["prompt"]) > 500 else p["prompt"]}
            for p in prompts
        ], f, indent=2, ensure_ascii=False)
    print(f"\nGenerated {len(prompts)} prompts → {prompts_path}")

    return prompts


def apply_classifications(
    prompts: list[dict],
    llm_outputs: dict[str, str],    # bioproject → LLM JSON response
    out_path: str | Path,
) -> dict[str, InterventionRegistryRow]:
    """Apply LLM classification outputs to build intervention_registry.csv.

    Args:
        prompts: List of prompt dicts from classify_projects().
        llm_outputs: Dict mapping bioproject → LLM JSON response text.
        out_path: Output CSV path.

    Returns:
        Dict of bioproject → InterventionRegistryRow.
    """
    out_path = Path(out_path)

    # Load existing
    existing: dict[str, InterventionRegistryRow] = {}
    if out_path.exists():
        with open(out_path, "r", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                bp = row.get("bioproject", "").strip()
                if bp:
                    existing[bp] = InterventionRegistryRow(
                        bioproject=bp,
                        study_design=row.get("study_design", ""),
                        intervention_type=row.get("intervention_type", ""),
                        arms=row.get("arms", ""),
                        control_type=row.get("control_type", ""),
                        randomization=row.get("randomization", "no"),
                        blinding=row.get("blinding", ""),
                        sample_timing=row.get("sample_timing", ""),
                        primary_pmid=row.get("primary_pmid", ""),
                        primary_pmc=row.get("primary_pmc", ""),
                        notes=row.get("notes", ""),
                        is_one_time=row.get("is_one_time", ""),
                        sample_timing_schedule=row.get("sample_timing_schedule", ""),
                        days_since_intervention=row.get("days_since_intervention", ""),
                        classification_confidence=row.get("classification_confidence", "good"),
                    )

    # Apply new classifications
    results: dict[str, InterventionRegistryRow] = dict(existing)
    errors_by_project: dict[str, list[str]] = {}

    for p in prompts:
        bp = p["bioproject"]
        output = llm_outputs.get(bp, "")
        if not output:
            print(f"  WARN: {bp} has no LLM output — skipping")
            continue

        row, errs = parse_classification_output(
            output, bp,
            pmid=p.get("pmid", ""),
            pmcid=p.get("pmcid", ""),
        )

        if row:
            results[bp] = row
            print(f"  OK: {bp} → {row.study_design}/{row.intervention_type} "
                  f"(confidence={row.classification_confidence})")
        if errs:
            errors_by_project[bp] = errs
            for e in errs[:3]:
                print(f"    Validation: {e}")

    # Write output CSV
    fieldnames = [
        "bioproject", "study_design", "intervention_type", "arms",
        "control_type", "randomization", "blinding", "sample_timing",
        "primary_pmid", "primary_pmc", "notes",
        "is_one_time", "sample_timing_schedule", "days_since_intervention",
        "classification_confidence",
    ]

    with open(out_path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for bp in sorted(results.keys()):
            row = results[bp]
            writer.writerow(asdict(row))

    print(f"\nWrote {len(results)} rows to {out_path}")
    if errors_by_project:
        print(f"\n{sum(len(v) for v in errors_by_project.values())} validation issues "
              f"across {len(errors_by_project)} projects")

    return results


# ============================================================
# CLI
# ============================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Phase 3: Paper-verified study design classification"
    )
    parser.add_argument(
        "--registry",
        required=True,
        help="Path to project_registry.csv",
    )
    parser.add_argument(
        "--out",
        default="intervention_registry.csv",
        help="Output CSV path",
    )
    parser.add_argument(
        "--no-fetch",
        action="store_true",
        help="Skip PMC full text fetch (only generate prompts)",
    )
    parser.add_argument(
        "--paper-dir",
        default="paper_texts",
        help="Directory to cache paper text files",
    )
    parser.add_argument(
        "--apply",
        help="Path to JSON file mapping bioproject → LLM output (apply classifications)",
    )

    args = parser.parse_args()

    if args.apply:
        # Apply mode: read LLM outputs and write registry
        with open(args.apply, "r", encoding="utf-8") as f:
            llm_outputs = json.load(f)
        prompts_path = Path(args.out).with_suffix(".prompts.json")
        prompts = []
        if prompts_path.exists():
            with open(prompts_path, "r", encoding="utf-8") as f:
                prompts = json.load(f)
        apply_classifications(prompts, llm_outputs, args.out)
    else:
        # Generate mode: prepare prompts and fetch papers
        classify_projects(
            registry_path=args.registry,
            out_path=args.out,
            fetch_papers=not args.no_fetch,
            paper_dir=args.paper_dir,
        )