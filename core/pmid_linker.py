# -*- coding: utf-8 -*-
"""PMID → BioProject linker (Gap #2).

Given a PubMed ID, find associated BioProject accessions via:
  - NCBI E-utilities elink (pubmed → bioproject)
  - NCBI E-utilities elink (pubmed → biosample → bioproject, if direct link fails)

This bridges the gap between "I have a paper" and "I have ENA data."
"""

import os
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET
from typing import Optional

NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")


def _eutils_url(endpoint: str, **params) -> str:
    """Build NCBI E-utilities URL with API key if available."""
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    qs = "&".join(f"{k}={urllib.request.quote(str(v))}" for k, v in params.items())
    return f"{NCBI_EUTILS_BASE}/{endpoint}?{qs}"


def pmid_to_bioprojects(pmid: str, timeout: int = 15) -> list[str]:
    """Find BioProject accessions associated with a PubMed ID.

    Strategy:
      1. pubmed → bioproject (direct elink)
      2. pubmed → biosample → bioproject (two-hop elink)
      3. pubmed → sra → bioproject (via SRA experiments)
      4. ENA Portal API text search by PMID

    Note: NCBI E-utilities cross-database links are not always maintained.
    Many papers (especially non-PMC) have no pubmed→bioproject links.
    In those cases, manual lookup via SRA Run Selector or paper Data
    Availability statement is required.

    Args:
        pmid: PubMed ID (e.g., '32620143').
        timeout: HTTP timeout per API call.

    Returns:
        List of unique BioProject accessions (e.g., ['PRJNA594393']).
    """
    pmid = pmid.strip()
    bioprojects: set[str] = set()

    # Strategy 1: Direct pubmed → bioproject elink
    direct = _elink_direct(pmid, timeout)
    if direct:
        bioprojects.update(direct)
        if bioprojects:
            return sorted(bioprojects)

    # Strategy 2: pubmed → biosample → bioproject
    via_biosample = _elink_via_biosample(pmid, timeout)
    if via_biosample:
        bioprojects.update(via_biosample)

    # Strategy 3: pubmed → SRA → bioproject
    via_sra = _elink_via_sra(pmid, timeout)
    if via_sra:
        bioprojects.update(via_sra)

    # Strategy 4: ENA Portal API text search
    via_ena = _search_ena_by_pmid(pmid, timeout)
    if via_ena:
        bioprojects.update(via_ena)

    return sorted(bioprojects)


def _elink_direct(pmid: str, timeout: int) -> list[str]:
    """Direct elink: pubmed → bioproject."""
    url = _eutils_url("elink.fcgi", dbfrom="pubmed", db="bioproject", id=pmid,
                       linkname="pubmed_bioproject")

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PhenotypeCuration/0.1"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            xml_text = resp.read().decode("utf-8")
    except Exception:
        return []

    root = ET.fromstring(xml_text)
    bioprojects = []

    for link_set in root.iter("LinkSetDb"):
        db_to = link_set.findtext("DbTo", "")
        if "bioproject" in db_to.lower():
            for link in link_set.iter("Id"):
                if link.text:
                    bioprojects.append(link.text.strip())

    return bioprojects


def _elink_via_biosample(pmid: str, timeout: int) -> list[str]:
    """Two-hop: pubmed → biosample → bioproject."""
    bioprojects: set[str] = set()

    # Step 1: pubmed → biosample
    url1 = _eutils_url("elink.fcgi", dbfrom="pubmed", db="biosample", id=pmid)

    try:
        req = urllib.request.Request(url1, headers={"User-Agent": "PhenotypeCuration/0.1"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            xml_text = resp.read().decode("utf-8")
    except Exception:
        return []

    root = ET.fromstring(xml_text)
    biosample_ids: list[str] = []

    for link_set in root.iter("LinkSetDb"):
        db_to = link_set.findtext("DbTo", "")
        if "biosample" in db_to.lower():
            for link in link_set.iter("Id"):
                if link.text:
                    biosample_ids.append(link.text.strip())

    if not biosample_ids:
        return []

    # Step 2: biosample → bioproject (batch with first 100)
    sample_batch = ",".join(biosample_ids[:100])
    url2 = _eutils_url("elink.fcgi", dbfrom="biosample", db="bioproject",
                       id=sample_batch)

    try:
        req = urllib.request.Request(url2, headers={"User-Agent": "PhenotypeCuration/0.1"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            xml_text = resp.read().decode("utf-8")
    except Exception:
        return []

    root2 = ET.fromstring(xml_text)
    for link_set in root2.iter("LinkSetDb"):
        db_to = link_set.findtext("DbTo", "")
        if "bioproject" in db_to.lower():
            for link in link_set.iter("Id"):
                if link.text:
                    bioprojects.add(link.text.strip())

    return sorted(bioprojects)


def _elink_via_sra(pmid: str, timeout: int) -> list[str]:
    """Via SRA: pubmed → sra → bioproject."""
    bioprojects: set[str] = set()

    # Step 1: pubmed → sra
    url1 = _eutils_url("elink.fcgi", dbfrom="pubmed", db="sra", id=pmid)

    try:
        req = urllib.request.Request(url1, headers={"User-Agent": "PhenotypeCuration/0.1"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            xml_text = resp.read().decode("utf-8")
    except Exception:
        return []

    root = ET.fromstring(xml_text)
    sra_ids: list[str] = []

    for link_set in root.iter("LinkSetDb"):
        for link in link_set.iter("Id"):
            if link.text:
                sra_ids.append(link.text.strip())

    if not sra_ids:
        return []

    # Step 2: sra → bioproject (SRA records link to their BioProject)
    sra_batch = ",".join(sra_ids[:100])
    url2 = _eutils_url("elink.fcgi", dbfrom="sra", db="bioproject", id=sra_batch)

    try:
        req = urllib.request.Request(url2, headers={"User-Agent": "PhenotypeCuration/0.1"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            xml_text = resp.read().decode("utf-8")
    except Exception:
        return []

    root2 = ET.fromstring(xml_text)
    for link_set in root2.iter("LinkSetDb"):
        db_to = link_set.findtext("DbTo", "")
        if "bioproject" in db_to.lower():
            for link in link_set.iter("Id"):
                if link.text:
                    bioprojects.add(link.text.strip())

    return sorted(bioprojects)


def _search_ena_by_pmid(pmid: str, timeout: int = 15) -> list[str]:
    """Search PubMed abstract and PMC full text for BioProject accessions.

    Strategy 4: Many papers list their data accession codes (PRJNA/PRJEB/PRJDB)
    in the abstract, methods, or data availability statement. This extracts them
    via regex from:
      1. PubMed abstract (efetch)
      2. PMC full-text XML (efetch), if the paper is in PMC

    This is a heuristic — it may find accessions for related but different
    projects. Results should be manually verified.
    """
    import re
    import os

    bioprojects: set[str] = set()
    api_key = os.environ.get("NCBI_API_KEY", "")
    api_suffix = f"&api_key={api_key}" if api_key else ""

    # Pattern: PRJ + (E|N|D) + (A|B|E|_)? + 4-10 digits
    prj_pattern = re.compile(r'(PRJ[END][A-Z]?\d{4,})', re.IGNORECASE)

    # Step 1: Search PubMed abstract
    abstract_url = (
        f"{NCBI_EUTILS_BASE}/efetch.fcgi"
        f"?db=pubmed&id={pmid}&rettype=abstract&retmode=text{api_suffix}"
    )
    try:
        req = urllib.request.Request(abstract_url, headers={"User-Agent": "PhenotypeCuration/0.1"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            text = resp.read().decode("utf-8", errors="replace")
            bioprojects.update(prj_pattern.findall(text))
    except Exception:
        pass

    # Step 2: Try PMC full text (if available)
    pmcid = _resolve_pmcid(pmid)
    if pmcid:
        pmc_url = (
            f"{NCBI_EUTILS_BASE}/efetch.fcgi"
            f"?db=pmc&id={pmcid}&rettype=xml{api_suffix}"
        )
        try:
            req = urllib.request.Request(pmc_url, headers={"User-Agent": "PhenotypeCuration/0.1"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                xml_text = resp.read().decode("utf-8", errors="replace")
                bioprojects.update(prj_pattern.findall(xml_text))
        except Exception:
            pass

    return sorted(bioprojects)


def _resolve_pmcid(pmid: str) -> str | None:
    """Resolve PMCID from PMID using NCBI elink.

    Returns PMCID (e.g., 'PMC7386241') or None.
    """
    url = _eutils_url("elink.fcgi", dbfrom="pubmed", db="pmc", id=pmid)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PhenotypeCuration/0.1"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            xml_text = resp.read().decode("utf-8")
        root = ET.fromstring(xml_text)
        # Without linkname, first LinkSetDb has the direct PMID→PMCID mapping
        for link_set in root.iter("LinkSetDb"):
            db_to = link_set.findtext("DbTo", "")
            if "pmc" in db_to.lower():
                for link in link_set.iter("Id"):
                    if link.text:
                        pmc_id = link.text.strip()
                        # Return the first non-reference ID (shorter, direct mapping)
                        if not pmc_id.startswith("PMC"):
                            return f"PMC{pmc_id}"
                        return pmc_id
    except Exception:
        pass
    return None


def batch_pmid_to_bioprojects(pmids: list[str], timeout: int = 15) -> dict[str, list[str]]:
    """Batch lookup: multiple PMIDs → their BioProjects.

    Args:
        pmids: List of PubMed IDs.
        timeout: HTTP timeout per call.

    Returns:
        Dict of PMID → list of bioproject accessions.
    """
    import time

    results: dict[str, list[str]] = {}
    for i, pmid in enumerate(pmids):
        results[pmid] = pmid_to_bioprojects(pmid, timeout)
        if (i + 1) % 3 == 0:
            time.sleep(0.35)  # NCBI rate limit

    return results


def link_paper_to_data(pmid: str) -> dict:
    """Full paper→data resolution for one PMID.

    Returns a dict with all discovered identifiers.

    Args:
        pmid: PubMed ID.

    Returns:
        Dict with keys: pmid, pmcid, bioprojects, title (if found).
    """
    result = {
        "pmid": pmid,
        "pmcid": None,
        "bioprojects": [],
        "title": None,
        "error": None,
    }

    # Find PMCID
    result["pmcid"] = _resolve_pmcid(pmid)

    # Find bioprojects
    result["bioprojects"] = pmid_to_bioprojects(pmid)
    return result


# ============================================================
# CLI
# ============================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="PMID → BioProject linker"
    )
    parser.add_argument("pmid", help="PubMed ID (single or comma-separated)")
    parser.add_argument("--batch", action="store_true",
                       help="Treat input as comma-separated list of PMIDs")
    parser.add_argument("--json", action="store_true", help="Output JSON")

    args = parser.parse_args()

    if args.batch:
        pmids = [p.strip() for p in args.pmid.split(",") if p.strip()]
        results = batch_pmid_to_bioprojects(pmids)
    else:
        # Single PMID, full resolution
        results = {"pmid": args.pmid}
        results["bioprojects"] = pmid_to_bioprojects(args.pmid)

    if args.json:
        import json
        print(json.dumps(results, indent=2))
    else:
        if isinstance(results, dict) and "pmid" in results:
            bp = results.get("bioprojects", [])
            if bp:
                print(f"PMID {results['pmid']} → BioProjects: {', '.join(bp)}")
            else:
                print(f"PMID {results['pmid']} → No BioProjects found")
        else:
            for pmid, bps in results.items():
                if bps:
                    print(f"PMID {pmid} → {', '.join(bps)}")
                else:
                    print(f"PMID {pmid} → (no bioprojects)")