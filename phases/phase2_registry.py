# -*- coding: utf-8 -*-
"""Phase 2: Auto-generate project registry from ENA metadata.

Input:  ena_metadata_master.csv (Phase 1 output)
Output: project_registry.csv (auto-generated, ready for manual review)

Scans ENA metadata per bioproject and auto-populates:
- Row counts (runs, biosamples, unique sample names)
- Metadata signal detection (intervention, timepoint, disease)
- Study design candidate classification
- Training readiness candidate
- Available host fields and common attribute keys

This is deterministic — no paper reading required.
Paper verification happens in Phase 3.
"""

import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional

_PARENT = Path(__file__).resolve().parent.parent
if str(_PARENT) not in sys.path:
    sys.path.insert(0, str(_PARENT))

from core.normalizer import (
    detect_available_fields,
    detect_subject_id_fields,
    detect_timepoint_fields,
    detect_host_fields,
)
from core.types import ProjectRegistryRow


# ============================================================
# Constants
# ============================================================

# ENA filereport columns (Strategy 1)
ENA_COLUMNS = [
    "run_accession", "sample_accession", "experiment_accession", "study_accession",
    "sample_alias", "sample_title", "sample_description",
    "library_name", "library_strategy", "library_source", "library_selection",
    "library_layout", "instrument_platform", "instrument_model",
    "fastq_ftp", "fastq_bytes", "fastq_md5",
    "host_organism", "host_sex", "host_body_site", "host_status",
    "host_phenotype", "host_genotype", "host_scientific_name", "host_tax_id",
    "collection_date", "collection_date_start", "collection_date_end",
    "geo_country", "geo_location", "scientific_name", "tax_id",
    "ena_first_public", "ena_last_updated", "read_count", "base_count",
    "bioproject",
]

# ENA structured host fields (column-based, not XML attributes)
STRUCTURED_HOST_FIELDS = [
    "host_sex", "host_body_site", "host_status",
    "host_phenotype", "host_genotype",
]

# Patterns that suggest intervention in sample aliases
INTERVENTION_ALIAS_PATTERNS = [
    r'\b(FMT|placebo|treatment|treated|post[-_]?treatment)\b',
    r'\b(dose|dosing|mg/kg|mg_?per_?kg)\b',
    r'\b(anti[-_]?PD[-_]?1|anti[-_]?CTLA|nivolumab|pembrolizumab|ipilimumab)\b',
    r'\b(prebiotic|probiotic|synbiotic|inulin|FOS|GOS|XOS)\b',
    r'\b(diet|feeding|fed|meal|intervention|supplement(?:ed)?)\b',
    r'\b(group|cohort|arm)[-_\s]*[A-Z0-9]',
    r'\b(responder|non[-_]?responder|R|NR)\b',
    r'\b(baseline|pre|post|before|after|follow[-_]?up)\b',
    r'\b(run[-_]?in|washout|wash[-_]?out)\b',
]

# Patterns that suggest longitudinal/multiple timepoints
TIMEPOINT_ALIAS_PATTERNS = [
    r'[TtWwVv]\d+',           # T0, T1, W0, W2, V1, v2
    r'[._-](day|d|week|wk|month|m|timepoint|tp)\d*[._-]',
    r'\b(baseline|follow[-_]?up|run[-_]?in)\b',
    r'(pre|post)[-_]',
    r'screening',
]

# Keywords in sample descriptions/titles that strongly suggest intervention studies
INTERVENTION_TEXT_KEYWORDS = [
    "anti-PD-1", "anti-PD1", "anti-CTLA", "anti-CTLA4",
    "nivolumab", "pembrolizumab", "ipilimumab", "atezolizumab",
    "immunotherapy", "immune checkpoint", "checkpoint inhibitor",
    "FMT", "fecal microbiota transplant", "faecal microbiota transplant",
    "probiotic", "prebiotic", "synbiotic", "inulin", "FOS", "GOS",
    "placebo", "randomized", "randomised", "treatment arm",
    "intervention", "dose escalation", "clinical trial",
    "responder", "non-responder", "progression",
    "RCT", "randomized controlled",
    "PD-1", "PD-L1", "CTLA-4",
    "vaginal seeding", "vaginal microbiome transfer",
    "RCC", "NSCLC", "melanoma",  # cancer types often studied in immunotherapy context
]

# Sample types that suggest gut metagenome (trainable)
GUT_SAMPLE_TYPES = [
    "feces", "stool", "faeces", "fecal", "faecal",
    "gut", "intestinal", "colon", "rectal",
    "human gut metagenome", "gut metagenome",
]


# ============================================================
# Scanners
# ============================================================

def _count_intervention_signals(aliases: list[str]) -> int:
    """Count how many sample aliases match intervention patterns."""
    count = 0
    for alias in aliases:
        for pattern in INTERVENTION_ALIAS_PATTERNS:
            if re.search(pattern, alias, re.IGNORECASE):
                count += 1
                break
    return count


def _count_timepoint_signals(aliases: list[str]) -> int:
    """Count how many sample aliases match timepoint patterns."""
    count = 0
    for alias in aliases:
        for pattern in TIMEPOINT_ALIAS_PATTERNS:
            if re.search(pattern, alias, re.IGNORECASE):
                count += 1
                break
    return count


def _detect_intervention_from_text(texts: list[str]) -> tuple[bool, list[str]]:
    """Check text fields (descriptions, titles) for intervention keywords.

    Returns:
        (has_intervention_signal, list of matched keywords)
    """
    matched = []
    for text in texts:
        if not text:
            continue
        text_lower = text.lower()
        for kw in INTERVENTION_TEXT_KEYWORDS:
            if kw.lower() in text_lower:
                if kw not in matched:
                    matched.append(kw)
    return bool(matched), matched


def _classify_study_design(
    has_intervention_meta: bool,
    has_timepoint_meta: bool,
    alias_intervention_ratio: float,
    alias_timepoint_ratio: float,
    raw_groups: str,
    has_text_intervention: bool = False,
) -> str:
    """Auto-classify study design from metadata signals.

    Returns one of: 'intervention_or_longitudinal_candidate',
    'observational_or_unclear', ''
    """
    # Strong signal: text descriptions mention intervention keywords
    if has_text_intervention:
        return "intervention_or_longitudinal_candidate"

    # Strong signal: explicit intervention group metadata
    if has_intervention_meta and raw_groups:
        return "intervention_or_longitudinal_candidate"

    # Strong signal: high ratio of intervention-like sample names
    if alias_intervention_ratio > 0.3:
        return "intervention_or_longitudinal_candidate"

    # Moderate signal: has timepoint metadata + some intervention pattern
    if has_timepoint_meta and alias_intervention_ratio > 0.1:
        return "intervention_or_longitudinal_candidate"

    # Longitudinally-sampled but unclear if intervention
    if alias_timepoint_ratio > 0.5 and has_timepoint_meta:
        return "intervention_or_longitudinal_candidate"

    # Default: no clear intervention signal
    return "observational_or_unclear"


def _classify_training_readiness(
    sample_titles: list[str],
    sample_descriptions: list[str],
    body_sites: list[str],
) -> str:
    """Classify training readiness based on sample type.

    Returns 'conditionally_ready' for gut/fecal samples,
    'pretraining_only_candidate' otherwise.
    """
    # Check sample titles and descriptions
    all_text = " ".join(
        str(t) for t in sample_titles + sample_descriptions + body_sites
        if t
    ).lower()

    for gut_type in GUT_SAMPLE_TYPES:
        if gut_type in all_text:
            return "conditionally_ready"

    # Check structured host_body_site
    for site in body_sites:
        if site and any(gt in site.lower() for gt in GUT_SAMPLE_TYPES):
            return "conditionally_ready"

    return "pretraining_only_candidate"


# ============================================================
# Registry Builder
# ============================================================

def _load_sample_attrs(sample_attrs_dir: Path, bioproject: str) -> list[dict]:
    """Load sample_attrs CSV for a project if it exists."""
    paths = list(sample_attrs_dir.glob(f"{bioproject}*sample_attrs*"))
    if not paths:
        return []
    rows = []
    for path in paths:
        try:
            with open(path, "r", encoding="utf-8-sig") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    rows.append(row)
        except Exception:
            continue
    return rows


def _parse_all_attributes_json(json_str: str) -> dict:
    """Parse all_attributes_json field from sample_attrs CSV."""
    if not json_str:
        return {}
    try:
        return json.loads(json_str)
    except (json.JSONDecodeError, TypeError):
        return {}


def build_registry(
    metadata_path: str | Path,
    sample_attrs_dir: Optional[str | Path] = None,
    out_path: Optional[str | Path] = None,
) -> list[ProjectRegistryRow]:
    """Generate project registry from ENA metadata.

    Args:
        metadata_path: Path to ena_metadata_master.csv or directory of per-project files.
        sample_attrs_dir: Directory containing *_sample_attrs.csv files (Strategy 2 output).
        out_path: Optional output CSV path. If None, writes to project_registry.csv in cwd.

    Returns:
        List of ProjectRegistryRow dataclasses.
    """
    metadata_path = Path(metadata_path)
    sample_attrs_dir = Path(sample_attrs_dir) if sample_attrs_dir else None

    # --- Step 1: Load all metadata rows ---
    print(f"Loading metadata from {metadata_path}...")
    all_rows: list[dict] = []

    if metadata_path.is_dir():
        for f in sorted(metadata_path.glob("*.csv")):
            with open(f, "r", encoding="utf-8-sig") as fh:
                reader = csv.DictReader(fh)
                all_rows.extend(list(reader))
    else:
        with open(metadata_path, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            all_rows = list(reader)

    print(f"  Loaded {len(all_rows)} rows")

    # --- Step 2: Group by bioproject ---
    by_project: dict[str, list[dict]] = defaultdict(list)
    for row in all_rows:
        bp = row.get("bioproject", "").strip()
        if bp:
            by_project[bp].append(row)

    print(f"  Found {len(by_project)} bioprojects")

    # --- Step 3: Build registry rows ---
    registry_rows: list[ProjectRegistryRow] = []

    for bioproject in sorted(by_project.keys()):
        rows = by_project[bioproject]
        print(f"\n  {bioproject}: {len(rows)} runs")

        # Counts
        run_count = len(rows)
        biosamples = set()
        sample_names = set()
        sample_titles = []
        sample_descriptions = []
        sample_aliases = []
        body_sites = []
        host_sexes = set()

        for row in rows:
            bs = row.get("sample_accession", "").strip()
            if bs:
                biosamples.add(bs)
            sn = row.get("sample_alias", "").strip()
            if sn:
                sample_names.add(sn)
                sample_aliases.append(sn)
            st = row.get("sample_title", "").strip()
            if st:
                sample_titles.append(st)
            sd = row.get("sample_description", "").strip()
            if sd:
                sample_descriptions.append(sd)
            hs = row.get("host_sex", "").strip()
            if hs:
                host_sexes.add(hs)
            hbs = row.get("host_body_site", "").strip()
            if hbs:
                body_sites.append(hbs)

        biosample_count = len(biosamples)
        sample_name_count = len(sample_names)

        # Metadata signal checks from ENA structured columns
        has_intervention_meta = bool(host_sexes)  # ENA structured host fields as signal
        has_timepoint_structured = any(
            row.get("collection_date", "").strip() for row in rows
        )

        # --- Step 3a: Load sample attributes (Strategy 2) if available ---
        sample_attrs_rows = []
        all_raw_attrs: dict[str, int] = Counter()  # key → count of samples with that key
        raw_group_values_set = set()
        raw_intervention_values_set = set()
        raw_disease_values_set = set()

        if sample_attrs_dir and sample_attrs_dir.exists():
            sample_attrs_rows = _load_sample_attrs(sample_attrs_dir, bioproject)

        if sample_attrs_rows:
            print(f"    Loaded {len(sample_attrs_rows)} sample attribute rows")

            # Collect all attribute keys
            for sa in sample_attrs_rows:
                raw = _parse_all_attributes_json(sa.get("all_attributes_json", ""))
                for k in raw:
                    all_raw_attrs[k] += 1

            # Detect available standard fields via normalizer
            # Use the first non-empty attributes dict
            first_raw = {}
            for sa in sample_attrs_rows:
                raw = _parse_all_attributes_json(sa.get("all_attributes_json", ""))
                if raw:
                    first_raw = raw
                    break

            if first_raw:
                detected = detect_available_fields(first_raw)
                if "intervention_group" in detected:
                    has_intervention_meta = True
                if "timepoint" in detected:
                    has_timepoint_structured = True

            # Collect raw values from sample attrs
            for sa in sample_attrs_rows:
                group_val = sa.get("group_raw", "").strip()
                if group_val:
                    raw_group_values_set.add(group_val)
                disease_val = sa.get("disease_raw", "").strip()
                if disease_val:
                    raw_disease_values_set.add(disease_val)
        else:
            # Fallback: detect signals from sample aliases
            for alias in sample_aliases:
                for pattern in INTERVENTION_ALIAS_PATTERNS:
                    m = re.search(pattern, alias, re.IGNORECASE)
                    if m:
                        raw_group_values_set.add(m.group(0).lower())

        # --- Step 3b: Compute signal metrics ---
        alias_intervention_count = _count_intervention_signals(sample_aliases)
        alias_timepoint_count = _count_timepoint_signals(sample_aliases)
        total_aliases = len(sample_aliases) or 1

        alias_intervention_ratio = alias_intervention_count / total_aliases
        alias_timepoint_ratio = alias_timepoint_count / total_aliases

        raw_groups_str = "; ".join(sorted(raw_group_values_set)[:15])

        # --- Step 3c: Classify ---
        # Detect intervention signals from text descriptions/titles
        has_text_intv, text_intv_keywords = _detect_intervention_from_text(
            sample_descriptions + sample_titles
        )

        study_design = _classify_study_design(
            has_intervention_meta=has_intervention_meta,
            has_timepoint_meta=has_timepoint_structured,
            alias_intervention_ratio=alias_intervention_ratio,
            alias_timepoint_ratio=alias_timepoint_ratio,
            raw_groups=raw_groups_str,
            has_text_intervention=has_text_intv,
        )

        training_readiness = _classify_training_readiness(
            sample_titles=sample_titles,
            sample_descriptions=sample_descriptions,
            body_sites=body_sites,
        )

        # --- Step 3d: Detect available fields for registry ---
        common_subject_id_fields = ""
        common_timepoint_fields = ""
        available_host_fields = ""

        if sample_attrs_rows:
            first_raw = {}
            for sa in sample_attrs_rows:
                raw = _parse_all_attributes_json(sa.get("all_attributes_json", ""))
                if raw:
                    first_raw = raw
                    break

            if first_raw:
                subject_fields = detect_subject_id_fields(first_raw)
                timepoint_fields = detect_timepoint_fields(first_raw)
                host_fields = detect_host_fields(first_raw)

                common_subject_id_fields = ", ".join(subject_fields[:5])
                common_timepoint_fields = ", ".join(timepoint_fields[:5])
                available_host_fields = ", ".join(host_fields)

        # --- Step 3e: Build registry row ---
        registry_row = ProjectRegistryRow(
            bioproject=bioproject,
            keep_flag="TRUE",
            scope_basis="auto_generated_from_ena",
            primary_pmid="",    # To be filled manually or from PMID source
            additional_pmids="",
            paper_count="",
            run_n_local=str(run_count),
            biosample_n_local=str(biosample_count),
            sample_name_n_local=str(sample_name_count),
            has_intervention_metadata=str(has_intervention_meta),
            has_timepoint_metadata=str(has_timepoint_structured),
            raw_group_values=raw_groups_str,
            raw_intervention_values="; ".join(sorted(raw_intervention_values_set)[:10]),
            raw_disease_values="; ".join(sorted(raw_disease_values_set)[:10]),
            study_design_candidate=study_design,
            training_readiness_candidate=training_readiness,
            manual_review_status="pending_paper_verification",
            supplement_status="to_be_verified",
            common_subject_id_fields=common_subject_id_fields,
            common_timepoint_fields=common_timepoint_fields,
            available_host_fields=available_host_fields,
        )
        registry_rows.append(registry_row)

        print(f"    runs={run_count}, biosamples={biosample_count}, "
              f"design={study_design}, readiness={training_readiness}")

    # --- Step 4: Write output ---
    if out_path is None:
        out_path = Path("project_registry.csv")

    out_path = Path(out_path)
    _write_registry_csv(registry_rows, out_path)

    return registry_rows


def _write_registry_csv(rows: list[ProjectRegistryRow], path: Path):
    """Write registry rows to CSV."""
    from dataclasses import asdict

    fieldnames = [
        "bioproject", "keep_flag", "scope_basis",
        "primary_pmid", "additional_pmids", "paper_count",
        "run_n_local", "biosample_n_local", "sample_name_n_local",
        "has_intervention_metadata", "has_timepoint_metadata",
        "raw_group_values", "raw_intervention_values", "raw_disease_values",
        "study_design_candidate", "training_readiness_candidate",
        "manual_review_status", "supplement_status",
        "common_subject_id_fields", "common_timepoint_fields",
        "available_host_fields",
    ]

    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            d = asdict(row)
            writer.writerow({k: d.get(k, "") for k in fieldnames})

    print(f"\nWrote {len(rows)} registry rows to {path}")


# ============================================================
# CLI
# ============================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Phase 2: Auto-generate project registry from ENA metadata"
    )
    parser.add_argument(
        "metadata",
        help="Path to ena_metadata_master.csv or directory of per-project CSVs",
    )
    parser.add_argument(
        "--sample-attrs-dir",
        help="Directory containing *_sample_attrs.csv files (from Strategy 2 download)",
        default=None,
    )
    parser.add_argument(
        "--out",
        default="project_registry.csv",
        help="Output CSV path (default: project_registry.csv)",
    )

    args = parser.parse_args()
    build_registry(
        metadata_path=args.metadata,
        sample_attrs_dir=args.sample_attrs_dir,
        out_path=args.out,
    )