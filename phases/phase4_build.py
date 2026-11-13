# -*- coding: utf-8 -*-
"""Phase 4: Build phenotype table from ENA metadata + intervention registry.

Input:
  - ena_metadata_master.csv (Phase 1)
  - intervention_registry.csv (Phase 3, paper-verified)
  - *_sample_attrs.csv (Phase 1 Strategy 2, optional but recommended)

Output:
  - phenotype_clean.csv (36 columns per PhenotypeRow schema)

Architecture:
  Registry-driven. The intervention_registry.csv is the single source of truth
  for each project's study design, intervention type, disease label, and
  timepoint phase mapping rules. Per-project custom logic is limited to
  subject ID parsing and timepoint value extraction — everything else is
  determined by registry fields.
"""

import csv
import json
import re
import sys
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import Optional

_PARENT = Path(__file__).resolve().parent.parent
if str(_PARENT) not in sys.path:
    sys.path.insert(0, str(_PARENT))

from core.types import (
    PhenotypeRow, InterventionRegistryRow, StudyDesign, InterventionType,
    IsIntervention,
)
from core.normalizer import normalize_value, get_field_value


# ============================================================
# Registry Loading
# ============================================================

def load_intervention_registry(path: str | Path) -> dict[str, InterventionRegistryRow]:
    """Load Phase 3 intervention registry into a dict keyed by bioproject."""
    path = Path(path)
    registry: dict[str, InterventionRegistryRow] = {}

    with open(path, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            bp = row.get("bioproject", "").strip()
            if not bp:
                continue
            registry[bp] = InterventionRegistryRow(
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

    return registry


# ============================================================
# Metadata Loading
# ============================================================

def load_ena_metadata(metadata_path: str | Path) -> dict[str, list[dict]]:
    """Load ENA metadata master CSV and group by bioproject."""
    metadata_path = Path(metadata_path)
    by_project: dict[str, list[dict]] = defaultdict(list)

    load_paths = []
    if metadata_path.is_dir():
        load_paths = sorted(metadata_path.glob("*.csv"))
    else:
        load_paths = [metadata_path]

    for lp in load_paths:
        with open(lp, "r", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                bp = row.get("bioproject", "").strip()
                if bp:
                    by_project[bp].append(row)

    return dict(by_project)


def load_sample_attrs(sample_attrs_dir: Path, bioproject: str) -> dict[str, dict]:
    """Load per-sample attributes from Strategy 2 download.

    Returns dict of sample_accession → normalized attributes dict.
    """
    attrs: dict[str, dict] = {}
    paths = list(sample_attrs_dir.glob(f"{bioproject}*sample_attrs*"))
    if not paths:
        return attrs

    for path in paths:
        try:
            with open(path, "r", encoding="utf-8-sig") as f:
                for row in csv.DictReader(f):
                    acc = row.get("sample_accession", "").strip()
                    if acc:
                        attrs[acc] = dict(row)
        except Exception:
            continue

    return attrs


# ============================================================
# Subject ID Parsing
# ============================================================

def parse_subject_id(
    sample_alias: str,
    sample_attrs: Optional[dict] = None,
    bioproject: str = "",
) -> tuple[str, str, str]:
    """Extract subject ID from sample alias and attributes.

    Returns:
        (subject_id, source_field, confidence)
    """
    alias = sample_alias.strip() if sample_alias else ""

    # Try sample attributes first
    if sample_attrs:
        sid = sample_attrs.get("subject_id_raw", "").strip()
        if sid and sid.lower() not in ("unknown", "na", "not collected", ""):
            return sid, "ena_host_subject_id", "parsed"

    # Parse from sample alias by common patterns
    if not alias:
        return alias, "unavailable", "cannot_parse"

    # Pattern 1: PRJEB75958 — "host subject id:P01;Sample name:R_P01_W00"
    m = re.search(r"host subject id:([^;]+)", alias)
    if m:
        return m.group(1).strip(), "sample_name_field", "parsed"

    # Pattern 2: PRJEB22863 — "NSCLC_B_001_T0" → "001"
    m = re.search(r"_(\d{3,4})_", alias)
    if m:
        return m.group(1).strip(), "sample_name_field", "parsed"

    # Pattern 3: PRJNA737472 — "host_assigned_donor_subject_id:D0044;host_subject_id:S0060;sample_ID:SS01218"
    m = re.search(r"host_subject_id:([^;]+)", alias)
    if m:
        return m.group(1).strip(), "sample_name_field", "parsed"

    # Pattern 4: PRJNA701480 — "baby_ID:GK9;host_subject_id:M06;mum_ID:M06"
    m = re.search(r"host_subject_id:([^;]+)", alias)
    if m:
        return m.group(1).strip(), "sample_name_field", "parsed"

    # Pattern 5: General — first segment before underscore if looks like ID
    parts = alias.split("_")
    if parts:
        # If first part looks like a subject ID (alphanumeric with digits)
        first = parts[0].strip()
        if re.search(r"\d", first) and 2 <= len(first) <= 12:
            return first, "sample_name_field", "parsed_numeric_only"
        # Short alpha-only is low confidence
        if len(first) <= 3 and first.isalpha():
            return first, "ena_host_subject_id", "low_quality_short_code"

    return alias, "ena_sample_name", "parsed"


# ============================================================
# Timepoint Extraction
# ============================================================

def extract_timepoint_index(sample_alias: str, sample_attrs: Optional[dict] = None) -> str:
    """Extract numeric timepoint index from sample alias."""
    alias = sample_alias.strip() if sample_alias else ""

    # From attributes
    if sample_attrs:
        tp = sample_attrs.get("timepoint_raw", "").strip()
        if tp:
            return tp

    # Pattern: T0, T1, T2, W0, W2, V1, v2, d0, day0
    m = re.search(r"[TWVDtwvd](\d+)", alias)
    if m:
        return str(int(m.group(1)))

    # Pattern: _number at end
    m = re.search(r"_(\d+)$", alias)
    if m:
        return str(int(m.group(1)))

    # Pattern: timepoint, tp, visit, week
    m = re.search(r"(?:timepoint|tp|visit|week|day)[_\s]*(\d+)", alias, re.IGNORECASE)
    if m:
        return str(int(m.group(1)))

    # Pattern: N_#### (subject N, sample ####)
    m = re.search(r"(\d{4,6})$", alias)
    if m:
        return str(int(m.group(1)))

    return ""


def extract_timepoint_label(sample_alias: str) -> str:
    """Extract descriptive timepoint label from sample alias."""
    alias = sample_alias.strip() if sample_alias else ""

    label_map = {
        "baseline": "baseline",
        "run_in": "run_in",
        "runin": "run_in",
        "pre": "baseline",
        "post": "post",
        "followup": "followup",
        "follow_up": "followup",
        "screening": "screening",
        "screen": "screening",
        "washout": "washout",
        "wash_out": "washout",
    }

    alias_lower = alias.lower().replace("_", " ").replace("-", " ")
    for key, label in label_map.items():
        if key.replace("_", " ") in alias_lower:
            return label

    # T0/W0/day0 → baseline
    if re.search(r"[TWVDtwvd]0\b", alias):
        return "baseline"

    # On-treatment signals
    if re.search(r"on[_\s]?treat|on[_\s]?therapy|during", alias_lower):
        return "on_treatment"

    # Post-treatment
    if re.search(r"post[_\s]?treat|post[_\s]?therapy|after", alias_lower):
        return "post_treatment"

    return "longitudinal_unspecified"


# ============================================================
# is_intervention Decision
# ============================================================

def determine_is_intervention(
    bioproject: str,
    timepoint_label: str,
    timepoint_index: str,
    registry: InterventionRegistryRow,
    phase_map: Optional[list] = None,
) -> str:
    """Determine whether a sample is from an active intervention period.

    Decision tree:
    1. Observational studies → always 'no'
    2. One-time interventions (single dose/procedure) → always 'no' (all samples are pre/post)
    3. Ongoing interventions → timepoint-dependent
       a. Use explicit phase_map if provided (e.g. FOS/PDX dose-escalation)
       b. Use registry sample_timing_schedule / days_since_intervention
       c. Fallback heuristic based on timepoint_label

    Args:
        bioproject: Project accession.
        timepoint_label: Normalized phase label.
        timepoint_index: Numeric index (as string).
        registry: InterventionRegistryRow for this project.
        phase_map: Optional list of (min_day, max_day, phase_label, is_intervention) tuples.

    Returns:
        'yes', 'no', or 'unknown'
    """
    # Rule 1: Observational → always no
    if registry.intervention_type == InterventionType.OBSERVATIONAL.value:
        return "no"

    # Rule 2: One-time intervention → always no
    if registry.is_one_time == "yes":
        return "no"

    # Rule 3: Use explicit phase map
    if phase_map:
        try:
            day = float(timepoint_index) if timepoint_index else 0
        except ValueError:
            day = 0

        for min_day, max_day, phase_label, is_intv in phase_map:
            if min_day <= day <= max_day:
                return is_intv

        return "unknown"

    # Rule 4: Heuristic based on timepoint label
    label = timepoint_label.lower()
    if label in ("baseline", "screening", "run_in", "pre"):
        return "no"
    if label in ("washout", "followup", "post_treatment", "post"):
        return "no"
    if label in ("feeding_1", "feeding_2", "on_treatment", "during_diet"):
        return "yes"
    if label in ("longitudinal_unspecified",):
        # Can't determine without more info → unknown
        return "unknown"

    return "unknown"


# ============================================================
# Phase Map Definitions (per-project, paper-verified)
# ============================================================

# Format: list of (min_day, max_day, phase_label, is_intervention)
# These are the phase maps verified from paper supplements.
# New projects add entries here during Phase 3 (paper verification).

PHASE_MAPS: dict[str, list] = {
    # PMC7386241: FOS/PDX dose-escalation (Creswell 2020)
    "PRJNA594393": [  # K001 FOS
        (float("-inf"), 0, "run_in", "no"),
        (1, 14, "feeding_1", "yes"),
        (15, 28, "feeding_2", "yes"),
        (29, float("inf"), "washout", "no"),
    ],
    "PRJNA594610": [  # K002 PDX
        (float("-inf"), 0, "run_in", "no"),
        (1, 14, "feeding_1", "yes"),
        (15, 28, "feeding_2", "yes"),
        (29, float("inf"), "washout", "no"),
    ],
    "PRJNA594613": [  # K003 FOS→PDX
        (float("-inf"), 0, "run_in", "no"),
        (1, 14, "feeding_1", "yes"),
        (15, 42, "washout_1", "no"),
        (43, 56, "feeding_2", "yes"),
        (57, float("inf"), "washout_2", "no"),
    ],
    "PRJNA594620": [  # K003 PDX→FOS
        (float("-inf"), 0, "run_in", "no"),
        (1, 14, "feeding_1", "yes"),
        (15, 42, "washout_1", "no"),
        (43, 56, "feeding_2", "yes"),
        (57, float("inf"), "washout_2", "no"),
    ],
}


# ============================================================
# Per-Project Timepoint Overrides
# ============================================================

# Maps timepoint_index (str) → (phase_label, is_intervention) for projects
# where the phase map by day doesn't apply (non-numeric timepoints).

TIMEPOINT_OVERRIDES: dict[str, dict] = {
    "PRJEB22863": {  # NSCLC anti-PD-1: 0=T0, 1=T1, 2=T2
        "0": ("baseline", "no"),
        "1": ("on_treatment", "yes"),
        "2": ("on_treatment", "yes"),
    },
    "PRJEB33500": {  # Mediterranean diet crossover
        "baseline": ("baseline", "no"),
        "longitudinal_unspecified": ("during_diet", "yes"),
    },
    "PRJNA541981": {  # Melanoma immunotherapy
        "baseline": ("baseline", "no"),
        "longitudinal_unspecified": ("on_treatment", "yes"),
    },
    "PRJNA637785": {  # Gut Bugs Trial single-dose FMT
        "baseline": ("baseline", "no"),
        "longitudinal_unspecified": ("post_fmt_followup", "no"),
    },
    "PRJNA701480": {  # Vaginal seeding at birth
        "longitudinal_unspecified": ("post_birth_followup", "no"),
    },
    "PRJNA762360": {  # Melanoma anti-PD-1 cross-sectional
        "longitudinal_unspecified": ("mixed_timing", "unknown"),
    },
    "PRJNA904790": {  # FMT multi-dose
        "pre": ("baseline", "no"),
        "longitudinal_unspecified": ("during_fmt", "yes"),
    },
    "PRJEB75958": {  # FMT capsules weekly x 8 weeks
        "0": ("baseline", "no"),
        "1": ("during_fmt", "yes"),
        "2": ("post_fmt_followup", "no"),
    },
}


# ============================================================
# Per-Project Disease Map
# ============================================================

# Paper-verified disease labels (from registry + manual curation).
# Populated during Phase 3. Defaults to "" (unknown).

DISEASE_MAP: dict[str, str] = {
    "PRJEB22863": "NSCLC",
    "PRJEB33500": "healthy",
    "PRJEB39610": "preterm_birth",
    "PRJEB48663": "healthy_infant",
    "PRJEB75958": "HIV",
    "PRJNA397450": "healthy",
    "PRJNA541981": "melanoma",
    "PRJNA594393": "healthy",
    "PRJNA594610": "healthy",
    "PRJNA594613": "healthy",
    "PRJNA594620": "healthy",
    "PRJNA637785": "obesity",
    "PRJNA701480": "healthy_infant",
    "PRJNA737472": "ulcerative_colitis",
    "PRJNA762360": "melanoma",
    "PRJNA895487": "autism_spectrum_disorder",
    "PRJNA904790": "healthy",
    "PRJNA929328": "hematologic_malignancy",
    "PRJNA929986": "healthy_child",
}

# Intervention agent descriptions
INTERVENTION_AGENT_MAP: dict[str, str] = {
    "PRJEB22863": "anti-PD-1 (nivolumab/pembrolizumab)",
    "PRJEB33500": "Mediterranean diet",
    "PRJEB75958": "FMT capsules (healthy donor stool)",
    "PRJNA541981": "anti-PD-1 / anti-CTLA-4",
    "PRJNA594393": "FOS (fructo-oligosaccharides)",
    "PRJNA594610": "PDX (polydextrose)",
    "PRJNA594613": "FOS + PDX crossover (FOS first)",
    "PRJNA594620": "PDX + FOS crossover (PDX first)",
    "PRJNA637785": "FMT capsules (lean donor stool)",
    "PRJNA701480": "vaginal seeding (maternal vaginal microbiota)",
    "PRJNA737472": "FMT capsules / ENMA",
    "PRJNA762360": "anti-PD-1",
    "PRJNA904790": "FMT (fecal microbiota transplantation)",
    "PRJNA929986": "BB-12 (Bifidobacterium animalis subsp. lactis)",
}


# ============================================================
# Main Builder
# ============================================================

def build_phenotype_table(
    metadata_path: str | Path,
    registry_path: str | Path,
    sample_attrs_dir: Optional[str | Path] = None,
    out_path: Optional[str | Path] = None,
) -> list[PhenotypeRow]:
    """Build phenotype_clean.csv from ENA metadata + intervention registry.

    Args:
        metadata_path: Path to ena_metadata_master.csv (Phase 1).
        registry_path: Path to intervention_registry.csv (Phase 3).
        sample_attrs_dir: Directory with *_sample_attrs.csv (Phase 1 Strategy 2).
        out_path: Output CSV path (default: phenotype_clean.csv).

    Returns:
        List of PhenotypeRow dataclasses.
    """
    # --- Load inputs ---
    print(f"Loading registry from {registry_path}...")
    registry = load_intervention_registry(registry_path)
    print(f"  Loaded {len(registry)} project registrations")

    print(f"Loading metadata from {metadata_path}...")
    by_project = load_ena_metadata(metadata_path)
    total_runs = sum(len(v) for v in by_project.values())
    print(f"  Loaded {total_runs} runs across {len(by_project)} projects")

    sample_attrs = {}
    if sample_attrs_dir:
        sample_attrs_dir = Path(sample_attrs_dir)

    # --- Build rows ---
    phenotype_rows: list[PhenotypeRow] = []
    project_count = 0
    skip_count = 0

    for bioproject in sorted(by_project.keys()):
        runs = by_project[bioproject]

        # Load sample attributes
        sa_dict: dict[str, dict] = {}
        if sample_attrs_dir and sample_attrs_dir.exists():
            sa_dict = load_sample_attrs(sample_attrs_dir, bioproject)

        # Get registry entry
        reg = registry.get(bioproject)
        if reg is None:
            print(f"  SKIP {bioproject}: not in registry ({len(runs)} runs)")
            skip_count += len(runs)
            continue

        project_count += 1
        print(f"  {bioproject}: {len(runs)} runs, "
              f"design={reg.study_design}, intv={reg.intervention_type}")

        # Get per-project config
        disease = DISEASE_MAP.get(bioproject, "")
        intervention_agent = INTERVENTION_AGENT_MAP.get(bioproject, "")
        phase_map = PHASE_MAPS.get(bioproject)
        tp_overrides = TIMEPOINT_OVERRIDES.get(bioproject, {})

        # --- Step 1: Group runs by biosample ---
        # One biosample can have multiple runs; phenotype is per-biosample
        by_biosample: dict[str, list[dict]] = defaultdict(list)
        for run in runs:
            sample_acc = run.get("sample_accession", "").strip()
            if sample_acc:
                by_biosample[sample_acc].append(run)

        print(f"    {len(by_biosample)} biosamples from {len(runs)} runs")

        # --- Step 2: Build one row per biosample ---
        sample_groups: dict[tuple, list[dict]] = defaultdict(list)

        for sample_acc, sample_runs in by_biosample.items():
            # Use the first run's metadata (same biosample shares most fields)
            run = sample_runs[0]
            run_acc = ";".join(
                sorted(set(r.get("run_accession", "").strip() for r in sample_runs
                           if r.get("run_accession", "").strip()))
            )
            sample_alias = run.get("sample_alias", "").strip()

            # Get sample attributes
            attrs = sa_dict.get(sample_acc, {})

            # Parse subject ID
            subject_id, sid_source, sid_confidence = parse_subject_id(
                sample_alias, attrs, bioproject
            )

            # Extract timepoint
            tp_index = extract_timepoint_index(sample_alias, attrs)
            tp_label = extract_timepoint_label(sample_alias)

            # Apply timepoint overrides
            is_intervention = "unknown"
            if tp_label in tp_overrides:
                tp_label, is_intervention = tp_overrides[tp_label]
            elif tp_index in tp_overrides:
                tp_label, is_intervention = tp_overrides[tp_index]
            else:
                # Use phase map or generic logic
                is_intervention = determine_is_intervention(
                    bioproject, tp_label, tp_index, reg, phase_map
                )

            # Timepoint unit inference
            tp_unit = "timepoint_index"  # default
            if re.search(r"day|d\b", sample_alias.lower()):
                tp_unit = "day"
            elif re.search(r"week|wk|w\b", sample_alias.lower()):
                tp_unit = "week"

            # Host metadata
            host_age = attrs.get("age_raw", "").strip() if attrs else ""
            host_sex = attrs.get("sex_raw", "").strip() if attrs else ""
            host_bmi = attrs.get("bmi_raw", "").strip() if attrs else ""

            # Fallback to ENA structured fields
            if not host_sex:
                host_sex = run.get("host_sex", "").strip()
                host_sex = normalize_value("host_sex", host_sex)

            # Intervention group from attributes
            intervention_group = attrs.get("group_raw", "").strip() if attrs else ""

            # Build phenotype row
            row = PhenotypeRow(
                bioproject=bioproject,
                biosample_id=sample_acc,
                run_accession=run_acc,
                subject_id=subject_id,
                sample_id_local=sample_alias,
                sample_type=run.get("library_source", ""),
                collection_date=run.get("collection_date", ""),
                library_layout=run.get("library_layout", ""),
                platform=run.get("instrument_platform", ""),
                disease=disease,
                intervention_type=reg.intervention_type,
                intervention_agent=intervention_agent,
                intervention_group=intervention_group,
                intervention_source="registry",
                host_age=host_age,
                host_sex=host_sex,
                host_bmi=host_bmi,
                ethnicity="",
                clinical_outcome="",
                phenotype_source="ena",
                timepoint_label=tp_label,
                timepoint_index=tp_index,
                timepoint_unit=tp_unit,
                is_intervention=is_intervention,
                paired_pre_sample="",
                paired_post_sample="",
                sample_seq="",
                total_samples="",
                pairing_confidence="unpaired",
                parallel_samples="",
                context_role="",
                context_id="",
                family_id="",
                notes="",
                subject_id_status="parsed_candidate",
                subject_id_confidence=sid_confidence,
            )
            phenotype_rows.append(row)

            # Track for pairing
            sample_groups[(bioproject, subject_id)].append(
                {"row": row, "tp_index": tp_index, "tp_label": tp_label}
            )

        # --- Pairing pass ---
        _pair_samples(sample_groups)

        print(f"    built {len(by_biosample)} phenotype rows from {len(runs)} runs")

    # --- Pairing pass (global) ---

    # --- Write output ---
    if out_path is None:
        out_path = Path("phenotype_clean.csv")
    out_path = Path(out_path)

    _write_phenotype_csv(phenotype_rows, out_path)

    print(f"\nWrote {len(phenotype_rows)} rows to {out_path}")
    if skip_count:
        print(f"  ({skip_count} runs skipped: projects not in registry)")

    return phenotype_rows


# ============================================================
# Sample Pairing
# ============================================================

def _pair_samples(sample_groups: dict[tuple, list[dict]]):
    """Compute longitudinal pair relationships within each subject."""
    for (bioproject, subject_id), samples in sample_groups.items():
        if len(samples) < 2:
            continue

        # Try to sort by numeric timepoint_index
        numeric = []
        for s in samples:
            try:
                numeric.append((float(s["tp_index"]), s))
            except (ValueError, TypeError):
                pass

        if len(numeric) < 2:
            # Try string sort
            sorted_samples = sorted(
                samples, key=lambda s: s["tp_index"] or ""
            )
            # Assign sequence numbers
            for i, s in enumerate(sorted_samples):
                s["row"].sample_seq = str(i + 1)
                s["row"].total_samples = str(len(samples))
                s["row"].pairing_confidence = "candidate"
            continue

        # Sort by numeric timepoint
        numeric.sort(key=lambda x: x[0])
        sorted_samples = [ns[1] for ns in numeric]

        for i, s in enumerate(sorted_samples):
            s["row"].sample_seq = str(i + 1)
            s["row"].total_samples = str(len(sorted_samples))
            s["row"].pairing_confidence = "candidate"

            if i > 0:
                s["row"].paired_pre_sample = sorted_samples[i - 1]["row"].biosample_id
            if i < len(sorted_samples) - 1:
                s["row"].paired_post_sample = sorted_samples[i + 1]["row"].biosample_id


# ============================================================
# CSV Output
# ============================================================

def _write_phenotype_csv(rows: list[PhenotypeRow], path: Path):
    """Write phenotype rows to CSV."""
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=PhenotypeRow.OUTPUT_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


# ============================================================
# CLI
# ============================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Phase 4: Build phenotype table from registry + metadata"
    )
    parser.add_argument(
        "metadata",
        help="Path to ena_metadata_master.csv",
    )
    parser.add_argument(
        "--registry",
        required=True,
        help="Path to intervention_registry.csv (Phase 3 output)",
    )
    parser.add_argument(
        "--sample-attrs-dir",
        help="Directory with *_sample_attrs.csv (Phase 1 Strategy 2)",
        default=None,
    )
    parser.add_argument(
        "--out",
        default="phenotype_clean.csv",
        help="Output CSV path",
    )

    args = parser.parse_args()
    build_phenotype_table(
        metadata_path=args.metadata,
        registry_path=args.registry,
        sample_attrs_dir=args.sample_attrs_dir,
        out_path=args.out,
    )