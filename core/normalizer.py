# -*- coding: utf-8 -*-
"""Field name normalization for ENA/NCBI metadata attributes.

Addresses Gap #3: ENA and NCBI BioSample use inconsistent field names
for the same biological concepts (age, sex, timepoint, etc.).

This module provides:
- A curated synonym table mapping raw attribute keys → standard phenotype fields
- Auto-detection of which fields are available in a project's metadata
- Value-level normalization (sex → male/female, age → numeric, etc.)
- Both ENA SAMPLE_ATTRIBUTE and NCBI Attribute element naming conventions

Design principle:
    The synonym table is the single source of truth. Each entry in SYNONYMS
    maps one standard field to a priority-ordered list of candidate raw keys.
    Keys are matched case-insensitively with underscores/spaces normalized.
"""

import re
from typing import Optional


# ============================================================
# Synonym Table
# ============================================================

# Each group maps a standard phenotype field name to a priority-ordered
# list of raw ENA/NCBI attribute keys. Keys are matched case-insensitively,
# with underscores, spaces, and hyphens treated as equivalent.
# The FIRST match in priority order wins.

SYNONYMS: dict[str, list[str]] = {
    # --- Identity ---
    "subject_id": [
        # Most-specific keys first
        "subject",
        "host_subject_id",
        "subject_id",
        "host_subject",
        "individual",
        "individual_id",
        "patient",
        "patient_id",
        "donor",
        "donor_id",
        "participant",
        "participant_id",
        "host",
        "host_id",
        "sample_name",          # fallback — often contains subject prefix
        "sample_title",         # fallback — sometimes structured
    ],

    # --- Demographics ---
    "host_age": [
        "host_age",
        "age",
        "age_years",
        "age_year",
        "age_at_collection",
        "age_collected",
        "host_age_years",
        "AGE",
        "Age",
    ],

    "host_sex": [
        "host_sex",
        "sex",
        "gender",
        "Sex",
        "Gender",
        "biological_sex",
        "host_gender",
        "Sex_gender",
    ],

    "host_bmi": [
        "host_body_mass_index",
        "body_mass_index",
        "bmi",
        "BMI",
        "host_bmi",
    ],

    "ethnicity": [
        "ethnicity",
        "host_ethnicity",
        "race",
        "host_race",
    ],

    # --- Clinical ---
    "disease": [
        "host_disease",
        "disease",
        "disease_status",
        "host_disease_status",
        "host_health_state",
        "host_status",
        "host_health",
        "health_status",
        "diagnosis",
        "host_condition",
        "condition",
        "host_phenotype",
        "subject_status",
        "patient_status",
        "clinical_status",
        "disease_state",
        "cohort",          # Sometimes cohort = disease group (e.g. "healthy", "NSCLC")
    ],

    # --- Intervention ---
    "intervention_group": [
        "treatment",
        "intervention",
        "group",
        "arm",
        "study_arm",
        "assigned_group",
        "cohort",          # Sometimes cohort = treatment arm
        "dose",
        "dose_group",
        "regimen",
        "diet_group",
        "fodmap_group",
        "treatment_group",
        "experimental_group",
        "donor_type",
        "inoculum",
        "fmt_source",
    ],

    "intervention_agent": [
        "intervention_agent",
        "drug",
        "medication",
        "drug_name",
        "treatment_name",
        "immunotherapy_agent",
        "probiotic_strain",
        "diet_type",
        "antibiotic",
    ],

    # --- Temporal ---
    "timepoint": [
        # Ordered: structured → descriptive → generic
        "timepoint",
        "sampling_timepoint",
        "collection_timepoint",
        "time_point",
        "study_day",
        "timepoint_day",
        "day",
        "visit",
        "visit_number",
        "week",
        "period",
        "time",
        "sampling_time",
        "collection_time",
        "sampling_date",
        "collection_date",
        "phase",
        "study_phase",
        "sample_time",
        "month",
        "time_point_number",
        "tp",
        "days_post_intervention",
        "days_after",
        "days_from_start",
    ],

    # --- Sample ---
    "body_site": [
        "isolation_source",
        "body_site",
        "sample_body_site",
        "anatomical_site",
        "tissue",
        "sample_type",
        "source",
        "host_body_site",
        "host_tissue",
        "specimen_source",
        "body_site_sampled",
        "sample_source",
    ],

    "sample_type": [
        "sample_type",
        "sample_material",
        "specimen_type",
        "material_type",
        "biosample_type",
        "source_material",
    ],

    # --- Geography ---
    "geo_location": [
        "geo_loc_name",
        "geographic_location",
        "geographic_location_country",
        "country",
        "geo_country",
        "location",
        "collection_site",
        "sampling_site",
        "region",
        "city",
        "state",
        "geo_location",
    ],

    # --- Clinical Outcome ---
    "clinical_outcome": [
        "clinical_outcome",
        "outcome",
        "response",
        "clinical_response",
        "treatment_response",
        "best_response",
        "RECIST",
        "recist",
        "PFS",
        "OS",
        "overall_survival",
        "progression_free_survival",
        "patient_outcome",
        "outcome_status",
        "responder",
        "clinical_benefit",
    ],
}

# Aliases: standard field → list of other standard fields that mean the same thing
# Used when merging mappings from different sources.
FIELD_ALIASES: dict[str, list[str]] = {
    "subject_id": ["sample_id", "participant_id"],
    "host_age": ["age", "age_years"],
    "host_sex": ["sex", "gender"],
    "host_bmi": ["bmi"],
    "body_site": ["isolation_source", "sample_site"],
    "timepoint": ["time_point", "visit", "time"],
    "intervention_group": ["group", "treatment", "arm"],
    "geo_location": ["geo_loc_name", "country"],
}


# ============================================================
# Value Normalizers
# ============================================================

def _normalize_key(key: str) -> str:
    """Normalize a raw key for comparison: lowercase, collapse separators."""
    return re.sub(r"[_\-\s]+", "_", key.strip().lower())


def _normalize_sex(value: str) -> str:
    """Normalize sex/gender values to 'male', 'female', or ''."""
    if not value:
        return ""
    v = value.strip().lower()
    if v in ("male", "m", "man", "boy", "masculine", "1"):
        return "male"
    if v in ("female", "f", "woman", "girl", "feminine", "2", "femme"):
        return "female"
    return value.strip()


def _normalize_age(value: str) -> str:
    """Normalize age: strip units, keep numeric or range."""
    if not value:
        return ""
    v = value.strip()
    # Remove common unit suffixes
    v = re.sub(r"\s*(years?|yrs?|y/o|years?_?old)\s*$", "", v, flags=re.IGNORECASE).strip()
    return v


def _normalize_bmi(value: str) -> str:
    """Normalize BMI: strip units, keep numeric."""
    if not value:
        return ""
    v = value.strip()
    v = re.sub(r"\s*kg/m[²2]\s*$", "", v, flags=re.IGNORECASE).strip()
    return v


def _normalize_timepoint(value: str) -> str:
    """Normalize timepoint value: standardize common patterns."""
    if not value:
        return ""
    v = value.strip()
    # Standardize common patterns
    v = v.replace("_", " ")
    v = re.sub(r"\s+", " ", v)
    return v


def normalize_value(field: str, value: str) -> str:
    """Normalize a raw value for a standard field.

    Args:
        field: Standard field name (e.g. 'host_sex', 'host_age').
        value: Raw value from metadata.

    Returns:
        Cleaned value string.
    """
    if not value:
        return ""

    if field in ("host_sex", "sex", "gender", "biological_sex"):
        return _normalize_sex(value)

    if field in ("host_age", "age", "age_years"):
        return _normalize_age(value)

    if field in ("host_bmi", "bmi", "body_mass_index"):
        return _normalize_bmi(value)

    if field == "timepoint":
        return _normalize_timepoint(value)

    return value.strip()


# ============================================================
# Field Detection & Matching
# ============================================================

def _build_lookup(raw_attrs: dict[str, str]) -> dict[str, str]:
    """Build a normalized-key → original-value lookup dict."""
    return {_normalize_key(k): v for k, v in raw_attrs.items()}


def detect_available_fields(raw_attrs: dict[str, str]) -> dict[str, str]:
    """Detect which standard phenotype fields are available in raw attributes.

    Args:
        raw_attrs: Raw attribute dict (key → value from ENA/NCBI).

    Returns:
        Dict of standard_field → matched_raw_key (only found fields).
    """
    lookup = _build_lookup(raw_attrs)
    found: dict[str, str] = {}

    for std_field, candidates in SYNONYMS.items():
        for candidate in candidates:
            nk = _normalize_key(candidate)
            if nk in lookup and lookup[nk]:
                value = lookup[nk].strip()
                if value and value.lower() not in (
                    "unknown", "not collected", "not provided",
                    "not applicable", "na", "n/a", "none", "missing",
                ):
                    found[std_field] = candidate  # original key, not normalized
                    break

    return found


def get_field_value(
    raw_attrs: dict[str, str],
    std_field: str,
    sample_alias: str = "",
    normalize: bool = True,
) -> str:
    """Extract the best available value for a standard field from raw attributes.

    Args:
        raw_attrs: Raw attribute dict (key → value).
        std_field: Standard field name (key in SYNONYMS).
        sample_alias: Fallback value used for subject_id when no match found.
        normalize: Apply value normalization.

    Returns:
        Best-matched value string, or '' if not found.
    """
    if std_field not in SYNONYMS:
        return ""

    lookup = _build_lookup(raw_attrs)

    for candidate in SYNONYMS[std_field]:
        nk = _normalize_key(candidate)
        if nk in lookup:
            value = lookup[nk].strip()
            if value and value.lower() not in (
                    "unknown", "not collected", "not provided",
                    "not applicable", "na", "n/a", "none", "missing",
            ):
                if normalize:
                    value = normalize_value(std_field, value)
                return value

    # Fallback for subject_id from sample_alias (the ENA sample_alias or sample_title)
    if std_field == "subject_id" and sample_alias:
        return sample_alias.strip()

    return ""


def get_field_value_with_source(
    raw_attrs: dict[str, str],
    std_field: str,
    sample_alias: str = "",
    normalize: bool = True,
) -> tuple[str, str]:
    """Like get_field_value but returns (value, matched_raw_key).

    Returns:
        Tuple of (value, matched_raw_key). matched_raw_key is the original
        ENA/NCBI attribute key that was used.
    """
    if std_field not in SYNONYMS:
        return "", ""

    lookup = _build_lookup(raw_attrs)

    for candidate in SYNONYMS[std_field]:
        nk = _normalize_key(candidate)
        if nk in lookup:
            value = lookup[nk].strip()
            if value and value.lower() not in (
                    "unknown", "not collected", "not provided",
                    "not applicable", "na", "n/a", "none", "missing",
            ):
                if normalize:
                    value = normalize_value(std_field, value)
                return value, candidate

    if std_field == "subject_id" and sample_alias:
        return sample_alias.strip(), "sample_alias"

    return "", ""


# ============================================================
# Registry Field Detection
# ============================================================

def detect_subject_id_fields(raw_attrs: dict[str, str]) -> list[str]:
    """Which raw keys likely contain subject IDs? (for registry auto-population)."""
    lookup = _build_lookup(raw_attrs)
    candidates = SYNONYMS.get("subject_id", [])[:10]  # first 10 are most specific
    found = []
    for c in candidates:
        nk = _normalize_key(c)
        if nk in lookup and lookup[nk]:
            found.append(c)
    return found


def detect_timepoint_fields(raw_attrs: dict[str, str]) -> list[str]:
    """Which raw keys likely contain timepoint info? (for registry auto-population)."""
    lookup = _build_lookup(raw_attrs)
    candidates = SYNONYMS.get("timepoint", [])
    found = []
    for c in candidates:
        nk = _normalize_key(c)
        if nk in lookup and lookup[nk]:
            val = lookup[nk].strip()
            if val and val.lower() not in ("unknown", "not collected", "not provided", "na", "n/a", ""):
                found.append(c)
    return found


def detect_host_fields(raw_attrs: dict[str, str]) -> list[str]:
    """Which ENA structured host fields are populated? (for registry auto-population)."""
    host_keys = [
        "host_sex", "host_age", "host_body_mass_index",
        "host_disease", "host_health_state", "host_body_site",
    ]
    found = []
    for key in host_keys:
        nk = _normalize_key(key)
        lookup = _build_lookup(raw_attrs)
        if nk in lookup and lookup[nk]:
            found.append(key)
    return found


# ============================================================
# Full Normalization (→ SampleAttributes)
# ============================================================

def normalize_attributes(
    raw_attrs: dict[str, str],
    sample_alias: str = "",
    sample_accession: str = "",
) -> dict:
    """Normalize a full set of raw ENA/NCBI attributes into standard fields.

    This is the primary entry point for Phase 1 → Phase 4 bridge.
    Returns a dict compatible with SampleAttributes dataclass.

    Args:
        raw_attrs: Raw key→value dict from ENA XML or NCBI BioSample.
        sample_alias: ENA sample_alias (fallback for subject_id).
        sample_accession: The sample accession (SAMN/SAMEA/SRS).

    Returns:
        Dict with keys matching SampleAttributes fields.
    """
    import json

    result = {
        "sample_accession": sample_accession,
        "subject_id_raw": get_field_value(raw_attrs, "subject_id", sample_alias),
        "age_raw": get_field_value(raw_attrs, "host_age"),
        "sex_raw": get_field_value(raw_attrs, "host_sex"),
        "bmi_raw": get_field_value(raw_attrs, "host_bmi"),
        "timepoint_raw": "",
        "timepoint_field": "",
        "group_raw": get_field_value(raw_attrs, "intervention_group"),
        "disease_raw": get_field_value(raw_attrs, "disease"),
        "body_site_raw": get_field_value(raw_attrs, "body_site"),
        "geo_raw": get_field_value(raw_attrs, "geo_location"),
        "all_attributes_json": json.dumps(raw_attrs, ensure_ascii=False),
    }

    # Timepoint: special handling — capture both value and which field was used
    tp_val, tp_field = get_field_value_with_source(raw_attrs, "timepoint")
    if not tp_val:
        # Try collection_date from ENA structured field
        collection_date = raw_attrs.get("collection_date", "").strip()
        if collection_date:
            tp_val = collection_date
            tp_field = "collection_date"
    result["timepoint_raw"] = _normalize_timepoint(tp_val)
    result["timepoint_field"] = tp_field

    return result


# ============================================================
# Synonym Table Maintenance
# ============================================================

def add_synonyms(field: str, new_keys: list[str], priority: str = "append"):
    """Programmatically extend the synonym table.

    Args:
        field: Standard field name (must exist in SYNONYMS).
        new_keys: New raw keys to add.
        priority: 'prepend' to put before existing, 'append' to put after.
    """
    if field not in SYNONYMS:
        raise KeyError(f"Unknown standard field: {field}")

    existing = SYNONYMS[field]
    for key in new_keys:
        nk = _normalize_key(key)
        # Avoid duplicates
        existing_norm = [_normalize_key(k) for k in existing]
        if nk not in existing_norm:
            if priority == "prepend":
                existing.insert(0, key)
            else:
                existing.append(key)


def export_synonyms() -> dict[str, list[str]]:
    """Export current synonym table (for debugging / documentation)."""
    return dict(SYNONYMS)


def print_synonym_summary():
    """Print a summary of the synonym table."""
    for field, keys in SYNONYMS.items():
        print(f"  {field}: {len(keys)} keys")
        for k in keys[:5]:
            print(f"    - {k}")
        if len(keys) > 5:
            print(f"    ... and {len(keys) - 5} more")