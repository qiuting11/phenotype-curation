# -*- coding: utf-8 -*-
"""Core dataclasses for phenotype curation pipeline."""

from dataclasses import dataclass, field, asdict
from datetime import datetime
from enum import Enum
from typing import Optional


# ============================================================
# Enums
# ============================================================

class StudyDesign(str, Enum):
    OBSERVATIONAL_COHORT = "observational_cohort"
    RCT = "RCT"
    RCT_CROSSOVER = "RCT_crossover"
    DOSE_ESCALATION = "dose_escalation"


class InterventionType(str, Enum):
    OBSERVATIONAL = "observational"
    IMMUNOTHERAPY = "immunotherapy"
    DIET = "diet"
    DIET_PREBIOTIC = "diet_prebiotic"
    FMT = "fmt"
    PROBIOTIC = "probiotic"
    VAGINAL_SEEDING = "vaginal_seeding"


class EvidenceQuality(str, Enum):
    VERIFIED = "verified"
    GOOD = "good"
    MODERATE = "moderate"
    LOW = "low"


class IsIntervention(str, Enum):
    YES = "yes"
    NO = "no"
    UNKNOWN = "unknown"


class FetchStrategy(str, Enum):
    STRATEGY1 = "strategy1"   # filereport only (fast)
    STRATEGY2 = "strategy2"   # + per-sample XML (slow, detailed)


# ============================================================
# Phase 1: ENA Metadata
# ============================================================

@dataclass
class EnaRunRecord:
    """One row from ENA read_run filereport."""
    bioproject: str = ""
    run_accession: str = ""
    sample_accession: str = ""
    experiment_accession: str = ""
    study_accession: str = ""
    sample_alias: str = ""
    sample_title: str = ""
    sample_description: str = ""
    # Library
    library_name: str = ""
    library_strategy: str = ""
    library_source: str = ""
    library_selection: str = ""
    library_layout: str = ""
    instrument_platform: str = ""
    instrument_model: str = ""
    # FASTQ
    fastq_ftp: str = ""
    fastq_bytes: str = ""
    fastq_md5: str = ""
    # Host (ENA structured — often empty)
    host_organism: str = ""
    host_sex: str = ""
    host_body_site: str = ""
    host_status: str = ""
    host_phenotype: str = ""
    host_genotype: str = ""
    host_scientific_name: str = ""
    host_tax_id: str = ""
    # Temporal
    collection_date: str = ""
    collection_date_start: str = ""
    collection_date_end: str = ""
    # Geography
    geo_country: str = ""
    geo_location: str = ""
    # Taxonomy
    scientific_name: str = ""
    tax_id: str = ""
    # Timestamps
    ena_first_public: str = ""
    ena_last_updated: str = ""
    # Stats
    read_count: str = ""
    base_count: str = ""


@dataclass
class SampleAttributes:
    """Extracted from ENA/NCBI BioSample XML attributes."""
    sample_accession: str = ""
    subject_id_raw: str = ""
    age_raw: str = ""
    sex_raw: str = ""
    bmi_raw: str = ""
    timepoint_raw: str = ""
    timepoint_field: str = ""
    group_raw: str = ""
    disease_raw: str = ""
    body_site_raw: str = ""
    geo_raw: str = ""
    all_attributes_json: str = ""   # full raw attributes for reference


# ============================================================
# Phase 2: Project Registry
# ============================================================

@dataclass
class ProjectRegistryRow:
    """project_registry.csv — auto-generated from metadata, manually reviewed."""
    bioproject: str = ""
    keep_flag: str = "TRUE"
    scope_basis: str = ""
    primary_pmid: str = ""
    additional_pmids: str = ""
    paper_count: str = ""
    run_n_local: str = ""
    biosample_n_local: str = ""
    sample_name_n_local: str = ""
    has_intervention_metadata: str = ""
    has_timepoint_metadata: str = ""
    raw_group_values: str = ""
    raw_intervention_values: str = ""
    raw_disease_values: str = ""
    study_design_candidate: str = ""
    training_readiness_candidate: str = ""
    manual_review_status: str = "pending_paper_verification"
    supplement_status: str = "to_be_verified"
    # Auto-detected fields
    common_subject_id_fields: str = ""
    common_timepoint_fields: str = ""
    available_host_fields: str = ""


# ============================================================
# Phase 3: Intervention Registry
# ============================================================

@dataclass
class InterventionRegistryRow:
    """project_intervention_registry.csv — paper-verified study design."""
    bioproject: str = ""
    study_design: str = ""
    intervention_type: str = ""
    arms: str = ""
    control_type: str = ""
    randomization: str = "no"
    blinding: str = ""
    sample_timing: str = ""
    primary_pmid: str = ""
    primary_pmc: str = ""
    notes: str = ""
    # Phase-3-added fields
    is_one_time: str = ""             # "yes" / "no" — single dose/procedure?
    sample_timing_schedule: str = ""  # human-readable: "day -14, 0, 7, 14, 28, 56"
    days_since_intervention: str = "" # "[0, 14, 28]" — numeric days from intervention start
    classification_confidence: str = EvidenceQuality.GOOD.value


# ============================================================
# Phase 4: Phenotype Table
# ============================================================

@dataclass
class PhenotypeRow:
    """phenotype_clean.csv — one row per biosample."""
    # Identity
    bioproject: str = ""
    biosample_id: str = ""
    run_accession: str = ""
    subject_id: str = ""
    sample_id_local: str = ""
    sample_type: str = ""
    collection_date: str = ""
    library_layout: str = ""
    platform: str = ""
    # Phenotype
    disease: str = ""
    intervention_type: str = ""
    intervention_agent: str = ""
    intervention_group: str = ""
    intervention_source: str = ""
    host_age: str = ""
    host_sex: str = ""
    host_bmi: str = ""
    ethnicity: str = ""
    clinical_outcome: str = ""
    phenotype_source: str = ""
    # Timeline
    timepoint_label: str = ""
    timepoint_index: str = ""
    timepoint_unit: str = ""
    is_intervention: str = IsIntervention.UNKNOWN.value
    # Lineage
    paired_pre_sample: str = ""
    paired_post_sample: str = ""
    sample_seq: str = ""
    total_samples: str = ""
    pairing_confidence: str = ""
    parallel_samples: str = ""
    # Context
    context_role: str = ""
    context_id: str = ""
    family_id: str = ""
    notes: str = ""
    # Provenance
    subject_id_status: str = ""
    subject_id_confidence: str = ""

    # Computed columns for display order
    OUTPUT_COLUMNS = [
        "bioproject", "biosample_id", "run_accession", "subject_id", "sample_id_local",
        "sample_type", "collection_date", "library_layout", "platform",
        "disease", "intervention_type", "intervention_agent", "intervention_group",
        "intervention_source", "host_age", "host_sex", "host_bmi", "ethnicity",
        "clinical_outcome", "phenotype_source",
        "timepoint_label", "timepoint_index", "timepoint_unit", "is_intervention",
        "paired_pre_sample", "paired_post_sample", "sample_seq", "total_samples",
        "pairing_confidence", "parallel_samples",
        "context_role", "context_id", "family_id", "notes",
        "subject_id_status", "subject_id_confidence",
    ]


# ============================================================
# Phase 5: Traceability
# ============================================================

@dataclass
class TraceabilityRow:
    """is_intervention_traceability.csv — decision audit trail."""
    bioproject: str = ""
    intervention_type: str = ""
    study_design: str = ""
    paper_pmcid: str = ""
    key_design_finding: str = ""
    is_intervention_rationale: str = ""
    is_intervention_yes: str = ""
    is_intervention_no: str = ""
    is_intervention_unknown: str = ""
    evidence_quality: str = EvidenceQuality.GOOD.value


# ============================================================
# Infrastructure Types
# ============================================================

@dataclass
class Checkpoint:
    """Progress checkpoint for resumable downloads."""
    bioproject: str = ""
    strategy: str = FetchStrategy.STRATEGY1.value
    total_samples: int = 0
    completed: list = field(default_factory=list)       # list of sample accessions
    failed: dict = field(default_factory=dict)           # {accession: error_reason}
    last_updated: str = field(default_factory=lambda: datetime.now().isoformat())
    output_dir: str = ""

    def done_count(self) -> int:
        return len(self.completed)

    def failed_count(self) -> int:
        return len(self.failed)

    def remaining(self) -> int:
        return max(0, self.total_samples - self.done_count())

    def is_complete(self) -> bool:
        return self.done_count() >= self.total_samples > 0


@dataclass
class FetchTask:
    """What to download for a single bioproject."""
    bioproject: str = ""
    strategy: str = FetchStrategy.STRATEGY1.value
    max_samples: int = 0    # 0 = all
    resume: bool = False


@dataclass
class FetchResult:
    """Result of a fetch operation (one project)."""
    bioproject: str = ""
    strategy: str = ""
    runs_downloaded: int = 0
    samples_downloaded: int = 0
    samples_failed: int = 0
    output_file: str = ""
    errors: list = field(default_factory=list)
    success: bool = False


@dataclass
class ValidationReport:
    """Phase 5 output."""
    total_rows: int = 0
    expected_rows: int = 0
    checks: list = field(default_factory=list)   # [{name, passed, details}]
    is_intervention_distribution: dict = field(default_factory=dict)
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    passed: bool = False