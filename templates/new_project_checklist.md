# New Project Onboarding Checklist

How to add a new metagenomics bioproject to the phenotype curation pipeline.
Time estimate: 30-60 minutes per project (depending on paper complexity).

---

## Prerequisites

- [ ] **BioProject accession** (e.g., `PRJNA1234567`) or **PMID** of the main paper
- [ ] Python 3.10+ with `pip install -e .` in `E:\ncbi\phenotype-curation\`
- [ ] `NCBI_API_KEY` env var (optional but recommended — raises rate limit from 1→10 req/s)
  - Get one: https://ncbiinsights.ncbi.nlm.nih.gov/2018/01/03/api-keys-for-enhanced-access-to-e-utilities/
- [ ] 1-2 GB free disk space (for ENA metadata downloads)

---

## Step 1: Download ENA Metadata (Phase 1)

```bash
# Fast run report (no sample attributes, ~30 seconds)
phenotype-curation download --projects PRJNA1234567 --out-dir ena_metadata/

# OR: Full sample attributes (slower, ~5 min per 500 samples, but detailed)
phenotype-curation download --projects PRJNA1234567 --sample-attrs --out-dir ena_metadata/
```

**What you get:**
- `ena_metadata/ena_metadata_master.csv` — 37-column run report
- `ena_metadata/PRJNA1234567_sample_attrs.csv` — per-sample attributes (Strategy 2 only)

**Common issues:**
- PRJNA projects may need Strategy 2 for sample attributes (ENA filereport `result=sample` returns 400 for NCBI-mirrored projects)
- If download fails mid-way, re-run with `--resume` — progress is checkpointed
- `NCBI_API_KEY` env var raises rate limit for NCBI fallback

---

## Step 2: Auto-Generate Registry (Phase 2)

```bash
phenotype-curation register ena_metadata/ena_metadata_master.csv \
  --sample-attrs-dir ena_metadata/ \
  --out project_registry.csv
```

**Review the output:**
Open `project_registry.csv` and verify:
- [ ] `study_design_candidate` is reasonable (check paper title for hints)
- [ ] `raw_group_values` contains expected group names
- [ ] `has_intervention_metadata` and `has_timepoint_metadata` are correct

**Manual additions:**
- [ ] Fill in `primary_pmid` (PubMed ID of the main paper)
- [ ] Fill in `primary_pmc` (PMC ID, if available — looks like `PMC1234567`)
- [ ] Set `scope_basis` to the reason this project is included
- [ ] Update `manual_review_status` to `pending_paper_verification`

---

## Step 3: Paper Verification (Phase 3)

```bash
# Generate classification prompts + fetch PMC full text
phenotype-curation classify \
  --registry project_registry.csv \
  --out intervention_registry.csv \
  --paper-dir paper_texts/
```

This does two things:
1. Fetches PMC full text for each paper (cached in `paper_texts/`)
2. Generates structured prompts in `intervention_registry.prompts.json`

### What you need to read from the paper

| Question | Where to look | Why it matters |
|----------|--------------|----------------|
| What was the intervention? | Abstract / Methods | Classify `intervention_type` |
| Single dose or repeated? | Methods "Study design" | Determines `is_one_time` |
| What is the sampling schedule? | Methods / Table 1 / Figure 1 | Map timepoints to phases |
| Washout/follow-up periods? | Methods / Figure 1 | Separate intervention from post |
| What are the study arms? | Methods / Table 1 | Populate `arms` and `control_type` |

### Classification decision tree

```
Is this an intervention study?
├── NO (observational cohort)
│   → study_design = "observational_cohort"
│   → intervention_type = "observational"
│   → is_one_time = ""
│   → ALL samples will be is_intervention = "no"
│
└── YES (the researchers administered something)
    │
    ├── One-time intervention (single dose/procedure)
    │   → is_one_time = "yes"
    │   → ALL samples are is_intervention = "no"
    │   → Examples: single-dose FMT, vaginal seeding at birth, HSCT transplant
    │
    └── Ongoing intervention (daily/weekly dosing, continuous diet)
        → is_one_time = "no"
        → Baseline/run-in = "no"
        → During active dosing = "yes"
        → Washout/follow-up = "no"
        → Need day-level sampling schedule
```

### Output the classification

For each project, output a JSON object matching the schema:

```json
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
```

### Apply classifications

```bash
# After filling in all LLM outputs → intervention_registry_llm_outputs.json:
phenotype-curation classify \
  --registry project_registry.csv \
  --apply intervention_registry_llm_outputs.json \
  --out intervention_registry.csv
```

---

## Step 4: Build Phenotype Table (Phase 4)

```bash
phenotype-curation build ena_metadata/ena_metadata_master.csv \
  --registry intervention_registry.csv \
  --sample-attrs-dir ena_metadata/ \
  --out phenotype_clean.csv
```

**What happens:**
- Reads registry to get `intervention_type`, `study_design`, `is_one_time`
- Extracts subject IDs from sample names using pattern matching
- Applies timepoint phase mapping (per-project phase maps)
- Determines `is_intervention` per sample
- Pairs longitudinal samples within subjects
- Outputs 36-column phenotype table

**Per-project customization needed:**
If the project has non-standard sample name conventions, add parsing rules to:
- `phase4_build.py` → `parse_subject_id()` for subject ID extraction
- `phase4_build.py` → `extract_timepoint_index()` for day/week extraction
- `phase4_build.py` → `PHASE_MAPS` for day→phase mapping (only for day-level precision)
- `phase4_build.py` → `TIMEPOINT_OVERRIDES` for label-level overrides
- `phase4_build.py` → `DISEASE_MAP` and `INTERVENTION_AGENT_MAP` for metadata backfill

---

## Step 5: Validate (Phase 5)

```bash
phenotype-curation validate phenotype_clean.csv \
  --registry intervention_registry.csv \
  --out validation_report.json
```

**7 automated checks:**
1. Row count consistency (phenotype vs expected from metadata)
2. Observational studies → all `is_intervention=no`
3. One-time interventions → all `is_intervention=no`
4. Ongoing RCTs → must have some `yes` samples
5. Phase label ↔ `is_intervention` consistency
6. No duplicate biosamples
7. Age/BMI value range sanity

**If checks fail:**
- Check 1: verify ENA metadata wasn't truncated
- Check 2: review intervention_type in registry
- Check 3: verify is_one_time classification
- Check 4: phase map may need adjustment (timepoint labeling)
- Check 5: timepoint labels not matching expected patterns
- Check 6: multiple runs per biosample (should be semicolon-joined by builder)
- Check 7: data entry errors in sample attributes

---

## Common Pitfalls

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| All samples `is_intervention=unknown` | Missing sample attribute data | Re-run Step 1 with `--sample-attrs` |
| Wrong subject IDs | Non-standard sample names | Add parsing rule to `parse_subject_id()` |
| Timepoints all "longitudinal_unspecified" | No day/week info in aliases | Need sample attrs or supplement data |
| Observational study has `is_intervention=yes` | Bug in Phase 4 logic | Check TIMEPOINT_OVERRIDES |
| Duplicate biosamples in validator | Builder not aggregating runs | Verify Phase 4 biosample dedup logic |
| PMC full text not found | Paper not in PMC (open-access) | Fall back to PubMed abstract or manual reading |
| Phase map wrong (FOS/PDX) | Day→phase mapping needs adjustment | Update PHASE_MAPS in phase4_build.py |

---

## Quick Reference: ENA API endpoints

| Purpose | Endpoint |
|---------|----------|
| Run report | `https://www.ebi.ac.uk/ena/portal/api/filereport?accession={PRJ}&result=read_run&fields=...&format=tsv` |
| Sample XML | `https://www.ebi.ac.uk/ena/browser/api/xml/{SAMN/SAMEA}` |
| NCBI fallback | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id={SAMN}&retmode=xml` |

## Quick Reference: NCBI E-utilities

| Purpose | Command |
|---------|---------|
| PMID→PMCID | `elink.fcgi?dbfrom=pubmed&db=pmc&id={PMID}` |
| PMCID→Full Text | `efetch.fcgi?db=pmc&id={PMCID}&rettype=xml` |
| PMID→Abstract | `efetch.fcgi?db=pubmed&id={PMID}&rettype=abstract&retmode=text` |