# Phenotype Curation Pipeline

**5-phase pipeline for curating biosample phenotype metadata from metagenomics studies.**

Given a list of BioProject accessions (and optionally PMIDs), produces a validated phenotype CSV where every biosample has standardized subject IDs, age, sex, BMI, disease metadata, granular timepoint phase labels (`run_in` / `feeding_1` / `washout` / ...), `is_intervention` classification, and full traceability.

## Quick Start

```bash
# Install
cd phenotype-curation
pip install -e .

# Full pipeline (with paper reading)
phenotype-curation download --projects PRJNA1234567 --sample-attrs
phenotype-curation register ena_metadata/ena_metadata_master.csv \
  --sample-attrs-dir ena_metadata/ --out project_registry.csv
phenotype-curation classify --registry project_registry.csv --out intervention_registry.csv
# → Read paper, fill intervention_registry.csv manually or via LLM
phenotype-curation build ena_metadata/ena_metadata_master.csv \
  --registry intervention_registry.csv --out phenotype_clean.csv
phenotype-curation validate phenotype_clean.csv \
  --registry intervention_registry.csv --out validation_report.json
```

## Architecture

```
[ENA/NCBI APIs]                    [Paper + Human/AI]
       │                                  │
       ▼                                  ▼
   Phase 1 ──→ Phase 2 ──→ Phase 3 ──→ Phase 4 ──→ Phase 5
   (下载)      (注册)       (论文验证)    (构建表型)    (7项校验)
```

Each phase has a **JSON Schema data contract** — output must validate before the next phase consumes it.

## CLI Commands

| Command | Phase | Description |
|---------|-------|-------------|
| `download` | 1 | Download ENA/NCBI run + sample metadata |
| `register` | 2 | Auto-generate project registry from metadata |
| `classify` | 3 | Generate paper-classification prompts (LLM-assisted) |
| `build` | 4 | Build phenotype table from registry + metadata |
| `validate` | 5 | 7-check automated validation |
| `run` | 1-3 | Full pipeline (download → classify) |
| `supplement` | util | Parse supplementary tables + detect timepoint units |

## Key Design Decisions

- **Registry-driven**: `intervention_registry.csv` is the single source of truth. Phase 4 reads it to determine `is_intervention` for every sample.
- **Biosample deduplication**: Multiple ENA runs per biosample → aggregated with semicolon-joined run_accessions.
- **Resilience**: Token bucket rate limiting (ENA 3 req/s, NCBI 1 req/s), exponential backoff with jitter, atomic checkpoint/resume for downloads.
- **Phase 3 is human-in-the-loop by design**: Paper reading requires judgment (study design classification, timing schedule extraction). The pipeline generates structured prompts, a human or LLM fills them in.

## Core Modules

| Module | Purpose |
|--------|---------|
| `core/http_client.py` | Rate-limited HTTP with exponential backoff |
| `core/checkpoint.py` | Progress tracking for resume |
| `core/schemas.py` | JSON Schema validation (4 schemas) |
| `core/normalizer.py` | 14-field synonym table for ENA/NCBI attributes |
| `core/pmid_linker.py` | PMID → BioProject resolution (4 strategies) |
| `core/supplement_parser.py` | Parse xlsx/csv/pdf supplement tables |
| `phases/phase1_fetch.py` | ENA filereport + per-sample XML download |
| `phases/phase2_registry.py` | Auto-generate project_registry.csv |
| `phases/phase3_classify.py` | PMC full-text fetch + LLM prompt generation |
| `phases/phase4_build.py` | Phenotype table construction + biosample dedup |
| `phases/phase5_validate.py` | 7 automated validation checks |

## Tested On

- **19 metagenomics studies** across diverse designs:
  - RCTs (dietary fiber, FMT, Mediterranean diet)
  - Crossover trials (FOS/PDX dose escalation)
  - Observational cohorts (cancer immunotherapy, preterm birth)
  - One-time interventions (vaginal seeding, HSCT transplant)
- **9,142 ENA runs** → **6,592 unique biosamples**
- **6/7 validation checks pass** (1 warning from missing Strategy 2 data)

## Requirements

- Python 3.10+
- `jsonschema>=4.0` (required)
- `openpyxl>=3.0`, `pdfplumber>=0.10` (optional, for `supplement` command)
- `NCBI_API_KEY` env var (optional, raises rate limit from 1→10 req/s for NCBI APIs)

## Known Limitations

- **`core/types.py` shadows stdlib `types`**: Running `python core/*.py` directly fails. Use `python -I` or the pip-installed entry point.
- **PMID → BioProject ~50% coverage**: Many papers don't list their accessions in article text. Falls back to manual SRA Run Selector lookup.
- **Phase 3 requires human judgment**: Study design classification and timepoint schedule extraction cannot be fully automated.
- **5 ongoing RCTs have 0 `is_intervention=yes`**: Needs Strategy 2 (per-sample XML) data for timepoint labels. Run `phenotype-curation download --sample-attrs`.

## License

MIT