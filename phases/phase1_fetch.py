# -*- coding: utf-8 -*-
"""Phase 1: ENA / NCBI Metadata Downloader (v2).

Downloads run-level + sample-level metadata from ENA Portal API and NCBI
BioSample for a list of bioproject accessions.

TWO STRATEGIES (auto-selected):
  Strategy 1 (fast, ~1 API call/project): ENA read_run filereport
    → run_accession, sample_accession, library info, fastq links
    → host_sex, host, host_body_site, host_status (if submitter filled them)
    → collection_date, sample_alias (often encodes subject+timepoint)
  Strategy 2 (slow, 1 API call/sample): ENA → NCBI BioSample XML attributes
    → full sample attributes (subject ID, timepoint, age, sex, BMI, treatment, etc.)
    → resumable via checkpoint file
    → automatic fallback to NCBI if ENA XML fails

RESILIENCE:
  - Exponential backoff with jitter on transient failures
  - 429 detection → Retry-After header parsing
  - Checkpoint/resume: kill mid-download, continue from where you left off
  - NCBI fallback when ENA sample XML is unavailable

Usage as library:
    from phases.phase1_fetch import fetch_metadata
    result = fetch_metadata(FetchTask(bioproject="PRJNA594393", strategy="strategy2"))

Usage as CLI:
    python phase1_fetch.py PRJNA904790 --sample-attrs
    python phase1_fetch.py --all --sample-attrs --resume
    python phase1_fetch.py --list projects.csv --out-dir ena_data/
"""

import csv
import json
import sys
import time
import argparse
import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict
from typing import Optional

# Allow running this file as a script from any directory
_PARENT = Path(__file__).resolve().parent.parent
if str(_PARENT) not in sys.path:
    sys.path.insert(0, str(_PARENT))

from core.http_client import (
    RateLimitedSession, ena_session, ncbi_session,
    FallbackFetcher, RetryConfig, ENA_PROVIDER, NCBI_PROVIDER,
)
from core.checkpoint import CheckpointManager
from core.types import (
    EnaRunRecord, SampleAttributes, Checkpoint,
    FetchTask, FetchResult, FetchStrategy,
)


# ============================================================
# CONFIGURATION
# ============================================================

# Fields to fetch from ENA read_run filereport
ENA_READ_RUN_FIELDS = [
    'run_accession', 'sample_accession', 'experiment_accession', 'study_accession',
    'sample_alias', 'sample_title', 'sample_description',
    'library_name', 'library_strategy', 'library_source', 'library_selection',
    'library_layout', 'instrument_platform', 'instrument_model',
    'fastq_ftp', 'fastq_bytes', 'fastq_md5',
    'host', 'host_sex', 'host_body_site', 'host_status', 'host_phenotype',
    'host_genotype', 'host_scientific_name', 'host_tax_id',
    'collection_date', 'collection_date_start', 'collection_date_end',
    'country', 'location',
    'scientific_name', 'tax_id',
    'first_public', 'last_updated',
    'read_count', 'base_count',
]

# Map ENA field name → EnaRunRecord attribute
ENA_FIELD_MAP = {
    'run_accession': 'run_accession',
    'sample_accession': 'sample_accession',
    'experiment_accession': 'experiment_accession',
    'study_accession': 'study_accession',
    'sample_alias': 'sample_alias',
    'sample_title': 'sample_title',
    'sample_description': 'sample_description',
    'library_name': 'library_name',
    'library_strategy': 'library_strategy',
    'library_source': 'library_source',
    'library_selection': 'library_selection',
    'library_layout': 'library_layout',
    'instrument_platform': 'instrument_platform',
    'instrument_model': 'instrument_model',
    'fastq_ftp': 'fastq_ftp',
    'fastq_bytes': 'fastq_bytes',
    'fastq_md5': 'fastq_md5',
    'host': 'host_organism',
    'host_sex': 'host_sex',
    'host_body_site': 'host_body_site',
    'host_status': 'host_status',
    'host_phenotype': 'host_phenotype',
    'host_genotype': 'host_genotype',
    'host_scientific_name': 'host_scientific_name',
    'host_tax_id': 'host_tax_id',
    'collection_date': 'collection_date',
    'collection_date_start': 'collection_date_start',
    'collection_date_end': 'collection_date_end',
    'country': 'geo_country',
    'location': 'geo_location',
    'scientific_name': 'scientific_name',
    'tax_id': 'tax_id',
    'first_public': 'ena_first_public',
    'last_updated': 'ena_last_updated',
    'read_count': 'read_count',
    'base_count': 'base_count',
}


# ============================================================
# STRATEGY 1: ENA read_run filereport
# ============================================================

def fetch_read_runs(bioproject: str, session: Optional[RateLimitedSession] = None) -> list[dict]:
    """Fetch all run metadata for a bioproject from ENA Portal API.

    Returns list of dicts with raw ENA field names as keys.
    """
    if session is None:
        session = ena_session()

    fields_str = ','.join(ENA_READ_RUN_FIELDS)
    url = (f"https://www.ebi.ac.uk/ena/portal/api/filereport"
           f"?accession={bioproject}"
           f"&result=read_run"
           f"&fields={fields_str}"
           f"&format=tsv")

    tsv_text = session.get(url)
    rows = list(csv.DictReader(tsv_text.splitlines(), delimiter='\t'))
    return rows


# ============================================================
# STRATEGY 2: Per-sample XML attributes
# ============================================================

def _parse_ena_sample_xml(xml_text: str) -> dict:
    """Parse ENA sample XML into a flat dict of TAG→VALUE pairs."""
    attrs = {}
    try:
        root = ET.fromstring(xml_text)
        for attr in root.iter('SAMPLE_ATTRIBUTE'):
            tag_el = attr.find('TAG')
            value_el = attr.find('VALUE')
            if tag_el is not None and value_el is not None:
                attrs[tag_el.text] = value_el.text if value_el.text else ''
    except ET.ParseError:
        pass
    return attrs


def _parse_ncbi_biosample_xml(xml_text: str) -> dict:
    """Parse NCBI BioSample efetch XML into a flat dict.

    NCBI BioSample XML wraps attributes differently than ENA:
      <Attribute attribute_name="..." harmonized_name="...">value</Attribute>
    """
    attrs = {}
    try:
        root = ET.fromstring(xml_text)
        for attr in root.iter('Attribute'):
            attr_name = attr.attrib.get('attribute_name', '')
            if attr_name:
                attrs[attr_name] = attr.text if attr.text else ''
            # Also capture harmonized name
            harm = attr.attrib.get('harmonized_name', '')
            if harm and harm != attr_name:
                attrs[harm] = attr.text if attr.text else ''
    except ET.ParseError:
        pass
    return attrs


def normalize_sample_attrs(raw_attrs: dict, sample_alias: str = '') -> SampleAttributes:
    """Normalize raw attribute dict into SampleAttributes dataclass."""
    result = SampleAttributes()
    raw_lower = {k.lower().replace(' ', '_'): v for k, v in raw_attrs.items()}

    # Subject ID
    for key in ['subject', 'host_subject_id', 'subject_id', 'individual', 'patient', 'donor']:
        if key in raw_lower and raw_lower[key]:
            result.subject_id_raw = raw_lower[key]
            break
    if not result.subject_id_raw:
        result.subject_id_raw = sample_alias

    # Age
    for key in ['host_age', 'age', 'age_years', 'age_at_collection']:
        if key in raw_lower and raw_lower[key]:
            result.age_raw = raw_lower[key]
            break

    # Sex
    for key in ['host_sex', 'sex', 'gender']:
        if key in raw_lower and raw_lower[key]:
            result.sex_raw = raw_lower[key]
            break

    # BMI
    for key in ['host_body_mass_index', 'bmi', 'body_mass_index']:
        if key in raw_lower and raw_lower[key]:
            result.bmi_raw = raw_lower[key]
            break

    # Timepoint — prefer meaningful values over "unknown"
    timepoint_value = ''
    timepoint_field = ''
    for key in ['timepoint', 'time_point', 'visit', 'week', 'day',
                 'period', 'study_day', 'sampling_timepoint', 'time',
                 'collection_date', 'sampling_date']:
        if key in raw_lower:
            val = raw_lower[key]
            if val and val.lower() not in ('unknown', 'not collected', 'not provided', 'na', ''):
                timepoint_value = val
                timepoint_field = key
                break
            elif not timepoint_value:
                timepoint_value = val
                timepoint_field = key
    result.timepoint_raw = timepoint_value
    result.timepoint_field = timepoint_field

    # Treatment / group
    for key in ['treatment', 'intervention', 'group', 'arm', 'cohort']:
        if key in raw_lower and raw_lower[key]:
            result.group_raw = raw_lower[key]
            break

    # Disease
    for key in ['host_disease', 'disease', 'disease_status', 'host_health_state']:
        if key in raw_lower and raw_lower[key]:
            result.disease_raw = raw_lower[key]
            break

    # Body site
    for key in ['isolation_source', 'body_site', 'tissue', 'sample_type']:
        if key in raw_lower and raw_lower[key]:
            result.body_site_raw = raw_lower[key]
            break

    # Geography
    if 'geo_loc_name' in raw_lower:
        result.geo_raw = raw_lower['geo_loc_name']

    result.all_attributes_json = json.dumps(raw_attrs, ensure_ascii=False)

    return result


def fetch_sample_attributes_ena(sample_accession: str, session: RateLimitedSession) -> str:
    """Fetch sample XML from ENA Browser API. Returns XML text."""
    url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{sample_accession}"
    return session.get(url, RetryConfig(max_retries=3, base_delay=2.0), timeout=60)


def fetch_sample_attributes_ncbi(sample_accession: str, session: RateLimitedSession) -> str:
    """Fetch sample XML from NCBI E-utilities BioSample. Returns XML text."""
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=biosample&id={sample_accession}&rettype=xml&retmode=xml")
    return session.get(url, RetryConfig(max_retries=3, base_delay=2.0), timeout=60)


def fetch_sample_with_fallback(sample_accession: str,
                               ena_sess: RateLimitedSession,
                               ncbi_sess: RateLimitedSession) -> tuple[dict, str]:
    """Fetch sample attributes with ENA primary → NCBI fallback.

    Returns (raw_attrs_dict, source_label).
    """
    fetcher = FallbackFetcher(
        primary=lambda: _parse_ena_sample_xml(
            fetch_sample_attributes_ena(sample_accession, ena_sess)),
        fallback=lambda: _parse_ncbi_biosample_xml(
            fetch_sample_attributes_ncbi(sample_accession, ncbi_sess)),
        name=f"sample {sample_accession}",
    )
    return fetcher.fetch()


# ============================================================
# MERGE + OUTPUT
# ============================================================

def merge_run_and_sample(run_rows: list[dict],
                         sample_attrs: Optional[list[SampleAttributes]] = None) -> list[dict]:
    """Left-join run metadata with sample attributes."""
    sample_map = {}
    if sample_attrs:
        sample_map = {s.sample_accession: s for s in sample_attrs
                      if s.sample_accession}

    merged = []
    for run in run_rows:
        row = {'bioproject': run.get('study_accession', '')}
        for ena_field, attr_name in ENA_FIELD_MAP.items():
            row[attr_name] = run.get(ena_field, '')

        sa = sample_map.get(run.get('sample_accession', ''))
        if sa:
            row['subject_id_raw'] = sa.subject_id_raw
            row['age_raw'] = sa.age_raw
            row['sex_raw'] = sa.sex_raw
            row['bmi_raw'] = sa.bmi_raw
            row['timepoint_raw'] = sa.timepoint_raw
            row['timepoint_field'] = sa.timepoint_field
            row['group_raw'] = sa.group_raw
            row['disease_raw'] = sa.disease_raw
            row['body_site_raw'] = sa.body_site_raw
            row['geo_raw'] = sa.geo_raw
            row['all_attributes_json'] = sa.all_attributes_json

        merged.append(row)

    return merged


def write_merged_csv(rows: list[dict], output_dir: Path, bioproject: str,
                     include_sample_attrs: bool = False) -> Path:
    """Write per-project merged CSV. Returns output path."""
    out_path = output_dir / f"{bioproject}_merged.csv"

    base_fields = ['bioproject'] + list(ENA_FIELD_MAP.values())
    sample_fields = ['subject_id_raw', 'age_raw', 'sex_raw', 'bmi_raw',
                     'timepoint_raw', 'timepoint_field', 'group_raw',
                     'disease_raw', 'body_site_raw', 'geo_raw',
                     'all_attributes_json']
    if include_sample_attrs:
        base_fields.extend(sample_fields)

    with open(out_path, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=base_fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

    return out_path


# ============================================================
# MAIN FETCH ORCHESTRATOR
# ============================================================

def fetch_metadata(task: FetchTask, output_dir: Path) -> FetchResult:
    """Download metadata for one bioproject.

    Handles both strategies, checkpointing, and fallback.
    """
    result = FetchResult(bioproject=task.bioproject)
    ena_sess = ena_session()
    ncbi_sess = ncbi_session()
    ck_mgr = CheckpointManager(output_dir)

    try:
        # Phase 1a: fetch run metadata (always)
        print(f"  [Strategy 1] Fetching run metadata...")
        run_rows = fetch_read_runs(task.bioproject, ena_sess)
        result.runs_downloaded = len(run_rows)
        print(f"    → {len(run_rows)} runs")

        # Save raw TSV for reference
        raw_path = output_dir / f"{task.bioproject}_runs.tsv"
        # Already downloaded by fetch_read_runs — the TSV is in memory
        # We'll save it from the merged output

        # Phase 1b: fetch sample attributes (if strategy 2)
        sample_attrs = None
        if task.strategy == FetchStrategy.STRATEGY2.value:
            sample_accessions = [r['sample_accession'] for r in run_rows
                                 if r.get('sample_accession')]
            if task.max_samples:
                sample_accessions = sample_accessions[:task.max_samples]

            total = len(sample_accessions)
            result.samples_downloaded = 0

            # Load or create checkpoint
            ck = ck_mgr.load_or_create(
                task.bioproject,
                total_samples=total,
                strategy=FetchStrategy.STRATEGY2.value,
            )
            # Update total if grown (e.g., max-samples increased)
            if total > ck.total_samples:
                ck.total_samples = total

            if task.resume and ck.done_count() > 0:
                print(f"  [Strategy 2] Resuming: {ck.done_count()}/{total} already done, "
                      f"{ck.failed_count()} failed, {ck.remaining()} remaining")

            sample_attrs = []
            fetched_count = 0
            done_set: set = set(ck.completed)

            # On resume, load previously-saved sample attributes
            sat_csv_path = output_dir / f"{task.bioproject}_sample_attrs.csv"
            if task.resume and ck.done_count() > 0 and sat_csv_path.exists():
                with open(sat_csv_path, encoding='utf-8-sig') as f:
                    for row in csv.DictReader(f):
                        sa = SampleAttributes(
                            sample_accession=row.get('sample_accession', ''),
                            subject_id_raw=row.get('subject_id_raw', ''),
                            age_raw=row.get('age_raw', ''),
                            sex_raw=row.get('sex_raw', ''),
                            bmi_raw=row.get('bmi_raw', ''),
                            timepoint_raw=row.get('timepoint_raw', ''),
                            timepoint_field=row.get('timepoint_field', ''),
                            group_raw=row.get('group_raw', ''),
                            disease_raw=row.get('disease_raw', ''),
                            body_site_raw=row.get('body_site_raw', ''),
                            geo_raw=row.get('geo_raw', ''),
                            all_attributes_json=row.get('all_attributes_json', ''),
                        )
                        sample_attrs.append(sa)
                        fetched_count += 1
                print(f"    Loaded {fetched_count} previously-fetched sample attributes from checkpoint")

            for i, acc in enumerate(sample_accessions):
                # Skip if already completed (loaded from CSV or previously done)
                if acc in done_set:
                    continue

                if (i + 1) % 50 == 0 or i == 0:
                    print(f"    [{i+1}/{total}] fetching sample attributes...")

                try:
                    raw_attrs, source = fetch_sample_with_fallback(
                        acc, ena_sess, ncbi_sess)
                    sa = normalize_sample_attrs(raw_attrs,
                                               run_rows[i].get('sample_alias', ''))
                    sa.sample_accession = acc
                    sample_attrs.append(sa)

                    ck_mgr.mark_done(ck, acc)
                    fetched_count += 1

                    # Save checkpoint every 50 samples
                    if fetched_count % 50 == 0:
                        ck_mgr.save(ck)

                except Exception as e:
                    print(f"    WARN: Failed {acc} after fallback: {e}")
                    ck_mgr.mark_failed(ck, acc, str(e))
                    result.samples_failed += 1
                    # Add placeholder
                    sample_attrs.append(SampleAttributes(sample_accession=acc))

                # Rate limit courtesy
                time.sleep(0.2)

            # Final save
            ck_mgr.save(ck)
            result.samples_downloaded = fetched_count
            print(f"    → {fetched_count} sample attributes fetched "
                  f"({ck.failed_count()} failed)")

            # Save sample attributes CSV for reference
            sat_path = output_dir / f"{task.bioproject}_sample_attrs.csv"
            with open(sat_path, 'w', newline='', encoding='utf-8-sig') as f:
                fields = ['sample_accession', 'subject_id_raw', 'age_raw',
                          'sex_raw', 'bmi_raw', 'timepoint_raw', 'timepoint_field',
                          'group_raw', 'disease_raw', 'body_site_raw', 'geo_raw',
                          'all_attributes_json']
                writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
                writer.writeheader()
                writer.writerows([s.__dict__ for s in sample_attrs])

        # Merge and write
        merged = merge_run_and_sample(run_rows, sample_attrs)
        out_path = write_merged_csv(
            merged, output_dir, task.bioproject,
            include_sample_attrs=(sample_attrs is not None),
        )
        result.output_file = str(out_path)
        result.success = True
        print(f"    → {len(merged)} merged rows → {out_path.name}")

    except Exception as e:
        result.success = False
        result.errors.append(str(e))
        print(f"  ERROR: {e}")

    return result


# ============================================================
# CLI
# ============================================================

KNOWN_BIOPROJECTS = [
    'PRJEB22863', 'PRJEB33500', 'PRJEB39610', 'PRJEB48663', 'PRJEB75958',
    'PRJNA397450', 'PRJNA541981',
    'PRJNA594393', 'PRJNA594610', 'PRJNA594613', 'PRJNA594620',
    'PRJNA637785', 'PRJNA701480', 'PRJNA737472', 'PRJNA762360',
    'PRJNA895487', 'PRJNA904790', 'PRJNA929328', 'PRJNA929986',
]


def main():
    parser = argparse.ArgumentParser(
        description='Phase 1: ENA Metadata Downloader (v2)')
    parser.add_argument('bioprojects', nargs='*',
                        help='Bioproject accessions')
    parser.add_argument('--list', '-l',
                        help='CSV file with "bioproject" column')
    parser.add_argument('--all', action='store_true',
                        help=f'Download all {len(KNOWN_BIOPROJECTS)} known projects')
    parser.add_argument('--sample-attrs', '-s', action='store_true',
                        help='Strategy 2: per-sample XML attributes (resumable, with NCBI fallback)')
    parser.add_argument('--resume', '-r', action='store_true',
                        help='Resume from checkpoint (skip already-fetched samples)')
    parser.add_argument('--out-dir', '-o', default='ena_metadata',
                        help='Output directory')
    parser.add_argument('--max-samples', type=int, default=0,
                        help='Max samples per project (0=all, for testing)')

    args = parser.parse_args()

    # Resolve bioprojects
    bioprojects = args.bioprojects or []
    if args.list:
        with open(args.list, encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            bioprojects.extend(r.get('bioproject', r.get('bioproject_accession', ''))
                             for r in reader)
    if args.all:
        bioprojects = list(KNOWN_BIOPROJECTS)

    bioprojects = [bp.strip() for bp in bioprojects if bp.strip()]
    if not bioprojects:
        print("ERROR: No bioprojects. Use --list, --all, or pass accessions.")
        sys.exit(1)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    strategy = FetchStrategy.STRATEGY2.value if args.sample_attrs else FetchStrategy.STRATEGY1.value

    print(f"ENA Metadata Downloader v2")
    print(f"  Projects: {len(bioprojects)}")
    print(f"  Strategy: {'Strategy 2 (run + sample XML with NCBI fallback)'
                        if args.sample_attrs else 'Strategy 1 (run metadata only)'}")
    print(f"  Output:   {out_dir.absolute()}")
    if args.resume:
        print(f"  Resume:   enabled (will skip completed samples)")
    print()

    all_merged = []
    success, failed = 0, 0

    for bp in bioprojects:
        print(f"[{bp}]")
        task = FetchTask(
            bioproject=bp,
            strategy=strategy,
            max_samples=args.max_samples,
            resume=args.resume,
        )
        result = fetch_metadata(task, out_dir)
        if result.success:
            success += 1
            if result.output_file:
                # Read merged file back to accumulate for master
                merged_path = Path(result.output_file)
                if merged_path.exists():
                    with open(merged_path, encoding='utf-8-sig') as f:
                        all_merged.extend(list(csv.DictReader(f)))
        else:
            failed += 1
        print()

    # Write master CSV
    if all_merged:
        # Union of all field names
        all_fields = list(dict.fromkeys(
            k for row in all_merged for k in row))

        master_path = out_dir / "ena_metadata_master.csv"
        with open(master_path, 'w', newline='', encoding='utf-8-sig') as f:
            writer = csv.DictWriter(f, fieldnames=all_fields, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(all_merged)
        print(f"Master file: {master_path.absolute()} "
              f"({len(all_merged)} rows, {len(all_fields)} columns)")

    # Summary
    print(f"\n{'='*60}")
    print(f"Done: {success} succeeded, {failed} failed")
    print(f"Output: {out_dir.absolute()}")

    # Coverage report
    if all_merged:
        print(f"\nMetadata coverage ({len(all_merged)} runs):")
        coverage = {
            'host_sex': 'host_sex filled',
            'host_body_site': 'host_body_site filled',
            'collection_date': 'collection_date not "unknown"',
            'subject_id_raw': 'Subject ID extracted',
            'age_raw': 'Age extracted',
            'timepoint_raw': 'Timepoint extracted',
            'group_raw': 'Group/treatment extracted',
        }
        for field, desc in coverage.items():
            if field in all_merged[0]:
                count = sum(1 for r in all_merged
                          if r.get(field) and r[field] not in
                          ('', 'unknown', 'not collected', 'not provided'))
                pct = count / len(all_merged) * 100
                print(f"  {desc:35s}: {count:5d}/{len(all_merged)} ({pct:.0f}%)")


if __name__ == '__main__':
    main()