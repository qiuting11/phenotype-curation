#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phenotype Curation CLI — unified entry point for the 5-phase pipeline.

Usage:
    phenotype-curation download  [--projects ...] [--sample-attrs] [--resume]
    phenotype-curation register  [--metadata ...] [--sample-attrs-dir ...]
    phenotype-curation classify  [--registry ...] [--paper-dir ...]
    phenotype-curation build     [--metadata ...] [--registry ...]
    phenotype-curation validate  [--phenotype ...] [--registry ...]
    phenotype-curation run       [--metadata ...] [--registry ...]  # Full pipeline

Each subcommand wraps one phase module and enforces input/output schema contracts.
"""

import argparse
import sys
from pathlib import Path

# Add parent to path for import
_PARENT = Path(__file__).resolve().parent
if str(_PARENT) not in sys.path:
    sys.path.insert(0, str(_PARENT))


def _path(value: str) -> Path:
    """Argparse type: coerce to Path."""
    return Path(value)


# ============================================================
# Subcommand: download (Phase 1)
# ============================================================

def cmd_download(args):
    """Phase 1: Download ENA/NCBI metadata."""
    from phases.phase1_fetch import fetch_metadata, FetchTask, FetchStrategy

    # Determine projects
    projects = []
    if args.projects:
        projects = [p.strip() for p in args.projects.split(",") if p.strip()]
    elif args.all:
        # Would read from registry or config
        print("ERROR: --all requires a project list file or registry. Not yet implemented.")
        sys.exit(1)
    else:
        print("ERROR: Specify --projects or --all")
        sys.exit(1)

    out_dir = args.out_dir or Path(f"ena_metadata_{args.strategy or 'auto'}")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    strategy = args.strategy or (
        FetchStrategy.STRATEGY2.value if args.sample_attrs
        else FetchStrategy.STRATEGY1.value
    )

    for bp in projects:
        task = FetchTask(
            bioproject=bp,
            strategy=strategy,
            max_samples=args.max_samples or 0,
            resume=args.resume,
        )
        result = fetch_metadata(task, out_dir)
        status = "OK" if result.success else "FAILED"
        print(f"  {bp}: {status} — {result.runs_downloaded} runs, "
              f"{result.samples_downloaded} samples")
        if result.errors:
            for e in result.errors[:3]:
                print(f"    {e}")


# ============================================================
# Subcommand: register (Phase 2)
# ============================================================

def cmd_register(args):
    """Phase 2: Auto-generate project registry from metadata."""
    from phases.phase2_registry import build_registry

    build_registry(
        metadata_path=args.metadata,
        sample_attrs_dir=args.sample_attrs_dir,
        out_path=args.out or "project_registry.csv",
    )


# ============================================================
# Subcommand: classify (Phase 3)
# ============================================================

def cmd_classify(args):
    """Phase 3: Paper-verify study design (LLM-assisted).

    This phase REQUIRES human/AI interaction — it provides a structured
    prompt per paper for classification but does not automatically
    classify without paper reading.
    """
    from phases.phase3_classify import classify_projects, apply_classifications
    import json

    if args.apply:
        # Apply mode: read LLM outputs and write registry
        with open(args.apply, "r", encoding="utf-8") as f:
            llm_outputs = json.load(f)
        prompts_path = Path(args.out).with_suffix(".prompts.json") if args.out else Path("classify_prompts.json")
        prompts = []
        if prompts_path.exists():
            with open(prompts_path, "r", encoding="utf-8") as f:
                prompts = json.load(f)
        apply_classifications(prompts, llm_outputs, args.out or "intervention_registry.csv")
    else:
        # Generate mode: prepare prompts and optionally fetch paper text
        classify_projects(
            registry_path=args.registry,
            out_path=args.out or "intervention_registry.csv",
            fetch_papers=not args.no_fetch,
            paper_dir=args.paper_dir or "paper_texts",
        )


# ============================================================
# Subcommand: build (Phase 4)
# ============================================================

def cmd_build(args):
    """Phase 4: Build phenotype table from registry + metadata."""
    from phases.phase4_build import build_phenotype_table

    build_phenotype_table(
        metadata_path=args.metadata,
        registry_path=args.registry,
        sample_attrs_dir=args.sample_attrs_dir,
        out_path=args.out or "phenotype_clean.csv",
    )


# ============================================================
# Subcommand: validate (Phase 5)
# ============================================================

def cmd_validate(args):
    """Phase 5: 7-check automated validation."""
    from phases.phase5_validate import validate

    report = validate(
        phenotype_path=args.phenotype,
        registry_path=args.registry,
        out_path=args.out or "validation_report.json",
    )

    if not report.passed:
        sys.exit(1)


# ============================================================
# Subcommand: run (Full Pipeline)
# ============================================================

def cmd_run(args):
    """Run the full 5-phase pipeline end-to-end."""
    print("=" * 60)
    print("PHENOTYPE CURATION — FULL PIPELINE")
    print("=" * 60)

    # Phase 1: Download
    print("\n[Phase 1] Download ENA metadata...")
    # Delegate to cmd_download with same args
    cmd_download(args)

    # Phase 2: Register
    print("\n[Phase 2] Generate project registry...")
    cmd_register(args)

    # Phase 3: Classify (requires human)
    print("\n[Phase 3] Classification prompts generated.")
    print("  ⚠ Phase 3 requires human/AI paper reading.")
    print("  Complete intervention_registry.csv then re-run from Phase 4.")
    cmd_classify(args)

    # Phase 4+5 are run separately after Phase 3 completion
    print("\nAfter Phase 3 completion, run:")
    print("  phenotype-curation build --metadata <dir> --registry intervention_registry.csv")
    print("  phenotype-curation validate --phenotype phenotype_clean.csv --registry intervention_registry.csv")


# ============================================================
# Subcommand: supplement (Utility)
# ============================================================

def cmd_supplement(args):
    """Parse supplementary tables from papers."""
    from core.supplement_parser import SupplementParser, detect_timepoint_unit

    parser = SupplementParser()

    if args.all_tables:
        results = parser.parse_all_tables(
            args.path,
            header_row=args.header_row,
            sheet_name=args.sheet or "",
            skip_rows=args.skip_rows,
        )
        for i, r in enumerate(results):
            sheet_info = f" [{r.sheet_name}]" if r.sheet_name else f" [table {r.table_index}]"
            print(f"{i + 1}.{sheet_info} {r.summary()}")
    else:
        result = parser.parse(
            args.path,
            header_row=args.header_row,
            sheet_name=args.sheet or "",
            table_index=args.table_index,
            skip_rows=args.skip_rows,
        )
        if args.report:
            print(parser.export_detection_report(result))
        else:
            print(result.summary())

        # Unit detection for timepoint columns (Gap #9)
        if args.detect_units:
            tp_values = result.column_values("timepoint")
            if tp_values and any(v.strip() for v in tp_values):
                unit_detection = detect_timepoint_unit(tp_values)
                print(f"  Timepoint unit: {unit_detection.unit} "
                      f"({unit_detection.confidence}) — {unit_detection.evidence}")
            else:
                print("  Timepoint unit: no timepoint column detected")

        if args.out:
            parser.export_csv(result, args.out)
            print(f"Exported normalized CSV → {args.out}")


# ============================================================
# Main CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Phenotype Curation — 5-phase pipeline for metagenomics biosample curation",
        prog="phenotype-curation",
    )
    subparsers = parser.add_subparsers(dest="command", help="Phase to execute")

    # --- download ---
    p_dl = subparsers.add_parser("download", help="Phase 1: Download ENA/NCBI metadata")
    p_dl.add_argument("--projects", help="Comma-separated bioproject accessions")
    p_dl.add_argument("--all", action="store_true", help="Download all projects from registry")
    p_dl.add_argument("--sample-attrs", action="store_true",
                      help="Use Strategy 2 (per-sample XML, slow but detailed)")
    p_dl.add_argument("--strategy", choices=["strategy1", "strategy2"],
                      help="Override strategy selection")
    p_dl.add_argument("--resume", action="store_true", help="Resume from checkpoint")
    p_dl.add_argument("--max-samples", type=int, help="Max samples per project (for testing)")
    p_dl.add_argument("--out-dir", help="Output directory")
    p_dl.add_argument("--list", action="store_true", help="List available projects in registry")

    # --- register ---
    p_reg = subparsers.add_parser("register", help="Phase 2: Auto-generate project registry")
    p_reg.add_argument("metadata", help="Path to ena_metadata_master.csv or directory")
    p_reg.add_argument("--sample-attrs-dir", help="Directory with *_sample_attrs.csv")
    p_reg.add_argument("--out", help="Output CSV path")

    # --- classify ---
    p_cls = subparsers.add_parser("classify",
                                  help="Phase 3: Generate paper-classification prompts (LLM-assisted)")
    p_cls.add_argument("--registry", required=True, help="Path to project_registry.csv")
    p_cls.add_argument("--paper-dir", default="paper_texts", help="Directory to cache paper text")
    p_cls.add_argument("--no-fetch", action="store_true", help="Skip PMC full text fetch")
    p_cls.add_argument("--out", help="Output CSV path (default: intervention_registry.csv)")
    p_cls.add_argument("--apply", help="JSON file mapping bioproject→LLM output (apply mode)")

    # --- build ---
    p_bld = subparsers.add_parser("build", help="Phase 4: Build phenotype table")
    p_bld.add_argument("metadata", help="Path to ena_metadata_master.csv")
    p_bld.add_argument("--registry", required=True, help="Path to intervention_registry.csv")
    p_bld.add_argument("--sample-attrs-dir", help="Directory with *_sample_attrs.csv")
    p_bld.add_argument("--out", help="Output CSV path")

    # --- validate ---
    p_val = subparsers.add_parser("validate", help="Phase 5: 7-check validation")
    p_val.add_argument("phenotype", help="Path to phenotype_clean.csv")
    p_val.add_argument("--registry", help="Path to intervention_registry.csv")
    p_val.add_argument("--out", help="Output JSON path")
    p_val.add_argument("--quiet", action="store_true", help="Silent pass/fail")

    # --- run ---
    p_run = subparsers.add_parser("run", help="Full pipeline (phases 1-3)")
    p_run.add_argument("--projects", help="Comma-separated bioproject accessions")
    p_run.add_argument("--all", action="store_true", help="Process all projects")
    p_run.add_argument("--sample-attrs", action="store_true", help="Strategy 2 download")
    p_run.add_argument("--resume", action="store_true", help="Resume from checkpoint")
    p_run.add_argument("--metadata", default="ena_metadata",
                       help="Metadata path (default: ena_metadata/)")
    p_run.add_argument("--out-dir", default="output", help="Output directory")

    # --- supplement ---
    p_sup = subparsers.add_parser("supplement", help="Utility: Parse supplementary tables")
    p_sup.add_argument("path", help="Path to supplement file (.xlsx/.csv/.tsv/.pdf)")
    p_sup.add_argument("--header-row", type=int, default=0,
                      help="Row index for headers (0-based)")
    p_sup.add_argument("--sheet", default="", help="Sheet name (xlsx only)")
    p_sup.add_argument("--table-index", type=int, default=0,
                      help="Table index (PDF only)")
    p_sup.add_argument("--skip-rows", type=int, default=0,
                      help="Rows to skip before header")
    p_sup.add_argument("--out", help="Output CSV path for normalized export")
    p_sup.add_argument("--report", action="store_true",
                      help="Print column detection report")
    p_sup.add_argument("--all-tables", action="store_true",
                      help="Parse all tables/sheets")
    p_sup.add_argument("--detect-units", action="store_true",
                      help="Auto-detect timepoint units (days vs weeks)")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    # Dispatch
    commands = {
        "download": cmd_download,
        "register": cmd_register,
        "classify": cmd_classify,
        "build": cmd_build,
        "validate": cmd_validate,
        "run": cmd_run,
        "supplement": cmd_supplement,
    }

    commands[args.command](args)


if __name__ == "__main__":
    main()