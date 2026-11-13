# -*- coding: utf-8 -*-
"""Phase 5: 7-check automated validation of phenotype table.

Input:
  - phenotype_clean.csv (Phase 4 output)
  - intervention_registry.csv (Phase 3 output)

Output:
  - validation_report.json (structured results)

Checks:
  1. Row count consistency (phenotype vs expected from registry)
  2. Registry consistency — observational studies: all is_intervention=no
  3. One-time interventions: all samples is_intervention=no
  4. Ongoing intervention RCTs: must have at least some yes samples
  5. Phase label ↔ is_intervention consistency
  6. No duplicate biosamples
  7. Value range sanity (age, BMI)
  + Deal-breaker: observational study with is_intervention=yes
"""

import csv
import json
import re
import sys
from collections import Counter, defaultdict
from dataclasses import asdict, field
from pathlib import Path
from typing import Optional

_PARENT = Path(__file__).resolve().parent.parent
if str(_PARENT) not in sys.path:
    sys.path.insert(0, str(_PARENT))

from core.types import ValidationReport


# ============================================================
# Check Execution
# ============================================================

def validate(
    phenotype_path: str | Path,
    registry_path: Optional[str | Path] = None,
    out_path: Optional[str | Path] = None,
) -> ValidationReport:
    """Run all 7 validation checks.

    Args:
        phenotype_path: Path to phenotype_clean.csv.
        registry_path: Path to intervention_registry.csv (optional, enables checks 1-4).
        out_path: Output JSON path (default: validation_report.json).

    Returns:
        ValidationReport dataclass.
    """
    # Load data
    print(f"Loading phenotype from {phenotype_path}...")
    with open(phenotype_path, "r", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    print(f"  Loaded {len(rows)} rows")

    registry = {}
    if registry_path:
        print(f"Loading registry from {registry_path}...")
        with open(registry_path, "r", encoding="utf-8-sig") as f:
            registry = {r["bioproject"]: r for r in csv.DictReader(f)}
        print(f"  Loaded {len(registry)} registrations")

    # Run checks
    checks = []
    errors = []
    warnings = []

    # CHECK 1: Row count
    print("\n[CHECK 1] Row count consistency")
    c1 = _check_row_count(rows, registry)
    checks.append(c1)
    if not c1["passed"]:
        errors.append(c1["details"])

    # CHECK 2: Observational studies — all no
    if registry:
        print("\n[CHECK 2] Observational studies: all is_intervention=no?")
        c2 = _check_observational_all_no(rows, registry)
        checks.append(c2)
        if not c2["passed"]:
            errors.append(c2["details"])

    # CHECK 3: One-time interventions — all no
    if registry:
        print("\n[CHECK 3] One-time interventions: all samples is_intervention=no?")
        c3 = _check_onetime_all_no(rows, registry)
        checks.append(c3)
        if not c3["passed"]:
            errors.append(c3["details"])

    # CHECK 4: Ongoing interventions — has yes samples
    if registry:
        print("\n[CHECK 4] Ongoing-intervention RCTs: must have is_intervention=yes?")
        c4 = _check_ongoing_has_yes(rows, registry)
        checks.append(c4)
        if not c4["passed"]:
            warnings.append(c4["details"])

    # CHECK 5: Phase label ↔ is_intervention consistency
    print("\n[CHECK 5] Phase label ↔ is_intervention consistency")
    c5 = _check_phase_label_consistency(rows)
    checks.append(c5)
    if not c5["passed"]:
        errors.append(c5["details"])

    # CHECK 6: No duplicate biosamples
    print("\n[CHECK 6] Duplicate biosample check")
    c6 = _check_duplicate_biosamples(rows)
    checks.append(c6)
    if not c6["passed"]:
        errors.append(c6["details"])

    # CHECK 7: Value range sanity
    print("\n[CHECK 7] Value range sanity")
    c7 = _check_value_ranges(rows)
    checks.append(c7)
    if not c7["passed"]:
        warnings.append(c7["details"])

    # Deal-breaker
    if registry:
        dealbreaker = _check_dealbreaker(rows, registry)
        if dealbreaker:
            errors.insert(0, dealbreaker)

    # is_intervention distribution
    is_intv_dist = dict(Counter(r["is_intervention"] for r in rows))

    # Build report
    report = ValidationReport(
        total_rows=len(rows),
        expected_rows=sum(
            int(r.get("biosample_n_local", 0))
            for r in registry.values()
            if r.get("biosample_n_local", "").isdigit()
        ) if registry else len(rows),
        checks=checks,
        is_intervention_distribution=is_intv_dist,
        errors=errors,
        warnings=warnings,
        passed=len(errors) == 0,
    )

    # Print summary
    print(f"\n{'='*60}")
    print("VALIDATION SUMMARY")
    print(f"{'='*60}")
    print(f"  Total rows: {report.total_rows}")
    print(f"  is_intervention: {is_intv_dist}")
    print(f"  Checks passed: {sum(1 for c in checks if c['passed'])}/{len(checks)}")
    if errors:
        print(f"  ERRORS ({len(errors)}):")
        for e in errors[:5]:
            print(f"    {e}")
    if warnings:
        print(f"  WARNINGS ({len(warnings)}):")
        for w in warnings[:5]:
            print(f"    {w}")
    print(f"  OVERALL: {'PASS' if report.passed else 'FAIL'}")

    # Write output
    if out_path:
        out_path = Path(out_path)
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(_report_to_dict(report), f, indent=2, ensure_ascii=False)
        print(f"\nReport written to {out_path}")

    return report


# ============================================================
# Individual Checks
# ============================================================

def _check_row_count(rows: list[dict], registry: dict) -> dict:
    """Check 1: Row count matches expected from registry."""
    total = len(rows)
    # Expected from registry biosample counts
    expected = 0
    for bp, r in registry.items():
        n = r.get("biosample_n_local", "")
        if n.isdigit():
            expected += int(n)

    details = f"{total} rows, registry expects {expected}"
    if expected and abs(total - expected) > total * 0.1:
        return _check_result("row_count", False, details)
    return _check_result("row_count", True, details)


def _check_observational_all_no(rows: list[dict], registry: dict) -> dict:
    """Check 2: Observational studies must have ALL is_intervention=no."""
    observational_bps = {
        bp for bp, r in registry.items()
        if r.get("intervention_type") == "observational"
    }
    failures = []
    for bp in sorted(observational_bps):
        bpr = [r for r in rows if r["bioproject"] == bp]
        if not bpr:
            continue
        bad = Counter(r["is_intervention"] for r in bpr if r["is_intervention"] != "no")
        if bad:
            failures.append(f"{bp}: {dict(bad)}")

    if failures:
        return _check_result("observational_all_no", False,
                             f"{len(failures)} projects with non-no: {'; '.join(failures[:5])}")
    return _check_result("observational_all_no", True,
                         f"All {len(observational_bps)} observational projects pass")


def _check_onetime_all_no(rows: list[dict], registry: dict) -> dict:
    """Check 3: One-time intervention studies must have ALL is_intervention=no."""
    onetime_bps = {
        bp for bp, r in registry.items()
        if r.get("is_one_time") == "yes"
    }
    failures = []
    for bp in sorted(onetime_bps):
        bpr = [r for r in rows if r["bioproject"] == bp]
        if not bpr:
            continue
        yes_count = sum(1 for r in bpr if r["is_intervention"] == "yes")
        if yes_count > 0:
            failures.append(f"{bp}: {yes_count} yes samples")

    if failures:
        return _check_result("onetime_all_no", False,
                             f"{len(failures)} one-time projects have yes: {'; '.join(failures)}")
    return _check_result("onetime_all_no", True,
                         f"All {len(onetime_bps)} one-time projects pass")


def _check_ongoing_has_yes(rows: list[dict], registry: dict) -> dict:
    """Check 4: Ongoing intervention RCTs should have at least some yes samples."""
    ongoing_bps = {
        bp for bp, r in registry.items()
        if (r.get("intervention_type") != "observational"
            and r.get("is_one_time") != "yes"
            and r.get("study_design") in ("RCT", "RCT_crossover", "dose_escalation"))
    }
    warnings_list = []
    for bp in sorted(ongoing_bps):
        bpr = [r for r in rows if r["bioproject"] == bp]
        if not bpr:
            continue
        has_yes = any(r["is_intervention"] == "yes" for r in bpr)
        if not has_yes:
            warnings_list.append(f"{bp}: 0 yes samples (intervention window may be missed)")

    if warnings_list:
        return _check_result("ongoing_has_yes", False,
                             f"{len(warnings_list)} ongoing RCTs with zero yes: {'; '.join(warnings_list[:5])}")
    return _check_result("ongoing_has_yes", True,
                         f"All {len(ongoing_bps)} ongoing RCTs have yes samples")


def _check_phase_label_consistency(rows: list[dict]) -> dict:
    """Check 5: Phase labels must be consistent with is_intervention."""
    INTERVENTION_LABELS = (
        "feeding", "on_treatment", "during_", "fmt_week",
        "probiotic_feeding",
    )
    NON_INTERVENTION_LABELS = (
        "baseline", "run_in", "washout", "followup", "post_",
        "peri_transplant", "post_treatment",
    )

    errors = []
    for r in rows:
        label = r.get("timepoint_label", "")
        is_intv = r.get("is_intervention", "")
        bp = r.get("bioproject", "?")
        bs = r.get("biosample_id", "?")

        # Intervention labels → is_intervention should be yes
        if any(p in label for p in INTERVENTION_LABELS) and is_intv != "yes":
            errors.append(f"{bp}/{bs}: {label} → is_intervention={is_intv} (expected yes)")

        # Non-intervention labels → is_intervention should not be yes
        if any(p in label for p in NON_INTERVENTION_LABELS) and is_intv == "yes":
            errors.append(f"{bp}/{bs}: {label} → is_intervention=yes (expected no)")

    if errors:
        return _check_result("phase_label_consistency", False,
                             f"{len(errors)} inconsistencies")

    return _check_result("phase_label_consistency", True,
                         "All phase labels consistent with is_intervention")


def _check_duplicate_biosamples(rows: list[dict]) -> dict:
    """Check 6: No duplicate biosample_ids."""
    seen = defaultdict(list)
    for i, r in enumerate(rows):
        bs = r.get("biosample_id", "")
        if bs:
            seen[bs].append(i)

    duplicates = {bs: indices for bs, indices in seen.items() if len(indices) > 1}
    if duplicates:
        examples = list(duplicates.items())[:3]
        return _check_result("duplicate_biosamples", False,
                             f"{len(duplicates)} duplicate biosamples: {examples}")

    return _check_result("duplicate_biosamples", True,
                         f"No duplicates among {len(seen)} unique biosamples")


def _check_value_ranges(rows: list[dict]) -> dict:
    """Check 7: Age and BMI within reasonable ranges."""
    age_issues = []
    bmi_issues = []

    for r in rows:
        age_str = r.get("host_age", "").strip()
        if age_str:
            try:
                age = float(age_str)
                if age < 0 or age > 120:
                    age_issues.append(f"{r.get('biosample_id','?')}: age={age}")
            except ValueError:
                # Non-numeric age is acceptable (e.g. "60-69")
                pass

        bmi_str = r.get("host_bmi", "").strip()
        if bmi_str:
            try:
                bmi = float(bmi_str)
                if bmi < 10 or bmi > 80:
                    bmi_issues.append(f"{r.get('biosample_id','?')}: BMI={bmi}")
            except ValueError:
                pass

    warnings_list = []
    if age_issues:
        warnings_list.append(f"{len(age_issues)} age outliers")
    if bmi_issues:
        warnings_list.append(f"{len(bmi_issues)} BMI outliers")

    if warnings_list:
        return _check_result("value_ranges", False, "; ".join(warnings_list))

    return _check_result("value_ranges", True,
                         "Age and BMI within reasonable ranges")


def _check_dealbreaker(rows: list[dict], registry: dict) -> Optional[str]:
    """Deal-breaker: observational study with any is_intervention=yes."""
    for r in rows:
        bp = r.get("bioproject", "")
        rinfo = registry.get(bp, {})
        if (rinfo.get("intervention_type") == "observational"
                and r.get("is_intervention") == "yes"):
            return (f"DEAL-BREAKER: {bp} observational but "
                    f"is_intervention=yes at {r.get('biosample_id', '?')}")
    return None


# ============================================================
# Helpers
# ============================================================

def _check_result(name: str, passed: bool, details: str) -> dict:
    """Format a single check result."""
    status = "PASS" if passed else "FAIL"
    print(f"  [{status}] {details}")
    return {"name": name, "passed": passed, "details": details}


def _report_to_dict(report: ValidationReport) -> dict:
    """Convert ValidationReport to JSON-serializable dict."""
    return {
        "total_rows": report.total_rows,
        "expected_rows": report.expected_rows,
        "passed": report.passed,
        "is_intervention_distribution": report.is_intervention_distribution,
        "checks": report.checks,
        "errors": report.errors,
        "warnings": report.warnings,
    }


# ============================================================
# Convenience
# ============================================================

def quick_check(phenotype_path: str | Path, registry_path: Optional[str | Path] = None) -> bool:
    """Quick pass/fail check. Returns True if all non-warning checks pass."""
    report = validate(phenotype_path, registry_path)
    return report.passed


# ============================================================
# CLI
# ============================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Phase 5: 7-check validation of phenotype table"
    )
    parser.add_argument("phenotype", help="Path to phenotype_clean.csv")
    parser.add_argument(
        "--registry",
        help="Path to intervention_registry.csv (enables registry checks)",
        default=None,
    )
    parser.add_argument(
        "--out",
        default="validation_report.json",
        help="Output JSON path",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Only print pass/fail summary",
    )

    args = parser.parse_args()
    report = validate(
        phenotype_path=args.phenotype,
        registry_path=args.registry,
        out_path=args.out if not args.quiet else None,
    )

    sys.exit(0 if report.passed else 1)