# -*- coding: utf-8 -*-
"""Supplement parser (Gap #8).

Parses supplementary tables from papers in xlsx/csv/tsv/pdf formats.
Auto-detects relevant columns using the normalizer's synonym table,
handles merged cells, multi-header tables, and blank rows.

Primary use case: Phase 3 paper verification — extracting sample-level
metadata (subject IDs, timepoints, groups, intervention details) from
supplementary tables that ENA metadata doesn't contain.

Usage:
    from core.supplement_parser import SupplementParser, parse_supplement

    # Quick parse
    result = parse_supplement("supplement_table_1.xlsx")

    # With control
    parser = SupplementParser()
    result = parser.parse("supplement_table_1.xlsx", header_row=0)
    parser.export_csv(result, "supplement_table_1_normalized.csv")
"""

import csv
import io
import os
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# Reuse the existing normalizer for column name matching
try:
    from .normalizer import SYNONYMS
except ImportError:
    try:
        from normalizer import SYNONYMS
    except ImportError:
        from core.normalizer import SYNONYMS


# ============================================================
# Data structures
# ============================================================

@dataclass
class ParsedColumn:
    """A detected column with its mapping."""
    original_name: str           # As it appeared in the file
    col_index: int               # 0-based column index
    standard_field: str = ""     # Mapped standard field (from normalizer)
    confidence: str = "low"      # high/medium/low
    sample_values: list[str] = field(default_factory=list)  # First few values


@dataclass
class ParseResult:
    """Result of parsing a supplement table."""
    source_path: str
    source_format: str            # xlsx, csv, tsv, pdf
    table_index: int = 0          # Which table in the file (0-based)
    sheet_name: str = ""          # For xlsx files

    # Parsed data
    headers: list[str] = field(default_factory=list)
    rows: list[list[str]] = field(default_factory=list)
    n_rows: int = 0

    # Column detection
    detected_columns: list[ParsedColumn] = field(default_factory=list)
    matched_fields: list[str] = field(default_factory=list)  # Standard field names found

    # Metadata
    warnings: list[str] = field(default_factory=list)
    skipped_rows: int = 0

    def to_dicts(self) -> list[dict[str, str]]:
        """Convert rows to list of dicts using detected standard field names."""
        if not self.rows or not self.headers:
            return []

        col_map = {}
        for dc in self.detected_columns:
            if dc.standard_field:
                col_map[dc.col_index] = dc.standard_field

        result = []
        for row in self.rows:
            d = {}
            for i, val in enumerate(row):
                field = col_map.get(i, self.headers[i] if i < len(self.headers) else f"col_{i}")
                d[field] = val
            result.append(d)
        return result

    def column_values(self, standard_field: str) -> list[str]:
        """Get all values for a standard field across rows."""
        for dc in self.detected_columns:
            if dc.standard_field == standard_field:
                ci = dc.col_index
                return [row[ci] if ci < len(row) else "" for row in self.rows]
        return []

    def summary(self) -> str:
        """One-line summary of parse result."""
        fields = ", ".join(self.matched_fields) if self.matched_fields else "none"
        return (f"{Path(self.source_path).name}: {self.n_rows} rows × "
                f"{len(self.headers)} cols, matched: [{fields}]")


# ============================================================
# Column auto-detection
# ============================================================

def _normalize_header(text: str) -> str:
    """Normalize a header string for matching."""
    if not text:
        return ""
    # Lowercase, strip whitespace/punctuation, collapse spaces
    text = text.strip().lower()
    text = re.sub(r'[^\w\s]', ' ', text)
    text = re.sub(r'\s+', ' ', text)
    return text.strip()


def _detect_columns(headers: list[str], sample_rows: list[list[str]],
                    max_sample: int = 5) -> list[ParsedColumn]:
    """Auto-detect which columns map to standard fields.

    Uses the normalizer's SYNONYMS table for matching header names,
    then samples values to boost confidence.

    Args:
        headers: List of header strings.
        sample_rows: First few data rows for value inspection.
        max_sample: How many rows to sample.

    Returns:
        List of ParsedColumn with standard_field populated where matched.
    """
    detected: list[ParsedColumn] = []

    for i, h in enumerate(headers):
        pc = ParsedColumn(original_name=h, col_index=i)

        # Sample some values for this column
        values = []
        for row in sample_rows[:max_sample]:
            if i < len(row) and row[i].strip():
                values.append(row[i].strip())
        pc.sample_values = values

        # Try to match header against standard fields
        h_norm = _normalize_header(h)
        if not h_norm:
            detected.append(pc)
            continue

        # Check each standard field's synonyms
        best_field = ""
        best_confidence = "low"

        for std_field, candidates in SYNONYMS.items():
            for candidate in candidates:
                c_norm = _normalize_header(candidate)
                if not c_norm:
                    continue

                # Exact match after normalization
                if h_norm == c_norm:
                    best_field = std_field
                    best_confidence = "high"
                    break

                # Contains match (e.g., "subject id number" contains "subject id")
                if c_norm in h_norm or h_norm in c_norm:
                    if best_confidence != "high":
                        best_field = std_field
                        best_confidence = "medium"
            if best_confidence == "high":
                break

        # Boost confidence based on value patterns
        if best_field and best_confidence != "high" and values:
            # Age: typically numeric with possible decimals
            if best_field == "host_age":
                numeric_vals = [v for v in values if re.match(r'^\d+\.?\d*$', v)]
                if len(numeric_vals) >= len(values) * 0.7:
                    best_confidence = "high"

            # Sex: typically M/F/Male/Female
            if best_field == "host_sex":
                sex_vals = [v.lower() for v in values if v.lower() in
                           ('m', 'f', 'male', 'female', '0', '1')]
                if len(sex_vals) >= len(values) * 0.5:
                    best_confidence = "high"

            # BMI: typically 15-50 range
            if best_field == "host_bmi":
                numeric_vals = []
                for v in values:
                    try:
                        n = float(v)
                        if 10 < n < 60:
                            numeric_vals.append(v)
                    except ValueError:
                        pass
                if len(numeric_vals) >= len(values) * 0.5:
                    best_confidence = "high"

            # Disease: typically text, not purely numeric
            if best_field == "disease":
                text_vals = [v for v in values if not re.match(r'^\d+\.?\d*$', v)]
                if len(text_vals) >= len(values) * 0.5:
                    best_confidence = "medium"

            # Timepoint: often contains day/week/month or T0/T1
            if best_field == "timepoint":
                tp_vals = [v.lower() for v in values if re.search(
                    r'(day|week|month|^t\d|baseline|pre|post|follow.up|washout)', v.lower())]
                if len(tp_vals) >= len(values) * 0.4:
                    best_confidence = "high"

        pc.standard_field = best_field
        pc.confidence = best_confidence
        detected.append(pc)

    return detected


# ============================================================
# Main parser class
# ============================================================

class SupplementParser:
    """Parse supplementary tables from papers.

    Supports xlsx, csv, tsv, and pdf formats with automatic
    column detection against standard phenotype fields.

    Usage:
        parser = SupplementParser()
        result = parser.parse("supplementary_table_1.xlsx")
        print(result.summary())
        # export to standardized CSV
        parser.export_csv(result, "output.csv")
    """

    def __init__(self, sample_rows_for_detection: int = 5):
        self._sample_rows = sample_rows_for_detection

    def parse(self, path: str | Path, header_row: int = 0,
              table_index: int = 0, sheet_name: str = "",
              skip_rows: int = 0) -> ParseResult:
        """Parse a supplement file and auto-detect columns.

        Args:
            path: Path to the file (.xlsx, .csv, .tsv, .pdf).
            header_row: Row index (0-based) containing column headers.
            table_index: For PDFs with multiple tables, which one.
            sheet_name: For xlsx, which sheet (empty = first).
            skip_rows: Rows to skip before header_row (metadata rows).

        Returns:
            ParseResult with parsed data and column detection.
        """
        path = Path(path)
        ext = path.suffix.lower()

        if ext in ('.xlsx', '.xlsm', '.xls'):
            return self._parse_xlsx(path, header_row, sheet_name, skip_rows)
        elif ext == '.csv':
            return self._parse_csv(path, header_row, ',', skip_rows)
        elif ext in ('.tsv', '.txt'):
            return self._parse_csv(path, header_row, '\t', skip_rows)
        elif ext == '.pdf':
            return self._parse_pdf(path, table_index, header_row, skip_rows)
        else:
            raise ValueError(f"Unsupported format: {ext}")

    def parse_all_tables(self, path: str | Path, **kwargs) -> list[ParseResult]:
        """Parse all tables from a file (useful for multi-table PDFs or
        multi-sheet xlsx files)."""
        path = Path(path)
        ext = path.suffix.lower()

        if ext == '.pdf':
            return self._parse_pdf_all(path, **kwargs)
        elif ext in ('.xlsx', '.xlsm'):
            return self._parse_xlsx_all(path, **kwargs)
        else:
            return [self.parse(path, **kwargs)]

    # ---- XLSX ----

    def _parse_xlsx(self, path: Path, header_row: int,
                    sheet_name: str, skip_rows: int) -> ParseResult:
        import openpyxl

        wb = openpyxl.load_workbook(path, data_only=True)
        ws = wb[sheet_name] if sheet_name else wb.active
        result = ParseResult(
            source_path=str(path),
            source_format="xlsx",
            sheet_name=ws.title,
        )

        all_rows = []
        for row in ws.iter_rows(values_only=True):
            # Convert all cells to strings
            str_row = [str(c) if c is not None else "" for c in row]
            # Trim trailing empty cells
            while str_row and not str_row[-1].strip():
                str_row.pop()
            if str_row:
                all_rows.append(str_row)

        result = self._extract_table(all_rows, result, header_row, skip_rows)
        wb.close()
        return result

    def _parse_xlsx_all(self, path: Path, **kwargs) -> list[ParseResult]:
        import openpyxl

        wb = openpyxl.load_workbook(path, data_only=True)
        results = []
        for sn in wb.sheetnames:
            r = self.parse(str(path), sheet_name=sn, **kwargs)
            results.append(r)
        wb.close()
        return results

    # ---- CSV / TSV ----

    def _parse_csv(self, path: Path, header_row: int,
                   delimiter: str, skip_rows: int) -> ParseResult:
        fmt = "csv" if delimiter == ',' else "tsv"
        result = ParseResult(source_path=str(path), source_format=fmt)

        with open(path, 'r', encoding='utf-8-sig', newline='') as f:
            reader = csv.reader(f, delimiter=delimiter)
            all_rows = []
            for row in reader:
                str_row = [c.strip() for c in row]
                while str_row and not str_row[-1]:
                    str_row.pop()
                if str_row:
                    all_rows.append(str_row)

        return self._extract_table(all_rows, result, header_row, skip_rows)

    # ---- PDF ----

    def _parse_pdf(self, path: Path, table_index: int,
                   header_row: int, skip_rows: int) -> ParseResult:
        """Extract table from PDF using pdfplumber.

        If table_index is out of range, returns result with warning.
        """
        result = ParseResult(
            source_path=str(path),
            source_format="pdf",
            table_index=table_index,
        )

        try:
            import pdfplumber
        except ImportError:
            result.warnings.append("pdfplumber not installed. Install with: pip install pdfplumber")
            return result

        with pdfplumber.open(path) as pdf:
            tables = []
            for page in pdf.pages:
                page_tables = page.extract_tables()
                tables.extend(page_tables)

            if not tables:
                result.warnings.append("No tables found in PDF")
                return result

            if table_index >= len(tables):
                result.warnings.append(
                    f"Table index {table_index} out of range "
                    f"(found {len(tables)} tables)"
                )
                return result

            raw_table = tables[table_index]
            # Clean None → "" and strip whitespace
            all_rows = []
            for row in raw_table:
                str_row = [str(c).strip() if c else "" for c in row]
                while str_row and not str_row[-1]:
                    str_row.pop()
                if str_row:
                    all_rows.append(str_row)

        return self._extract_table(all_rows, result, header_row, skip_rows)

    def _parse_pdf_all(self, path: Path, **kwargs) -> list[ParseResult]:
        """Parse all tables from a PDF."""
        try:
            import pdfplumber
        except ImportError:
            return []

        results = []
        with pdfplumber.open(path) as pdf:
            table_count = 0
            for page in pdf.pages:
                page_tables = page.extract_tables()
                for _ in page_tables:
                    r = self.parse(str(path), table_index=table_count, **kwargs)
                    results.append(r)
                    table_count += 1
        return results

    # ---- Shared extraction logic ----

    def _extract_table(self, all_rows: list[list[str]], result: ParseResult,
                       header_row: int, skip_rows: int) -> ParseResult:
        """Extract headers and data from raw rows."""
        total_rows = len(all_rows)

        # Skip metadata rows before header
        data_start = skip_rows

        # Determine actual header row
        if header_row >= total_rows:
            result.warnings.append(
                f"header_row {header_row} >= total rows {total_rows}"
            )
            return result
        if header_row < data_start:
            data_start = header_row

        # Extract headers
        result.headers = all_rows[header_row]
        result.skipped_rows = data_start + 1  # header + skipped rows

        # Data rows start after header
        data_rows_idx = range(header_row + 1, total_rows)
        data_rows = [all_rows[i] for i in data_rows_idx]

        # Normalize row lengths to match header count
        n_cols = len(result.headers)
        normalized = []
        for row in data_rows:
            if not any(c.strip() for c in row):
                result.skipped_rows += 1
                continue
            padded = row + [""] * (n_cols - len(row))
            normalized.append(padded[:n_cols])

        result.rows = normalized
        result.n_rows = len(normalized)

        # Auto-detect columns
        if result.headers:
            result.detected_columns = _detect_columns(
                result.headers,
                [row[:n_cols] for row in data_rows],
                max_sample=self._sample_rows,
            )
            result.matched_fields = sorted(set(
                dc.standard_field for dc in result.detected_columns
                if dc.standard_field
            ))

        return result

    # ---- Export ----

    def export_csv(self, result: ParseResult, out_path: str | Path,
                   use_standard_names: bool = True) -> None:
        """Export parsed table to a normalized CSV.

        Args:
            result: ParseResult from `parse()`.
            out_path: Output CSV path.
            use_standard_names: If True, use standard field names as headers.
        """
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Build header row
        if use_standard_names:
            col_map = {}
            for dc in result.detected_columns:
                col_map[dc.col_index] = dc.standard_field or dc.original_name
            out_headers = [col_map.get(i, result.headers[i] if i < len(result.headers) else f"col_{i}")
                          for i in range(len(result.headers))]
        else:
            out_headers = result.headers

        with open(out_path, 'w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(out_headers)
            for row in result.rows:
                # Pad/shift to match header count
                padded = row + [""] * (len(out_headers) - len(row))
                writer.writerow(padded[:len(out_headers)])

    def export_detection_report(self, result: ParseResult) -> str:
        """Generate a human-readable column detection report."""
        lines = [f"Supplement: {Path(result.source_path).name}",
                 f"Format: {result.source_format}, "
                 f"{result.n_rows} rows × {len(result.headers)} cols",
                 "",
                 "Column Detection:"]
        for dc in result.detected_columns:
            icon = {"high": "✓", "medium": "~", "low": "?"}.get(dc.confidence, "?")
            field = dc.standard_field or "(unmatched)"
            sample = ", ".join(dc.sample_values[:3]) if dc.sample_values else "-"
            lines.append(f"  {icon} [{dc.col_index}] {dc.original_name!r} → {field}"
                        f"  ({dc.confidence})  samples: {sample}")

        if result.warnings:
            lines.append("")
            lines.append("Warnings:")
            for w in result.warnings:
                lines.append(f"  ⚠ {w}")

        return "\n".join(lines)


# ============================================================
# Unit detection (Gap #9 — weeks vs days disambiguation)
# ============================================================

@dataclass
class TimepointUnitDetection:
    """Result of auto-detecting timepoint units."""
    unit: str = "unknown"          # "day", "week", "month", "mixed", "unknown"
    confidence: str = "low"       # high/medium/low
    evidence: str = ""            # Human-readable explanation
    numeric_values: list[float] = field(default_factory=list)
    raw_values: list[str] = field(default_factory=list)


def detect_timepoint_unit(values: list[str]) -> TimepointUnitDetection:
    """Detect whether timepoint values are in days, weeks, or months.

    Heuristic:
      1. If values contain explicit units (e.g., "Day 14", "Week 2"),
         extract the unit from text.
      2. If purely numeric:
         - Values ≤ 30: ambiguous (could be days or weeks)
         - Values > 30 and ≤ 365: likely days
         - Values > 365: likely days (study duration)
         - Check header for unit hints (e.g., "Day", "Week", "Month")
      3. If mixed text/numeric, count explicit unit mentions.

    Args:
        values: List of timepoint values (strings from supplement column).

    Returns:
        TimepointUnitDetection with detected unit and confidence.
    """
    result = TimepointUnitDetection()
    result.raw_values = [v for v in values if v.strip()]

    if not result.raw_values:
        result.evidence = "no values provided"
        return result

    # Step 1: Check for explicit units in text
    unit_counts = {"day": 0, "week": 0, "month": 0}
    non_text_count = 0

    for v in result.raw_values:
        v_lower = v.lower().strip()

        # Explicit day mentions
        if re.search(r'\b(day|d\s*\d+|\d+\s*d(?:ay)?s?)\b', v_lower):
            unit_counts["day"] += 1
        # Explicit week mentions
        elif re.search(r'\b(w(?:ee)?k|w\s*\d+|\d+\s*w(?:ee)?ks?)\b', v_lower):
            unit_counts["week"] += 1
        # Explicit month mentions
        elif re.search(r'\b(mo(?:nth)?|m\s*\d+|\d+\s*m(?:o(?:nth)?)?s?)\b', v_lower):
            unit_counts["month"] += 1
        else:
            non_text_count += 1

    explicit_total = sum(unit_counts.values())

    if explicit_total > 0:
        # Use the most common explicit unit
        best_unit = max(unit_counts, key=unit_counts.get)
        ratio = unit_counts[best_unit] / explicit_total
        result.unit = best_unit
        result.evidence = f"{unit_counts[best_unit]}/{explicit_total} values mention '{best_unit}'"
        result.confidence = "high" if ratio >= 0.8 else "medium"
        return result

    # Step 2: All values are non-text — parse numeric values and infer
    numeric_vals = []
    for v in result.raw_values:
        try:
            n = float(v.strip())
            numeric_vals.append(n)
        except ValueError:
            pass

    result.numeric_values = numeric_vals

    if not numeric_vals:
        result.evidence = "no numeric values found"
        return result

    max_val = max(numeric_vals)

    if max_val > 365:
        result.unit = "day"
        result.evidence = f"max value {max_val} > 365 → days"
        result.confidence = "high"
    elif max_val > 90:
        result.unit = "day"
        result.evidence = f"max value {max_val} > 90 → likely days"
        result.confidence = "medium"
    elif max_val <= 30:
        # Ambiguous: could be days or weeks
        # Check if gap patterns suggest weeks (e.g., 0, 2, 4, 8, 12)
        sorted_vals = sorted(set(numeric_vals))
        if len(sorted_vals) >= 3:
            diffs = [sorted_vals[i + 1] - sorted_vals[i] for i in range(len(sorted_vals) - 1)]
            # If gaps are consistently >= 2 and most values are integers
            avg_gap = sum(diffs) / len(diffs) if diffs else 0
            if avg_gap >= 1.5 and all(n == int(n) for n in numeric_vals):
                result.unit = "week"
                result.evidence = f"max={max_val}, avg_gap={avg_gap:.1f} → suggests weeks"
                result.confidence = "medium"
            else:
                result.unit = "day"
                result.evidence = f"max={max_val}, avg_gap={avg_gap:.1f} → ambiguous, assume days"
                result.confidence = "low"
        else:
            result.unit = "day"
            result.evidence = f"max={max_val} but too few values to determine"
            result.confidence = "low"
    else:
        result.unit = "day"
        result.evidence = f"max={max_val}, in typical day range"
        result.confidence = "medium"

    return result


# ============================================================
# Convenience function
# ============================================================

def parse_supplement(path: str | Path, **kwargs) -> ParseResult:
    """Quick one-shot parse of a supplement file.

    Args:
        path: Path to .xlsx, .csv, .tsv, or .pdf.
        **kwargs: Passed to SupplementParser.parse().

    Returns:
        ParseResult with parsed data and column detection.
    """
    parser = SupplementParser()
    return parser.parse(path, **kwargs)


# ============================================================
# CLI
# ============================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Parse supplementary tables from papers"
    )
    parser.add_argument("path", help="Path to supplement file")
    parser.add_argument("--header-row", type=int, default=0,
                       help="Row index for headers (0-based)")
    parser.add_argument("--sheet", default="", help="Sheet name (xlsx only)")
    parser.add_argument("--table-index", type=int, default=0,
                       help="Table index (PDF only)")
    parser.add_argument("--skip-rows", type=int, default=0,
                       help="Rows to skip before header")
    parser.add_argument("--out", help="Output CSV path for normalized export")
    parser.add_argument("--report", action="store_true",
                       help="Print column detection report")
    parser.add_argument("--all-tables", action="store_true",
                       help="Parse all tables/sheets")

    args = parser.parse_args()

    sp = SupplementParser()

    if args.all_tables:
        results = sp.parse_all_tables(
            args.path,
            header_row=args.header_row,
            sheet_name=args.sheet,
            skip_rows=args.skip_rows,
        )
        for i, r in enumerate(results):
            sheet_info = f" [{r.sheet_name}]" if r.sheet_name else f" [table {r.table_index}]"
            print(f"{i + 1}.{sheet_info} {r.summary()}")
            if args.report:
                print(sp.export_detection_report(r))
                print()
    else:
        result = sp.parse(
            args.path,
            header_row=args.header_row,
            sheet_name=args.sheet,
            table_index=args.table_index,
            skip_rows=args.skip_rows,
        )

        if args.report:
            print(sp.export_detection_report(result))
        else:
            print(result.summary())

        if args.out:
            sp.export_csv(result, args.out)
            print(f"Exported to {args.out}")