# -*- coding: utf-8 -*-
"""Schema validation helpers for phenotype curation pipeline.

Loads JSON schemas from schemas/ and validates data at each phase boundary.
Uses jsonschema (draft-07) for structural validation and adds row-level
error collection with source context.
"""

import json
import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import jsonschema

_SCHEMA_DIR = Path(__file__).resolve().parent.parent / "schemas"
_SCHEMA_CACHE: dict[str, dict] = {}


# ============================================================
# Schema Loading
# ============================================================

def _resolve_schema(uri: str) -> dict:
    """Resolve a $id to its loaded schema (for cross-schema $ref)."""
    for schema in _SCHEMA_CACHE.values():
        if schema.get("$id") == uri:
            return schema
    raise jsonschema.RefResolutionError(f"Schema $id={uri} not found in cache")


def load_schema(name: str) -> dict:
    """Load a JSON schema by name (without .schema.json).

    Args:
        name: One of 'registry_row', 'intervention_registry_row',
              'phenotype_row', 'traceability_row'.

    Returns:
        Parsed schema dict.
    """
    if name in _SCHEMA_CACHE:
        return _SCHEMA_CACHE[name]

    path = _SCHEMA_DIR / f"{name}.schema.json"
    if not path.exists():
        raise FileNotFoundError(f"Schema not found: {path}")

    with open(path, "r", encoding="utf-8") as f:
        schema = json.load(f)

    _SCHEMA_CACHE[name] = schema

    # Also index by $id for cross-schema resolution
    if "$id" in schema:
        _SCHEMA_CACHE[schema["$id"]] = schema

    return schema


def preload_all_schemas() -> list[str]:
    """Load all schemas into cache. Returns list of loaded names."""
    names = []
    for path in sorted(_SCHEMA_DIR.glob("*.schema.json")):
        name = path.stem.replace(".schema", "")
        load_schema(name)
        names.append(name)
    return names


# ============================================================
# Validation Result Types
# ============================================================

@dataclass
class ValidationError:
    """A single validation failure."""
    row_index: int          # 0-based row index in CSV
    field: str = ""         # problematic field (empty for row-level errors)
    message: str = ""       # human-readable error description
    value: str = ""         # actual value that failed
    json_path: str = ""     # JSON path to the error (e.g. "$.study_design")

    def __str__(self) -> str:
        loc = f"row {self.row_index + 1}"
        if self.field:
            loc += f", field '{self.field}'"
        return f"[{loc}] {self.message} (got: {self.value!r})"


@dataclass
class ValidationResult:
    """Result of validating a CSV against a schema."""
    schema_name: str = ""
    total_rows: int = 0
    valid_rows: int = 0
    errors: list[ValidationError] = field(default_factory=list)
    warnings: list[ValidationError] = field(default_factory=list)

    @property
    def error_count(self) -> int:
        return len(self.errors)

    @property
    def passed(self) -> bool:
        return self.error_count == 0

    @property
    def pass_rate(self) -> float:
        if self.total_rows == 0:
            return 0.0
        return self.valid_rows / self.total_rows

    def summary(self) -> str:
        lines = [
            f"Schema: {self.schema_name}",
            f"Rows: {self.total_rows} total, {self.valid_rows} valid, "
            f"{self.error_count} errors, {len(self.warnings)} warnings",
            f"Pass rate: {self.pass_rate:.1%}",
        ]
        if self.errors:
            lines.append("Errors:")
            for e in self.errors[:20]:  # cap at 20 for readability
                lines.append(f"  {e}")
            if len(self.errors) > 20:
                lines.append(f"  ... and {len(self.errors) - 20} more")
        if self.warnings:
            lines.append("Warnings:")
            for w in self.warnings[:10]:
                lines.append(f"  {w}")
        return "\n".join(lines)


# ============================================================
# Row-Level Validation
# ============================================================

def _make_validator(schema: dict) -> jsonschema.Draft7Validator:
    """Create a Draft7Validator with local schema resolution."""
    resolver = jsonschema.RefResolver(
        base_uri="",
        referrer=schema,
        handlers={"": lambda uri: _resolve_schema(uri.strip("/"))},
    )
    return jsonschema.Draft7Validator(schema, resolver=resolver)


def validate_row(row: dict, schema_name: str) -> list[ValidationError]:
    """Validate a single dict row against a named schema.

    Args:
        row: Dict of field_name -> value.
        schema_name: Name of schema to validate against.

    Returns:
        List of ValidationError objects (empty = valid).
    """
    schema = load_schema(schema_name)
    validator = _make_validator(schema)
    errors: list[ValidationError] = []

    for js_error in validator.iter_errors(row):
        # Extract the field name from the JSON path
        field = ""
        json_path = ""
        if js_error.absolute_path:
            json_path = "$." + ".".join(str(p) for p in js_error.absolute_path)
            field = str(js_error.absolute_path[-1]) if js_error.absolute_path else ""

        value = js_error.instance if js_error.instance is not None else ""

        errors.append(ValidationError(
            row_index=-1,   # caller should set this
            field=field,
            message=js_error.message,
            value=str(value)[:200],
            json_path=json_path,
        ))

    # Additional check: warn on unknown fields not in schema
    known_fields = set(schema.get("properties", {}).keys())
    extra_fields = set(row.keys()) - known_fields
    for ef in sorted(extra_fields):
        if ef:  # skip empty keys
            errors.append(ValidationError(
                row_index=-1,
                field=ef,
                message=f"Unknown field (not in schema {schema_name})",
                value=str(row.get(ef, ""))[:200],
            ))

    return errors


# ============================================================
# CSV/Batch Validation
# ============================================================

def validate_csv(
    csv_path: str | Path,
    schema_name: str,
    encoding: str = "utf-8",
    max_errors: int = 200,
) -> ValidationResult:
    """Validate all rows in a CSV file against a schema.

    Args:
        csv_path: Path to CSV file.
        schema_name: Name of schema to validate against.
        encoding: File encoding (default utf-8).
        max_errors: Stop after this many errors (0 = no limit).

    Returns:
        ValidationResult with all errors collected.
    """
    schema = load_schema(schema_name)
    result = ValidationResult(schema_name=schema_name)

    with open(csv_path, "r", encoding=encoding, newline="") as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            result.total_rows += 1
            row_errors = validate_row(dict(row), schema_name)
            if row_errors:
                for e in row_errors:
                    e.row_index = i
                result.errors.extend(row_errors)
            else:
                result.valid_rows += 1

            if max_errors and len(result.errors) >= max_errors:
                break

    return result


def validate_rows(
    rows: list[dict],
    schema_name: str,
    max_errors: int = 200,
) -> ValidationResult:
    """Validate a list of dict rows against a schema.

    Args:
        rows: List of dicts to validate.
        schema_name: Name of schema to validate against.
        max_errors: Stop after this many errors (0 = no limit).

    Returns:
        ValidationResult with all errors collected.
    """
    result = ValidationResult(schema_name=schema_name)
    result.total_rows = len(rows)

    for i, row in enumerate(rows):
        row_errors = validate_row(dict(row), schema_name)
        if row_errors:
            for e in row_errors:
                e.row_index = i
            result.errors.extend(row_errors)
        else:
            result.valid_rows += 1

        if max_errors and len(result.errors) >= max_errors:
            break

    return result


# ============================================================
# Convenience Functions
# ============================================================

def validate_registry(csv_path: str | Path) -> ValidationResult:
    """Validate a project_registry.csv."""
    return validate_csv(csv_path, "registry_row")


def validate_intervention_registry(csv_path: str | Path) -> ValidationResult:
    """Validate a project_intervention_registry.csv."""
    return validate_csv(csv_path, "intervention_registry_row")


def validate_phenotype(csv_path: str | Path) -> ValidationResult:
    """Validate a phenotype_clean.csv."""
    return validate_csv(csv_path, "phenotype_row")


def validate_traceability(csv_path: str | Path) -> ValidationResult:
    """Validate an is_intervention_traceability.csv."""
    return validate_csv(csv_path, "traceability_row")


# Preload on import for fast first-use
preload_all_schemas()