"""Write analysis results to a run directory: results.json and CSV tables."""

from __future__ import annotations

import csv
import json
from pathlib import Path

# CSV files written for each results key (when present), and the key to read.
CSV_TABLES = {
    "cycle_distribution.csv": "cycle_distribution",
    "cycle_count_distribution.csv": "cycle_count_distribution",
    "largest_cycle_distribution.csv": "largest_cycle_distribution",
    "null_cycle_distribution.csv": "null_cycle_distribution",
    "observed_cycle_structures.csv": "observed_cycle_structures_under_null",
    "rejoin_probabilities.csv": "rejoins",
}


def write_results(results: dict, directory: str | Path) -> list[Path]:
    """Write ``results.json`` and every applicable CSV table; return the paths
    written.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    written = [directory / "results.json"]
    with open(written[0], "w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    tables = dict(CSV_TABLES)
    completions = (results.get("observed") or {}).get("completions")
    if completions:
        rows = completions.get("distribution")
        if rows:
            written.append(_write_csv(directory / "completion_distribution.csv", rows))
    for file_name, key in tables.items():
        rows = results.get(key)
        if rows:
            written.append(_write_csv(directory / file_name, rows))
    return written


def flatten(row: dict, prefix: str = "") -> dict:
    """``{"a": {"b": 1}}`` → ``{"a_b": 1}``; lists become ``;``-joined strings."""
    out = {}
    for key, value in row.items():
        name = f"{prefix}{key}"
        if isinstance(value, dict):
            out.update(flatten(value, f"{name}_"))
        elif isinstance(value, list):
            out[name] = ";".join(map(str, value))
        else:
            out[name] = value
    return out


def _write_csv(path: Path, rows: list[dict]) -> Path:
    flat = [flatten(r) for r in rows]
    columns: list[str] = []
    for r in flat:
        columns += [c for c in r if c not in columns]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(flat)
    return path
