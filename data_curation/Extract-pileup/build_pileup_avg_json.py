#!/usr/bin/env python3
"""Build per-position pileup summary (median + IQR) for a virus folder.

Usage:
    python Extract-pileup/build_pileup_avg_json.py cov-229e

This reads all pileup depth files in:
    Extract-pileup/all_pileups/<virus_name>/
and writes:
    Extract-pileup/all_pileups/<virus_name>/<virus_name>_avg.json
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List, Tuple


def _is_number(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def _parse_depth_line(line: str) -> Tuple[str | None, int, float] | None:
    """Parse one pileup/depth line and return (label, position, depth)."""
    line = line.strip()
    if not line:
        return None

    parts = line.split()
    if len(parts) < 2:
        return None

    # Common formats:
    # 1) ref pos depth
    # 2) pos depth
    if len(parts) >= 3 and _is_number(parts[1]) and _is_number(parts[2]):
        return parts[0], int(float(parts[1])), float(parts[2])
    if _is_number(parts[0]) and _is_number(parts[1]):
        return None, int(float(parts[0])), float(parts[1])

    return None


def _iter_depth_file(path: Path):
    """Yield normalized (label, position, depth) tuples from one depth file."""
    with path.open("r", encoding="utf-8") as fh:
        for raw_line in fh:
            parsed = _parse_depth_line(raw_line)
            if parsed is None:
                continue
            label, pos, depth = parsed
            if pos <= 0:
                continue
            yield label, pos, depth


def _load_segment_ref_map(
    segment_ref_table: Path, virus_name: str
) -> tuple[Dict[str, str], set[str]]:
    """Map reference/segment labels to canonical segment names for one virus."""
    ref_to_segment: Dict[str, str] = {}
    segment_names: set[str] = set()
    if not segment_ref_table.exists():
        return ref_to_segment, segment_names

    with segment_ref_table.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh)
        for row in reader:
            if len(row) < 3:
                continue
            ref = row[0].strip()
            virus = row[1].strip()
            segment = row[2].strip()
            if not ref or not virus or not segment:
                continue
            if virus != virus_name:
                continue
            ref_to_segment[ref] = segment
            segment_names.add(segment)
    return ref_to_segment, segment_names


def _percentile_sorted(values: List[float], p: float) -> float:
    """Linear-interpolated percentile for sorted values, p in [0,1]."""
    if not values:
        return float("nan")
    if len(values) == 1:
        return values[0]

    idx = (len(values) - 1) * p
    lower = int(idx)
    upper = min(lower + 1, len(values) - 1)
    frac = idx - lower
    return values[lower] * (1.0 - frac) + values[upper] * frac


def _collect_input_files(folder: Path, virus_name: str) -> List[Path]:
    files = [
        p
        for p in sorted(folder.iterdir())
        if p.is_file()
        and not p.name.startswith(".")
        and p.suffix.lower() in {".txt", ".tsv", ".csv"}
        and p.name != f"{virus_name}_avg.json"
    ]
    return files


def build_avg_json(virus_name: str, input_root: Path, segment_ref_table: Path) -> Path:
    folder = input_root / virus_name
    if not folder.exists() or not folder.is_dir():
        raise FileNotFoundError(f"Virus folder not found: {folder}")

    input_files = _collect_input_files(folder, virus_name)
    if not input_files:
        raise FileNotFoundError(f"No pileup input files found in: {folder}")

    by_position: Dict[int, List[float]] = {}
    by_segment_position: Dict[str, Dict[int, List[float]]] = {}
    segment_file_counts: Dict[str, int] = {}
    used_files = 0
    unknown_ref_counts: Dict[str, int] = {}

    ref_to_segment, segment_names = _load_segment_ref_map(segment_ref_table, virus_name)
    segment_map_active = bool(ref_to_segment or segment_names)

    for path in input_files:
        file_has_any = False
        file_segments_seen: set[str] = set()
        for label, pos, depth in _iter_depth_file(path):
            if segment_map_active:
                if label is None:
                    continue
                if label in ref_to_segment:
                    segment = ref_to_segment[label]
                elif label in segment_names:
                    segment = label
                else:
                    unknown_ref_counts[label] = unknown_ref_counts.get(label, 0) + 1
                    continue

                by_segment_position.setdefault(segment, {}).setdefault(pos, []).append(depth)
                file_segments_seen.add(segment)
                file_has_any = True
            else:
                by_position.setdefault(pos, []).append(depth)
                file_has_any = True

        if not file_has_any:
            continue
        used_files += 1
        for segment in file_segments_seen:
            segment_file_counts[segment] = segment_file_counts.get(segment, 0) + 1

    if unknown_ref_counts:
        for ref_name, count in sorted(unknown_ref_counts.items()):
            print(f"[build_pileup_avg_json] Skipping unknown segment/reference '{ref_name}' ({count} rows)")

    if used_files == 0 or (not by_position and not by_segment_position):
        raise ValueError(f"No valid depth values parsed from files in: {folder}")

    if by_segment_position:
        segments_out: Dict[str, dict] = {}
        for segment in sorted(by_segment_position.keys()):
            positions = sorted(by_segment_position[segment].keys())
            median: List[float] = []
            q1: List[float] = []
            q3: List[float] = []
            counts: List[int] = []

            for pos in positions:
                vals = sorted(by_segment_position[segment][pos])
                median.append(_percentile_sorted(vals, 0.5))
                q1.append(_percentile_sorted(vals, 0.25))
                q3.append(_percentile_sorted(vals, 0.75))
                counts.append(len(vals))

            segments_out[segment] = {
                "positions": positions,
                "median": median,
                "q1": q1,
                "q3": q3,
                "n_samples_per_position": counts,
                "n_files": segment_file_counts.get(segment, 0),
            }

        output = {
            "virus": virus_name,
            "is_segmented": True,
            "segments": segments_out,
        }
    else:
        positions = sorted(by_position.keys())
        median = []
        q1 = []
        q3 = []
        counts = []

        for pos in positions:
            vals = sorted(by_position[pos])
            median.append(_percentile_sorted(vals, 0.5))
            q1.append(_percentile_sorted(vals, 0.25))
            q3.append(_percentile_sorted(vals, 0.75))
            counts.append(len(vals))

        output = {
            "virus": virus_name,
            "n_files": used_files,
            "positions": positions,
            "median": median,
            "q1": q1,
            "q3": q3,
            "n_samples_per_position": counts,
        }

    out_path = folder / f"{virus_name}_avg.json"
    with out_path.open("w", encoding="utf-8") as fh:
        json.dump(output, fh, ensure_ascii=True)

    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-position pileup median and IQR for all depth files in "
            "Extract-pileup/all_pileups/<virus_name>."
        )
    )
    parser.add_argument("virus_name", help="Virus folder name under all_pileups/")
    parser.add_argument(
        "--input-root",
        default="all_pileups",
        help="Root directory containing virus subfolders (default: all_pileups)",
    )
    parser.add_argument(
        "--segment-ref-table",
        default="defaults/flu_segments_ref.csv",
        help=(
            "CSV with columns: reference,virus,segment used for segmented viruses "
            "(default: defaults/flu_segments_ref.csv)"
        ),
    )
    args = parser.parse_args()

    out_path = build_avg_json(
        args.virus_name, Path(args.input_root), Path(args.segment_ref_table)
    )
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
