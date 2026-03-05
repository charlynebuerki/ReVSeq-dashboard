#!/usr/bin/env python3
"""Build per-sample pileup traces JSON for a virus folder.

Usage:
    python Extract-pileup/build_pileup_individual_json.py cov-229e

Reads all depth files in:
    Extract-pileup/all_pileups/<virus_name>/
Writes:
    Extract-pileup/all_pileups/<virus_name>/<virus_name>_indiv.json

Output format:
{
  "virus": "cov-229e",
  "positions": [1,2,3,...],
  "samples": {
    "sample_a": [12, 11, null, ...],
    "sample_b": [8, 7, 9, ...]
  },
  "n_samples": 42
}
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


def _collect_input_files(folder: Path, virus_name: str) -> List[Path]:
    skip_suffixes = {f"{virus_name}_avg.json", f"{virus_name}_indiv.json", f"{virus_name}_individual.json"}
    files = [
        p
        for p in sorted(folder.iterdir())
        if p.is_file()
        and not p.name.startswith(".")
        and p.suffix.lower() in {".txt", ".tsv", ".csv"}
        and p.name not in skip_suffixes
    ]
    return files


def _sample_id_from_path(path: Path) -> str:
    stem = path.stem
    # Typical files: <sample>_depth.tsv
    if stem.endswith("_depth"):
        return stem[: -len("_depth")]
    # Keep full stem for pileup_<...>.txt
    return stem


def _iter_depth_file(path: Path):
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


def build_individual_json(
    virus_name: str,
    input_root: Path,
    segment_ref_table: Path,
    fill_missing_with_zero: bool = False,
) -> Path:
    folder = input_root / virus_name
    if not folder.exists() or not folder.is_dir():
        raise FileNotFoundError(f"Virus folder not found: {folder}")

    input_files = _collect_input_files(folder, virus_name)
    if not input_files:
        raise FileNotFoundError(f"No pileup input files found in: {folder}")

    ref_to_segment, segment_names = _load_segment_ref_map(segment_ref_table, virus_name)
    segment_map_active = bool(ref_to_segment or segment_names)

    per_sample: Dict[str, Dict[int, float]] = {}
    all_positions = set()
    per_segment_sample: Dict[str, Dict[str, Dict[int, float]]] = {}
    per_segment_positions: Dict[str, set[int]] = {}
    unknown_ref_counts: Dict[str, int] = {}

    for path in input_files:
        sample_id = _sample_id_from_path(path)
        file_depth_map: Dict[int, float] = {}
        file_segment_depths: Dict[str, Dict[int, float]] = {}

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
                file_segment_depths.setdefault(segment, {})[pos] = depth
            else:
                file_depth_map[pos] = depth

        if segment_map_active:
            for segment, depth_map in file_segment_depths.items():
                if not depth_map:
                    continue
                per_segment_sample.setdefault(segment, {})[sample_id] = depth_map
                per_segment_positions.setdefault(segment, set()).update(depth_map.keys())
        else:
            if not file_depth_map:
                continue
            per_sample[sample_id] = file_depth_map
            all_positions.update(file_depth_map.keys())

    if unknown_ref_counts:
        for ref_name, count in sorted(unknown_ref_counts.items()):
            print(f"[build_pileup_individual_json] Skipping unknown segment/reference '{ref_name}' ({count} rows)")

    if not per_sample and not per_segment_sample:
        raise ValueError(f"No valid depth values parsed from files in: {folder}")

    if per_segment_sample:
        segments_out: Dict[str, dict] = {}
        for segment in sorted(per_segment_sample.keys()):
            positions = sorted(per_segment_positions.get(segment, set()))
            sample_traces: Dict[str, List[float | None]] = {}
            for sample_id in sorted(per_segment_sample[segment].keys()):
                depth_map = per_segment_sample[segment][sample_id]
                if fill_missing_with_zero:
                    trace = [float(depth_map.get(pos, 0.0)) for pos in positions]
                else:
                    trace = [float(depth_map[pos]) if pos in depth_map else None for pos in positions]
                sample_traces[sample_id] = trace

            segments_out[segment] = {
                "positions": positions,
                "samples": sample_traces,
                "n_samples": len(sample_traces),
            }

        output = {
            "virus": virus_name,
            "is_segmented": True,
            "segments": segments_out,
        }
    else:
        positions = sorted(all_positions)
        samples_out: Dict[str, List[float | None]] = {}
        for sample_id in sorted(per_sample.keys()):
            depth_map = per_sample[sample_id]
            if fill_missing_with_zero:
                trace = [float(depth_map.get(pos, 0.0)) for pos in positions]
            else:
                trace = [float(depth_map[pos]) if pos in depth_map else None for pos in positions]
            samples_out[sample_id] = trace

        output = {
            "virus": virus_name,
            "positions": positions,
            "samples": samples_out,
            "n_samples": len(samples_out),
        }

    out_path = folder / f"{virus_name}_indiv.json"
    with out_path.open("w", encoding="utf-8") as fh:
        json.dump(output, fh, ensure_ascii=True)

    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Merge pileup depth files into one per-sample traces JSON for interactive plotting "
            "(positions + samples)."
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
    parser.add_argument(
        "--fill-missing-with-zero",
        action="store_true",
        help="Fill missing sample/position values with 0 instead of null.",
    )
    args = parser.parse_args()

    out_path = build_individual_json(
        virus_name=args.virus_name,
        input_root=Path(args.input_root),
        segment_ref_table=Path(args.segment_ref_table),
        fill_missing_with_zero=args.fill_missing_with_zero,
    )
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
