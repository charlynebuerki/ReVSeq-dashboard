from __future__ import annotations

import hashlib
from pathlib import Path


def asset_is_stale(output_path: str | Path, input_paths: list[str | Path]) -> bool:
    """Return True when output is missing or older than any input file."""
    output_path = Path(output_path)
    if not output_path.exists():
        return True

    output_mtime = output_path.stat().st_mtime
    for input_path in input_paths:
        input_path = Path(input_path)
        if input_path.exists() and input_path.stat().st_mtime > output_mtime:
            return True
    return False


def asset_version(input_paths: list[str | Path]) -> str:
    """Build a short cache-busting token from input mtimes."""
    mtimes = []
    for input_path in input_paths:
        p = Path(input_path)
        mtimes.append(str(p.stat().st_mtime) if p.exists() else "missing")
    return hashlib.md5("|".join(mtimes).encode("utf-8")).hexdigest()[:10]
