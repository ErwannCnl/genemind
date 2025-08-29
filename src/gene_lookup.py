# gene_lookup.py
from __future__ import annotations
import csv
from pathlib import Path
from typing import Optional, Dict

def _to_float(s: str) -> Optional[float]:
    """Best-effort numeric parse; returns None if not numeric."""
    try:
        return float(s.replace(",", ""))
    except Exception:
        return None

def get_gene_proportion(
    file_path: str | Path,
    gene_name: str,
    *,
    gene_col: int = 6,        # 0-based index; 7th column by default
    value_col: int = -1,      # last column by default
    delimiter: str = "\t",
    encoding: str = "utf-8",
    skip_header: bool = False,
    case_insensitive: bool = False,
) -> Optional[float]:
    """
    Return the queried gene's value divided by the maximum value in the file.
    - One streaming pass computes max and captures the gene's value.
    - Returns None if the gene isn't found or no numeric values exist.
    """
    p = Path(file_path)
    if not p.exists():
        raise FileNotFoundError(f"File not found: {p}")

    target = gene_name.casefold() if case_insensitive else gene_name
    max_val: Optional[float] = None
    gene_val: Optional[float] = None

    with p.open(newline="", encoding=encoding) as f:
        reader = csv.reader(f, delimiter=delimiter)
        if skip_header:
            next(reader, None)

        for row in reader:
            if not row:
                continue
            try:
                gsym = row[gene_col]
                val_raw = row[value_col]
            except IndexError:
                continue

            val = _to_float(val_raw)
            if val is not None:
                max_val = val if max_val is None or val > max_val else max_val

            if case_insensitive:
                if gsym.casefold() == target:
                    gene_val = val
            else:
                if gsym == target:
                    gene_val = val

    if gene_val is None or max_val is None:
        return None
    if max_val == 0:
        # All values are zero → define proportion of zero-valued genes as 0.0
        return 0.0
    return gene_val / max_val


def build_gene_proportion_lookup(
    file_path: str | Path,
    *,
    gene_col: int = 6,
    value_col: int = -1,
    delimiter: str = "\t",
    encoding: str = "utf-8",
    skip_header: bool = False,
    case_insensitive: bool = False,
    prefer: str = "first",    # "first" or "last" if duplicate symbols appear
) -> Dict[str, float]:
    """
    Build {gene_symbol: proportion} for all genes.
    Two passes: (1) find max, (2) compute each proportion.
    """
    p = Path(file_path)
    if not p.exists():
        raise FileNotFoundError(f"File not found: {p}")

    # Pass 1: find max
    max_val: Optional[float] = None
    with p.open(newline="", encoding=encoding) as f:
        reader = csv.reader(f, delimiter=delimiter)
        if skip_header:
            next(reader, None)
        for row in reader:
            if not row:
                continue
            try:
                val = _to_float(row[value_col])
            except IndexError:
                continue
            if val is not None:
                max_val = val if max_val is None or val > max_val else max_val

    if max_val is None or max_val == 0:
        # No numeric values or all zeros → everything maps to 0.0
        max_val = 0.0

    # Pass 2: build proportions
    props: Dict[str, float] = {}
    with p.open(newline="", encoding=encoding) as f:
        reader = csv.reader(f, delimiter=delimiter)
        if skip_header:
            next(reader, None)
        for row in reader:
            if not row:
                continue
            try:
                gsym = row[gene_col]
                val_raw = row[value_col]
            except IndexError:
                continue

            val = _to_float(val_raw)
            if val is None:
                continue

            key = gsym.casefold() if case_insensitive else gsym
            prop = 0.0 if max_val == 0.0 else (val / max_val)

            if prefer == "first":
                props.setdefault(key, prop)
            else:  # "last"
                props[key] = prop

    return props
