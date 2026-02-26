#!/usr/bin/env python3
"""
Reorder a FASTA to match the required convention:

1) SUPER_* groups first, sorted by TOTAL group length (main + all _unloc_*).
2) Within each group: main first, then its _unloc_* in lexicographic order.
3) All non-SUPER_* records afterwards, sorted by individual record length (desc).

- Preserves headers and sequence text exactly (copies record byte ranges; no re-wrapping).
- Memory-safe: one pass to index offsets and lengths; second pass to copy blocks.
"""

import os
import sys
import re
from typing import List, Tuple, Dict

UNLOC_RE = re.compile(r'^(?P<prefix>[^\s>]+)_unloc_\d+$')

Record = Tuple[bytes, int, int, int]  # (header_line, seq_len, start_offset, end_offset)


def die(msg: str, code: int = 1) -> None:
    sys.stderr.write(msg + "\n")
    sys.exit(code)


def index_fasta(path: str) -> List[Record]:
    """Return list of (header_line, seq_len, start, end) for each record."""
    records: List[Record] = []
    with open(path, "rb") as fh:
        pos = 0
        current_header: bytes = b""
        current_start: int = -1
        seqlen: int = 0

        def finalize(end_pos: int) -> None:
            nonlocal current_header, current_start, seqlen
            if current_header:
                records.append((current_header, seqlen, current_start, end_pos))
                current_header = b""
                current_start = -1
                seqlen = 0

        while True:
            line = fh.readline()
            if not line:
                finalize(pos)
                break
            line_start = pos
            pos += len(line)
            if line.startswith(b">"):
                finalize(line_start)
                current_header = line.rstrip(b"\r\n")
                current_start = line_start
                seqlen = 0
            else:
                # Count non-whitespace characters in sequence lines
                seqlen += len(line.strip())

    return records


def header_id(header_line: bytes) -> str:
    """FASTA ID = first token after '>' up to first whitespace."""
    text = header_line[1:].decode("utf-8", errors="ignore")
    return text.split()[0]


def compute_final_order(records: List[Record]) -> List[str]:
    """Compute ID order according to the required grouping rules."""
    lengths: Dict[str, int] = {}
    idx_bounds: Dict[str, Tuple[int, int]] = {}
    ids: List[str] = []

    for hdr, seqlen, start, end in records:
        rid = header_id(hdr)
        ids.append(rid)
        idx_bounds[rid] = (start, end)
        lengths[rid] = seqlen

    # Build SUPER_* groups (prefix -> {main, unlocs})
    groups: Dict[str, Dict[str, List[str] or str]] = {}
    for rid in ids:
        m = UNLOC_RE.match(rid)
        if m:
            prefix = m.group("prefix")
            groups.setdefault(prefix, {"main": None, "unlocs": []})
            groups[prefix]["unlocs"].append(rid)
        elif rid.startswith("SUPER_"):
            groups.setdefault(rid, {"main": None, "unlocs": []})
            groups[rid]["main"] = rid

    # Total length per group (main + all unlocs present)
    group_totals: List[Tuple[str, int]] = []
    for prefix, g in groups.items():
        total = 0
        main = g["main"]
        if main:
            total += lengths.get(main, 0)
        for u in g["unlocs"]:
            total += lengths.get(u, 0)
        group_totals.append((prefix, total))

    # Groups sorted by total length (desc), then by prefix
    groups_sorted = sorted(group_totals, key=lambda x: (-x[1], x[0]))

    final_ids: List[str] = []
    seen: set = set()

    print("sorted groups:", groups_sorted)

    for prefix, _ in groups_sorted:
        g = groups[prefix]
        main = g["main"]
        if main and main in lengths:
            final_ids.append(main)
            seen.add(main)
        for u in sorted(g["unlocs"]):
            if u in lengths:
                final_ids.append(u)
                seen.add(u)

    print("Final IDs:", final_ids, len(final_ids))
    print("Seen:", seen, len(seen))

    # Append non-SUPER_* records by individual length (desc)
    others = [(rid, lengths[rid]) for rid in ids if rid not in seen]
    others_sorted = sorted(others, key=lambda x: (-x[1], x[0]))
    final_ids.extend([rid for rid, _ in others_sorted])

    return final_ids


def write_by_offsets(in_path: str, out_path: str, records: List[Record], final_order: List[str]) -> None:
    # Map id -> (start, end)
    bounds = {header_id(h): (s, e) for (h, _, s, e) in records}
    tmp_path = f"{out_path}.tmp"
    with open(in_path, "rb") as src, open(tmp_path, "wb") as out:
        for rid in final_order:
            s, e = bounds[rid]
            src.seek(s)
            block = src.read(e - s)
            # ensure record ends with newline between records
            out.write(block)
            if not block.endswith(b"\n"):
                out.write(b"\n")
    # Sanity
    with open(tmp_path, "rb") as fh:
        if fh.read(1) != b">":
            os.remove(tmp_path)
            die("ERROR: Output contains no FASTA header; refusing to write empty output.")
    os.replace(tmp_path, out_path)


def main(argv: List[str]) -> None:
    if len(argv) != 3:
        die("Usage: fasta_reorder.py <input.fa> <output.fa>")
    in_path, out_path = argv[1], argv[2]
    if not os.path.isfile(in_path):
        die(f"ERROR: Input not found: {in_path}")
    if os.path.exists(out_path) and os.path.samefile(in_path, out_path):
        die("ERROR: Input and output paths are the same.")

    sys.stderr.write("Indexing lengths…\n")
    recs = index_fasta(in_path)
    sys.stderr.write(f"Found {len(recs)} records.\n")
    if not recs:
        die("ERROR: No FASTA records read from input; refusing to write an empty output.")

    final_order = compute_final_order(recs)
    sys.stderr.write(f"Writing: {out_path}\n")
    write_by_offsets(in_path, out_path, recs, final_order)


if __name__ == "__main__":
    main(sys.argv)
