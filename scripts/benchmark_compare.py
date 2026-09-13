#!/usr/bin/env python3
"""Compare medians across repeated chmath_bench process runs (ns/item)."""
import argparse
import json
from pathlib import Path
import re
from statistics import median


def read(path):
    values = {}
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        match = re.match(r"(.+?)\s+([0-9]+\.[0-9]+) ns/item\s", line)
        if match:
            values.setdefault(match[1].strip(), []).append(float(match[2]))
    if not values:
        raise SystemExit(f"No benchmark measurements: {path}")
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("before", type=Path)
    parser.add_argument("after", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    before, after = read(args.before), read(args.after)
    if before.keys() != after.keys():
        raise SystemExit("Before and after measurement names differ")
    result = {}
    for name, old in before.items():
        new = after[name]
        a, b = median(old), median(new)
        if min(a, b) <= 0:
            raise SystemExit(f"Nonpositive measurement: {name}")
        result[name] = {"before_ns": a, "after_ns": b, "speedup": a / b,
                        "time_reduction_percent": 100 * (a - b) / a,
                        "before_samples": old, "after_samples": new}
        print(f"{name:28s} {a:9.3f} -> {b:9.3f} ns/item  {a/b:6.3f}x")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
