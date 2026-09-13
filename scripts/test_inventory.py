#!/usr/bin/env python3
"""Audit independently named test cases; do not count loop iterations as cases."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum", type=int, default=1000)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--build-dir", type=Path, help="Verify executable --list output matches source inventory")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    cases = []
    names = set()
    for source in sorted((root / "tests").glob("test_*.cpp")):
        for number, line in enumerate(source.read_text(encoding="utf-8").splitlines(), 1):
            found = re.match(r"\s*CH_TEST\(([A-Za-z0-9_]+)\)", line)
            if found:
                name = found[1]
                if name in names:
                    raise SystemExit(f"Duplicate case: {name}")
                names.add(name)
                cases.append({"name": name, "source": source.relative_to(root).as_posix(), "line": number,
                              "category": next((category for prefix, category in (("perf_", "performance"), ("stress_", "stress"), ("safety_", "safety")) if name.startswith(prefix)), source.stem)})
    report = {"count": len(cases), "cases": cases}
    if args.build_dir:
        build = args.build_dir.resolve()
        for source in sorted({case["source"] for case in cases}):
            executable = build / (Path(source).stem + (".exe" if os.name == "nt" else ""))
            actual = subprocess.check_output([str(executable), "--list"], text=True).splitlines()
            expected = [case["name"] for case in cases if case["source"] == source]
            if actual != expected:
                raise SystemExit(f"Compiled test inventory mismatch: {executable}")
        print(f"Verified compiled --list output in {build}")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"{len(cases)} unique independent cases; minimum={args.minimum}")
    if len(cases) < args.minimum:
        raise SystemExit(1)

if __name__ == "__main__":
    main()
