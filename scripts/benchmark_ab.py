#!/usr/bin/env python3
"""Run two benchmark executables in alternating order, using identical inputs."""
import argparse
import json
from pathlib import Path
import platform
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("before", type=Path)
    parser.add_argument("after", type=Path)
    parser.add_argument("--count", type=int, default=65536)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if not 16 <= args.count <= 1 << 24 or args.repeats < 1:
        parser.error("count must be 16..16777216 and repeats must be positive")
    executables = {"before": args.before.resolve(), "after": args.after.resolve()}
    for executable in executables.values():
        if not executable.is_file():
            parser.error(f"Executable not found: {executable}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    logs = {side: args.output_dir / f"{side}.txt" for side in executables}
    for log in logs.values():
        log.write_text("", encoding="utf-8")
    for repeat in range(args.repeats):
        order = ("before", "after") if repeat % 2 == 0 else ("after", "before")
        for side in order:
            completed = subprocess.run(
                [str(executables[side]), str(args.count)],
                check=True, capture_output=True, text=True, timeout=600,
            )
            with logs[side].open("a", encoding="utf-8") as log:
                log.write(completed.stdout)
            if completed.stderr:
                print(completed.stderr, file=sys.stderr, end="")
        print(f"Completed pair {repeat + 1}/{args.repeats}", flush=True)
    metadata = {
        "platform": platform.platform(), "machine": platform.machine(),
        "executables": {k: str(v) for k, v in executables.items()},
        "count": args.count, "repeats": args.repeats,
    }
    (args.output_dir / "environment.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    subprocess.run([
        sys.executable, str(Path(__file__).with_name("benchmark_compare.py")),
        str(logs["before"]), str(logs["after"]),
        "--output", str(args.output_dir / "comparison.json"),
    ], check=True)


if __name__ == "__main__":
    main()
