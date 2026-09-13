#!/usr/bin/env python3
"""Merge Clang profiles and produce text, JSON and HTML coverage without extra packages."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("--profiles", type=Path, help="Profile directory (default: <build>/profiles)")
    parser.add_argument("--llvm-profdata", default=shutil.which("llvm-profdata") or "llvm-profdata-18")
    parser.add_argument("--llvm-cov", default=shutil.which("llvm-cov") or "llvm-cov-18")
    args = parser.parse_args()
    build = args.build.resolve()
    profiles = sorted((args.profiles or (build / "profiles")).glob("*.profraw"))
    if not profiles:
        parser.error("No profiles: run CTest with LLVM_PROFILE_FILE=<build>/profiles/%p-%m.profraw")
    binaries = []
    names = {p.stem if p.suffix == ".exe" else p.name for p in build.glob("test_*")
             if p.is_file() and p.suffix in ("", ".exe")}
    for name in sorted(names | {"chmath_odr", "chmath_example"}):
        candidate = build / name
        if not candidate.is_file():
            candidate = candidate.with_suffix(".exe")
        if candidate.is_file():
            binaries.append(candidate)
    for unified in (build / "chmath_coverage", build / "chmath_coverage.exe"):
        if unified.is_file():
            binaries = [unified]
            break
    if not binaries:
        parser.error("No test executables in build directory")
    output = build / "coverage"
    output.mkdir(exist_ok=True)
    merged = output / "merged.profdata"
    subprocess.run([args.llvm_profdata, "merge", "-sparse", *map(str, profiles), "-o", str(merged)], check=True)
    objects = [str(binaries[0])]
    for binary in binaries[1:]:
        objects += ["-object", str(binary)]
    common = [*objects, f"-instr-profile={merged}", "-ignore-filename-regex=(/tests/|/examples/|/usr/|/install-smoke/)"]
    report = subprocess.check_output([args.llvm_cov, "report", *common], text=True)
    (output / "report.txt").write_text(report, encoding="utf-8")
    exported = subprocess.check_output([args.llvm_cov, "export", *common], text=True)
    (output / "coverage.json").write_text(exported, encoding="utf-8")
    subprocess.run([args.llvm_cov, "show", *common, "-format=html", f"-output-dir={output / 'html'}", "-show-branches=count"], check=True, stdout=subprocess.DEVNULL)
    data = json.loads(exported)
    uncovered = []
    for unit in data["data"]:
        for file in unit["files"]:
            filename = file["filename"].replace("\\", "/")
            if "/src/chmath/" in filename:
                lines = sorted({s[0] for s in file["segments"] if s[2] == 0 and s[3] and s[4] and not s[5]})
                uncovered.append({"file": filename.split("/src/")[-1], "uncovered_region_start_lines": lines, "summary": file["summary"]})
    (output / "uncovered.json").write_text(json.dumps(uncovered, indent=2), encoding="utf-8")
    print(report)
    print(f"HTML: {output / 'html' / 'index.html'}")


if __name__ == "__main__":
    main()
