#!/usr/bin/env python3
"""Portable entry point for the current H15 proof checks (standard library)."""
from __future__ import annotations

import datetime
import platform
import subprocess
import sys
import time
from pathlib import Path

if sys.flags.optimize:
    raise SystemExit("ERROR: Python optimization disables required proof checks; run without -O/-OO or PYTHONOPTIMIZE.")

BASE = Path(__file__).resolve().parent


def main():
    start = time.monotonic()
    print("H15 proof checks: global Case 2 formula and exact mixed-state checks", flush=True)
    print("UTC:", datetime.datetime.now(datetime.timezone.utc).isoformat(), flush=True)
    print("Python:", sys.version.replace("\n", " "), flush=True)
    print("Platform:", platform.platform(), flush=True)
    steps = [
        ("Localized Cases 3/4 geometric cap-cap", "Geometric/verify_capcap_circles.py", []),
        ("Case 2 cap-cap: one global coefficient rule and fixed algebraic bounds", "Geometric/verify_case2_sine.py", []),
        ("Case 2 actual mixed: 12 words and all radial directions", "Case2/verify_case2_actual_mixed_prefix.py", ["case2_actual_mixed_prefix_data.json"]),
        ("Case 3 actual mixed: 12 words and all radial directions", "Case3/verify_case3_mixed.py", ["case3_mixed_certificate_data.json"]),
        ("Case 4 actual mixed: 12 words and all radial directions", "Case4/verify_case4_mixed.py", ["case4_mixed_certificate_data.json"]),
    ]
    for title, relative_path, arguments in steps:
        script = BASE / relative_path
        print(f"\n[{title}]", flush=True)
        subprocess.run([sys.executable, str(script), *arguments], cwd=script.parent, check=True)
    print("\nH15 PROOF CHECKS: PASS", flush=True)
    print("Elapsed seconds:", round(time.monotonic() - start, 3), flush=True)
    print("The canonical-wedge, actual-witness, symmetry, radial-radial, and analytic equality reductions remain mathematical dependencies in the manuscript.", flush=True)


if __name__ == "__main__":
    main()
