#!/usr/bin/env python3
"""Check the distributed snapshot against its current top-level SHA-256 inventory.

Additional locally generated build files and run logs are allowed. This checks
integrity of the recorded snapshot, not a digital signature or proof correctness.
"""
import hashlib
import json
from pathlib import Path


def main():
    root = Path(__file__).resolve().parent
    payload = json.loads((root / "SHA256SUMS.json").read_text(encoding="utf-8"))
    if payload.get("algorithm") != "sha256" or not isinstance(payload.get("files"), dict):
        raise SystemExit("FAIL: invalid snapshot inventory")
    files = payload["files"]
    if not files:
        raise SystemExit("FAIL: empty snapshot inventory")
    failures = []
    for name, expected in files.items():
        path = root / name
        if path.is_absolute() and not path.resolve().is_relative_to(root):
            failures.append(f"invalid path: {name}")
        elif not path.is_file():
            failures.append(f"missing: {name}")
        elif hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            failures.append(f"mismatch: {name}")
    if failures:
        raise SystemExit("FAIL: snapshot integrity\n" + "\n".join(failures))
    print(f"PASS: {len(files)} snapshot files match their SHA-256 hashes")


if __name__ == "__main__":
    main()
