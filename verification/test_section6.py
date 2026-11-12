#!/usr/bin/env python3
"""Control tests for Section6/verify_zarc_test.py.

The checker must accept the inequalities stated in the paper and must
reject altered ones for the intended reason:

  1. with the values for t = 3.7 (chi = 16.85, delta = 13.15) the test of
     Proposition 6.6 must fail at theta = beta for reflected copies, since
     the arc Z does fit in H_t for t = 3.7;
  2. with the single chord chi = 15 on the whole range [0, 180], the test
     for the pentagon must fail (two chords are needed in Remark 6.8);
  3. the triangle conv{O, p_0, p_30} (chi = 15, h = cos 15) must pass.

Only the Python standard library is used.  Run without -O.
"""
from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
from fractions import Fraction as F
from pathlib import Path

if sys.flags.optimize:
    raise SystemExit("ERROR: run without -O/-OO or PYTHONOPTIMIZE; the checks are assertions.")

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("verify_zarc_test", HERE / "Section6" / "verify_zarc_test.py")
checker = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(checker)


def quiet(function, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return function(*args, **kwargs)


def expect_failure(label, function, *args, expected_point, **kwargs):
    try:
        quiet(function, *args, **kwargs)
    except AssertionError as error:
        details = error.args[0] if error.args else ()
        point = details[1] if isinstance(details, tuple) and len(details) > 1 else None
        if point != expected_point:
            raise SystemExit(f"FAIL: {label}: rejected at {point!r}, expected {expected_point!r}")
        print(f"PASS rejection: {label}: fails at theta = {point}")
        return
    raise SystemExit(f"FAIL: {label}: altered input was accepted")


def main() -> int:
    # positive control: the distributed checker passes
    assert quiet(checker.main) == 0
    print("PASS: distributed Section 6 checks")

    # positive control: the triangle conv{O, p_0, p_30} (t = 0)
    for eps in (+1, -1):
        count, _ = quiet(checker.check_family, "H_0", F(15), eps, F(15))
        assert count == 6
    print("PASS: the arc Z does not fit in the triangle conv{O, p_0, p_30}")

    # rejection 1: t = 3.7
    quiet(checker.check_family, "H_3.7", F(337, 20), +1, F(263, 20))  # direct copies still excluded
    expect_failure("t = 3.7, reflected copies", checker.check_family,
                   "H_3.7", F(337, 20), -1, F(263, 20), expected_point="beta")

    # rejection 2: one chord is not enough for the pentagon
    expect_failure("pentagon with the chord chi = 15 only", checker.check_family,
                   "K10", 15, +1, 5, expected_point="210 - beta")

    print("PASS: 2 positive controls and 2 rejection controls")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
