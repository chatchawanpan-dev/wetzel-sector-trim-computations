#!/usr/bin/env python3
"""Numerical screening of a simplified chain model for the trims H_t.

This is the experiment described in Section 1.1 of the paper.  It is not
part of the proof of any statement.

Model.  One unit sector with apex O and radii of directions 0 and 30
degrees is fixed.  A chain visits eight points:

  * L1, U2, L3: on the lower, the upper and the lower radius, in this
    parameter order;
  * U1, L2, U3: on the upper, the lower and the upper radius, in this
    parameter order;
  * Xm, Xp: a point of the closed lower cap and of the closed upper cap
    of the trim H_t.

All 8!/(3! 3!) = 1120 orders of the eight labels that keep the order within
each triple are enumerated.  For each order the length of the polygonal
chain is minimized over the ten coordinates (six radii in [0, 1] and two
coordinates for each cap point) by L-BFGS-B from deterministic starting
points.  The objective is nonsmooth and nonconvex, so the values are local
minima found by the optimizer, not certified lower bounds, and the solver
often stops with a warning at a nonsmooth point.

Usage (NumPy and SciPy are needed):

    python3 screen_cap_branches.py --candidate H16 --restarts 1 --maxiter 300
    python3 screen_cap_branches.py --all-candidates --restarts 1 --maxiter 300

A candidate is written H<t> with 0 < t < 30, for example H15 or H15.5;
--all-candidates runs H15 and H16.  The full table of results is written
as JSON to the folder "output" beside this file (or to --json-out).
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import platform
import scipy
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.optimize import minimize


REPORTS = Path(__file__).resolve().parent / "output"

LOWER_TRIPLE = ("L1", "U2", "L3")
UPPER_TRIPLE = ("U1", "L2", "U3")
CAP_LABELS = ("Xm", "Xp")
LOWER_LABELS = {"L1", "L2", "L3"}
UPPER_LABELS = {"U1", "U2", "U3"}


def p(degrees: float) -> np.ndarray:
    radians = math.radians(degrees)
    return np.array([math.cos(radians), math.sin(radians)], dtype=float)


P30 = p(30.0)


@dataclass(frozen=True)
class Candidate:
    name: str
    retained_degrees: float

    @property
    def lower_cap_end(self) -> float:
        return 30.0 - self.retained_degrees

    @property
    def upper_cap_start(self) -> float:
        return self.retained_degrees

    @property
    def upper_normal(self) -> float:
        return 0.5 * (self.retained_degrees + 30.0)

    @property
    def upper_offset(self) -> float:
        return math.cos(math.radians(0.5 * (30.0 - self.retained_degrees)))

    @property
    def lower_normal(self) -> float:
        return 0.5 * self.lower_cap_end

    @property
    def lower_offset(self) -> float:
        return math.cos(math.radians(0.5 * self.lower_cap_end))


def interleavings(a: tuple[str, ...], b: tuple[str, ...]) -> list[tuple[str, ...]]:
    """All interleavings preserving the internal order of a and b."""
    out: list[tuple[str, ...]] = []
    for positions in itertools.combinations(range(len(a) + len(b)), len(a)):
        result: list[str] = []
        ia = ib = 0
        pos_set = set(positions)
        for idx in range(len(a) + len(b)):
            if idx in pos_set:
                result.append(a[ia])
                ia += 1
            else:
                result.append(b[ib])
                ib += 1
        out.append(tuple(result))
    return out


def all_branch_orders() -> list[tuple[str, ...]]:
    """All orders preserving both Lambda-triple orders."""
    labels = [*LOWER_TRIPLE, *UPPER_TRIPLE, *CAP_LABELS]
    out: list[tuple[str, ...]] = []
    for perm in itertools.permutations(labels):
        if tuple(x for x in perm if x in LOWER_TRIPLE) != LOWER_TRIPLE:
            continue
        if tuple(x for x in perm if x in UPPER_TRIPLE) != UPPER_TRIPLE:
            continue
        out.append(perm)
    return out


def ray_point(label: str, r: float) -> np.ndarray:
    if label in LOWER_LABELS:
        return np.array([r, 0.0], dtype=float)
    if label in UPPER_LABELS:
        return r * P30
    raise KeyError(label)


def cap_point(candidate: Candidate, label: str, angle_unit: float, radial_unit: float) -> np.ndarray:
    if label == "Xm":
        angle = candidate.lower_cap_end * angle_unit
        normal = candidate.lower_normal
        offset = candidate.lower_offset
    elif label == "Xp":
        angle = candidate.upper_cap_start + (30.0 - candidate.upper_cap_start) * angle_unit
        normal = candidate.upper_normal
        offset = candidate.upper_offset
    else:
        raise KeyError(label)

    denom = math.cos(math.radians(angle - normal))
    r_min = offset / denom
    # Numerical guard against tiny roundoff above 1 at endpoints.
    r_min = min(max(r_min, 0.0), 1.0)
    radius = r_min + radial_unit * (1.0 - r_min)
    return radius * p(angle)


def unpack(candidate: Candidate, order: tuple[str, ...], x: np.ndarray) -> dict[str, np.ndarray]:
    points: dict[str, np.ndarray] = {}
    cursor = 0
    for label in (*LOWER_TRIPLE, *UPPER_TRIPLE):
        points[label] = ray_point(label, x[cursor])
        cursor += 1
    points["Xm"] = cap_point(candidate, "Xm", x[cursor], x[cursor + 1])
    cursor += 2
    points["Xp"] = cap_point(candidate, "Xp", x[cursor], x[cursor + 1])
    return points


def objective(candidate: Candidate, order: tuple[str, ...], x: np.ndarray) -> float:
    points = unpack(candidate, order, x)
    return objective_from_points(order, points)


def objective_from_points(order: tuple[str, ...], points: dict[str, np.ndarray]) -> float:
    return float(
        sum(np.linalg.norm(points[b] - points[a]) for a, b in zip(order[:-1], order[1:]))
    )


def default_start(candidate: Candidate, order: tuple[str, ...], jitter: float = 0.0) -> np.ndarray:
    # Place side contacts at increasing radii and cap witnesses near the middle
    # of their cap intervals.  The small deterministic jitter helps avoid
    # repeatedly landing on a nonsmooth coincident-contact kink.
    base = {
        "L1": 0.22,
        "U2": 0.45,
        "L3": 0.68,
        "U1": 0.26,
        "L2": 0.50,
        "U3": 0.74,
    }
    x = np.array([base[label] for label in (*LOWER_TRIPLE, *UPPER_TRIPLE)] + [0.5, 0.2, 0.5, 0.2], dtype=float)
    if jitter:
        offsets = np.sin(np.arange(1, len(x) + 1) * (0.731 + jitter)) * 0.03
        x = np.clip(x + offsets, 0.0, 1.0)
    return x


def solve_order(
    candidate: Candidate,
    order: tuple[str, ...],
    restarts: int,
    maxiter: int,
) -> dict:
    bounds = [(0.0, 1.0)] * 10
    starts = [default_start(candidate, order, 0.0)]
    starts.extend(default_start(candidate, order, float(i)) for i in range(1, restarts))

    best = None
    for start in starts:
        result = minimize(
            lambda z: objective(candidate, order, z),
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": maxiter, "ftol": 1e-12, "gtol": 1e-9},
        )
        value = float(result.fun)
        if best is None or value < best["value"]:
            best = {
                "value": value,
                "success": bool(result.success),
                "message": str(result.message),
                "x": result.x.tolist(),
                "nit": int(result.nit),
            }
    assert best is not None
    variables = np.array(best["x"], dtype=float)
    points = unpack(candidate, order, variables)
    value_from_variables = objective(candidate, order, variables)
    value_from_points = objective_from_points(order, points)
    solver_value = best["value"]
    value = value_from_variables
    return {
        "order": list(order),
        "value": value,
        "solver_value": solver_value,
        "value_from_variables": value_from_variables,
        "value_from_points": value_from_points,
        "solver_recompute_delta": abs(solver_value - value_from_variables),
        "points_recompute_delta": abs(value_from_points - value_from_variables),
        "success": best["success"],
        "message": best["message"],
        "nit": best["nit"],
        "variables": best["x"],
        "points": {label: point.tolist() for label, point in points.items()},
    }


def summarize(results: list[dict], threshold: float, margin: float) -> dict:
    values = [r["value"] for r in results]
    above_threshold = [r for r in results if r["value"] > threshold + margin]
    open_rows = [r for r in results if r["value"] <= threshold + margin]
    solver_warnings = [r for r in results if not r["success"]]
    weakest = min(results, key=lambda r: r["value"])
    return {
        "num_branches": len(results),
        "num_output_above_threshold": len(above_threshold),
        "num_open": len(open_rows),
        "num_at_or_below_threshold": len(open_rows),
        "num_solver_warnings": len(solver_warnings),
        "threshold": threshold,
        "margin": margin,
        "min_value": min(values),
        "median_value": float(np.median(values)),
        "max_value": max(values),
        "weakest_order": weakest["order"],
        "weakest_value": weakest["value"],
        "screening_clears_threshold": len(open_rows) == 0,
    }


def provenance(args: argparse.Namespace) -> dict:
    return {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "argv": sys.argv,
        "parameters": {
            "all_candidates": bool(args.all_candidates),
            "candidate": args.candidate.name if hasattr(args.candidate, "name") else str(args.candidate),
            "restarts": args.restarts,
            "maxiter": args.maxiter,
            "threshold": args.threshold,
            "margin": args.margin,
        },
        "determinism": "Deterministic starts; no random number generator is used.",
    }


def candidate_from_name(name: str) -> Candidate:
    normalized = name.upper()
    if normalized == "H15":
        return Candidate("H15", 15.0)
    if normalized == "H16":
        return Candidate("H16", 16.0)
    try:
        t = float(normalized.lstrip("H"))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("candidate must be H<t> with 0 < t < 30, for example H15 or H12.5") from exc
    if not 0.0 < t < 30.0:
        raise argparse.ArgumentTypeError("candidate must be H<t> with 0 < t < 30")
    return Candidate(f"H{t:g}", t)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=candidate_from_name, default=Candidate("H16", 16.0))
    parser.add_argument("--all-candidates", action="store_true")
    parser.add_argument("--restarts", type=int, default=4)
    parser.add_argument("--maxiter", type=int, default=500)
    parser.add_argument("--threshold", type=float, default=1.0)
    parser.add_argument("--margin", type=float, default=1e-5)
    parser.add_argument("--json-out", type=Path, default=None)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    REPORTS.mkdir(parents=True, exist_ok=True)
    candidates = [Candidate("H15", 15.0), Candidate("H16", 16.0)] if args.all_candidates else [args.candidate]
    orders = all_branch_orders()

    payload = {
        "model": "two_lambda_triples_with_two_cap_witnesses",
        "model_notes": [
            "Two Lambda triples are represented by labels L1,U2,L3 and U1,L2,U3.",
            "The all_branches ledger contains every strict order preserving each triple order.",
            "Xm is the lower mirror-cap witness and Xp is the upper cap witness.",
            "This reduced nonlinear optimization layer is a screening experiment, not a global lower-bound certificate or proof.",
            "Each row value is recomputed from the serialized variables and points before writing.",
        ],
        "num_orders_per_candidate": len(orders),
        "provenance": provenance(args),
        "candidates": {},
    }

    t0 = time.time()
    for candidate in candidates:
        results = []
        for idx, order in enumerate(orders, start=1):
            row = solve_order(candidate, order, args.restarts, args.maxiter)
            row["id"] = f"{candidate.name}-{idx:04d}"
            results.append(row)
            if not args.quiet and (idx == 1 or idx % 100 == 0 or idx == len(orders)):
                print(f"{candidate.name}: solved {idx}/{len(orders)}; current min={min(r['value'] for r in results):.12f}")
        summary = summarize(results, args.threshold, args.margin)
        payload["candidates"][candidate.name] = {
            "retained_degrees": candidate.retained_degrees,
            "summary": summary,
            "all_branches": results,
            "open_branches": [
                r for r in results if r["value"] <= args.threshold + args.margin
            ],
            "solver_warning_branches": [
                r for r in results if not r["success"]
            ],
            "weakest_branches": sorted(results, key=lambda r: r["value"])[:20],
        }
        if not args.quiet:
            print(f"{candidate.name} summary: {json.dumps(summary, indent=2)}")
        weakest = min(results, key=lambda r: r["value"])
        chain = ", ".join(
            f"{label}=({weakest['points'][label][0]:.6f}, {weakest['points'][label][1]:.6f})"
            for label in weakest["order"]
        )
        print(f"{candidate.name}: smallest length found {weakest['value']:.9f}")
        print(f"{candidate.name}: attained by the chain {chain}")

    payload["elapsed_seconds"] = time.time() - t0
    output = args.json_out
    if output is None:
        suffix = "all" if args.all_candidates else args.candidate.name.lower()
        output = REPORTS / f"cap_branch_optimization_{suffix}.json"
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
