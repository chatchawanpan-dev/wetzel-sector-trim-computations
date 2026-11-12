# Reproducibility guide

All commands below are run from the supplement root, the directory containing `verification/`, `validation/`, and `verify_snapshot.py`.

## Exact verification

Use Python 3 with its standard library. The recorded Python version appears in the execution logs in `validation/`. No numerical optimizer, coefficient search, plotting package, or network access is required to check the proof. Keep `PYTHONOPTIMIZE` unset and do not use `-O` or `-OO`.

From the supplement root:

```sh
python3 verification/verify_proof_essential.py
```

The driver resolves files relative to its own location and starts each checker in that checker’s directory. An absolute invocation of the driver therefore works from any working directory. It stops at the first unsuccessful stage and refuses optimized Python. A complete successful run exits with status zero and prints `H15 PROOF CHECKS: PASS` after all five stages.

The execution order is Cases 3/4 cap–cap circles, Case 2 cap–cap sine rule, Case 2 mixed, Case 3 mixed, and Case 4 mixed. [PROOF_ESSENTIAL.md](PROOF_ESSENTIAL.md) describes every certificate and its mathematical dependencies. [PROOF_ESSENTIAL_FILES.json](PROOF_ESSENTIAL_FILES.json) supplies the nine required files and their hashes.

## Individual stages

These commands are for the supplement root. The mixed-data argument is explicit because the checkers’ default filenames are relative to the working directory.

```sh
python3 verification/Geometric/verify_capcap_circles.py
python3 verification/Geometric/verify_case2_sine.py
python3 verification/Case2/verify_case2_actual_mixed_prefix.py verification/Case2/case2_actual_mixed_prefix_data.json
python3 verification/Case3/verify_case3_mixed.py verification/Case3/case3_mixed_certificate_data.json
python3 verification/Case4/verify_case4_mixed.py verification/Case4/case4_mixed_certificate_data.json
```

The two geometric programs write their check summaries beside their source files. Such output summaries are not runtime inputs.

## Domain, arithmetic, and acceptance criteria

Each mixed case verifies the full closed rectangle `theta in [150,210]` degrees and `phi in [0,30]` degrees. Coefficient rectangles must cover the domain, including its boundaries. Coefficients are parsed as exact rational values from strings; symbolic `-sqrt(3)*lambda` entries preserve that formula. Interval arithmetic uses outward rational enclosures. Trigonometric bounds use quadrant reduction, Machin’s formula for pi, and alternating Taylor bounds with explicit remainders.

Each mixed case has twelve orders of nine labels, whose proper prefixes form thirteen distinct label sets; each order has eight proper prefixes. Acceptance requires the complete rectangle inventory, valid endpoints and coefficient conditions, positive multipliers, full coverage, and every required squared prefix norm at most one, either by interval certification or by an analytic case with verified hypotheses.

| Case | Required rectangles | Interval comparisons | Analytic comparisons |
| --- | ---: | ---: | ---: |
| 2 | 33 | 427 | 2 |
| 3 | 20 | 237 | 23 |
| 4 | 21 | 250 | 23 |
| Total | 74 | 914 | 48 |

The Case 2 sine checker verifies formal rational-polynomial identities and fixed algebraic bounds on the entire orientation interval, with no angular sampling or subdivision. All nine prefix norms are strictly below `49/50`, giving the closed-model length bound `50/49`. The Cases 3/4 circle checker performs sixteen circular-arc comparisons with norm strictly below `999/1000`, bounds the constant-prefix forms by `3/4`, and retains the analytic critical chord bound one. In the mixed argument the strict inequality for the length comes from positive cap and radial slack, also when a prefix norm equals one.

## Section 6.2

From the supplement root:

```sh
python3 verification/Section6/verify_zarc_test.py
```

This checker is independent of the five-stage proof route and of the data files. For an angle `chi` it encloses

```
F(chi, eps, theta) = f_chi(theta)/12 + sqrt(111) f_chi(theta + eps*beta)/12,
f_chi(phi) = sin(120 - chi)|sin(phi)| + cos(chi)|cos(phi - 120)| + (1/2)|cos(phi - chi)|,
```

where `cos(beta) = 6/sqrt(111)` and `sin(beta) = 5 sqrt(3)/sqrt(111)`, at the points `0, 30, 90 + chi` and `-eps*beta, 30 - eps*beta, 90 + chi - eps*beta` modulo 180 degrees. Sines and cosines of rational angles are enclosed as in the mixed checkers; `sqrt(3)` and `sqrt(111)` are enclosed from integer square roots, and angles containing `beta` are handled by the addition formulas. All end points are rational and rounded outward to 40 decimal places.

Acceptance requires:

1. `F(16.8, eps, theta) > cos(13.2)` at all twelve points, for `eps = 1` and `eps = -1` (Proposition 6.6 and Table 4);
2. `<Q_3, p_16.85> < cos(13.15)` and `|Q_3|^2 = 5317/5328 < 1` (Proposition 6.6);
3. `F(15, 1, theta) > cos(5)` at the five points in `[0, 140]` and at `140`, and `F(25, 1, theta) > cos(5)` at the one point in `[140, 180]` and at `140` and `180` (Remark 6.8). The program also checks that every other point lies strictly outside the interval in question.

A successful run prints `SECTION 6 CHECKS: PASS`.

## Control tests and execution records

From the supplement root:

```sh
python3 verification/test_geometric_certificates.py
python3 verification/test_mixed_certificates.py
python3 verification/test_section6.py
```

The geometric suite must report three positive and twenty-three rejection controls. It checks the Case 2 sine program under normal and optimized execution and confirms the circle program’s optimization refusal and its intended rejection conditions. This control coverage does not authorize optimized execution of the complete proof route.

The mixed suite must report all twelve intended rejections: deleting a rectangle, creating an angular coverage gap, making a weight negative, and changing the sign of the coefficient `c` of one rectangle, in each of Cases 2, 3, and 4. In the last control a prefix norm exceeds one, and the checker must fail at its bisection limit. A control succeeds only if the altered input fails for the expected reason. Both suites use temporary copies and preserve distributed source and data files.

The Section 6 suite must report two positive and two rejection controls. It confirms that the checker accepts the triangle `conv{O, p_0, p_30}`, that it rejects the values for `t = 3.7` at `theta = beta` for reflected copies, where the arc does fit, and that the single chord `chi = 15` does not suffice for the pentagon.

The distributed execution logs and machine-readable summary are in `validation/`. Consult them for the interpreter version, individual stage output, exit codes, and control results. When reproducing the checks, require the success conditions above rather than an exact match of timings or platform descriptions.

Run times depend on the machine and are not acceptance thresholds. Executable checks do not replace the manuscript’s geometric reductions or independent human review.

## Independent re-check

From the supplement root:

```sh
python3 verification/Independent/recheck_with_mpmath.py
```

The program needs NumPy and mpmath. Part A is rigorous: it reads the three coefficient files, checks the coverage of `[150,210] x [0,30]`, and bounds every squared prefix norm with mpmath intervals at 120 bits, bisecting the longer side when a bound exceeds one. It must report 427, 237 and 250 pairs settled by interval bounds and 2, 23 and 23 by Lemma 5.5, as in Table 3 of the paper. In the recorded run it also visits the same numbers of rectangles as the checkers (3529, 1681 and 2568) and reaches the same depths (13, 9 and 10). Part B is floating-point sampling and proves nothing: prefix norms on a 41 x 41 grid in each rectangle, the value of `Phi` on random chains, and the numbers of Appendix A.1, Figure 6 and Table 5. A successful run prints `INDEPENDENT RE-CHECK: PASS`.

## Experiments

The programs in `exploration/` are the numerical experiments of Section 1.1. They are not needed to check any statement of the paper; `exploration/README.md` describes them and records their output.

## Snapshot integrity

From the supplement root, run:

```sh
python3 verify_snapshot.py
```

This compares distributed files against `SHA256SUMS.json` in the supplement root. Extra locally generated outputs are allowed. The manifest itself is the comparison record; hash agreement is not a signature or a proof of mathematical correctness. The separate nine-file runtime inventory identifies exactly which source and data files the proof driver requires.
