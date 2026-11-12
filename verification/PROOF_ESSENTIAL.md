# Proof-essential verification guide

The proof of Theorem 1.1 combines the manuscript’s geometric reductions with one five-stage exact verification route. Run `python3 verification/verify_proof_essential.py` from the supplement root, without Python optimization. Only the Python standard library is required.

All commands and inventory paths below are relative to the supplement root unless stated otherwise.

## Five stages in execution order

| Stage | Checker | Certificate and acceptance |
| --- | --- | --- |
| 1 | `Geometric/verify_capcap_circles.py` | Cases 3 and 4 cap–cap: sixteen noncritical circular-arc comparisons give norm strictly below `999/1000`; four constant-prefix forms have norm at most `3/4`. The critical chord prefix has the analytic bound one. |
| 2 | `Geometric/verify_case2_sine.py` | Case 2 cap–cap: the single rule `a = -(2 + sin(theta))/3`, `c = -(1 + sin(theta))/3` works throughout `150 <= theta <= 210` degrees. The four orders of Lemma 4.5(c) give nine distinct sets of proper-prefix labels. Exact identities and fixed algebraic estimates give all nine norms strictly below `49/50`, yielding the closed-model chain-length bound `50/49 > 1`. |
| 3 | `Case2/verify_case2_actual_mixed_prefix.py` | Case 2 mixed: all 33 coefficient rectangles, 427 interval comparisons, and 2 analytic comparisons. |
| 4 | `Case3/verify_case3_mixed.py` | Case 3 mixed: all 20 coefficient rectangles, 237 interval comparisons, and 23 analytic comparisons. |
| 5 | `Case4/verify_case4_mixed.py` | Case 4 mixed: all 21 coefficient rectangles, 250 interval comparisons, and 23 analytic comparisons. |

Paths in the table are relative to `verification/`. Each mixed case covers the entire closed domain `[150,210]` degrees in `theta` times `[0,30]` degrees in `phi`. Each case has twelve orders of nine labels, whose proper prefixes form thirteen distinct label sets; each order has eight proper prefixes. The 74 rectangles account for 914 interval comparisons and 48 analytic comparisons, or 962 rectangle-prefix comparisons in total. Every rectangle is part of the proof.

## Complete runtime inventory

[PROOF_ESSENTIAL_FILES.json](PROOF_ESSENTIAL_FILES.json) gives the authoritative nine-file runtime inventory, byte sizes, and SHA-256 hashes. The files are:

1. `verification/verify_proof_essential.py`
2. `verification/Geometric/verify_capcap_circles.py`
3. `verification/Geometric/verify_case2_sine.py`
4. `verification/Case2/verify_case2_actual_mixed_prefix.py`
5. `verification/Case2/case2_actual_mixed_prefix_data.json`
6. `verification/Case3/verify_case3_mixed.py`
7. `verification/Case3/case3_mixed_certificate_data.json`
8. `verification/Case4/verify_case4_mixed.py`
9. `verification/Case4/case4_mixed_certificate_data.json`

The circle checker embeds the fixed `(a,c)` pairs from Appendix A.2:

| Case | Left half, `[150,180]` degrees | Right half, `[180,210]` degrees |
| --- | --- | --- |
| 3 | `(-13/20, -5/16)` | `(-1/2, -1/2)` |
| 4 | `(-3/4, -1/2)` | `(-1/2, -1/4)` |

It reads no coefficient data file. The mixed checkers read the three JSON files listed above. Section 5.3 of the paper and REPRODUCIBILITY.md explain the coverage, the coefficient strings, the symbolic entries and the prefix counts.

## Mathematical dependencies and strictness

Executable acceptance certifies the stated coefficient identities and prefix inequalities. The manuscript supplies the compact-cover reduction, supporting wedges and the order of the contacts (the Wedge Lemma of Section 4), the choice of witnesses, localization, symmetry and reflection arguments, radial–radial reduction, Case 1 unfolding, and the balanced-prefix length inequality. Those arguments connect the finite checks to universality for unit arcs.

The mixed certificates prove prefix norms at most one. Analytic norm-one cases are part of the proof. Strict length greater than one follows from positive radial or cap slack with `0 < lambda < 1` and positive cap multiplier `(1 - lambda)kappa`; it does not require a uniform strict prefix-norm margin. In particular, the two prefixes of Lemma 5.5 are handled analytically under their checked hypotheses.

The Case 2 cap–cap certificate separately has the uniform bound `49/50` and the resulting `50/49` chain-length bound. Figure 6 illustrates this certificate; the plotted curves are not part of its verification.

## Acceptance and controls

The driver must exit successfully after all five stages and print `H15 PROOF CHECKS: PASS`. It rejects optimized Python, because optimization can disable assertions needed by the complete route. Do not set `PYTHONOPTIMIZE` or pass `-O` or `-OO`.

`test_geometric_certificates.py` has three positive and twenty-three rejection controls, including the circle checker’s refusal of optimized execution. `test_mixed_certificates.py` has twelve rejection controls: a missing rectangle, an angular coverage gap, a negative weight, and a coefficient with the wrong sign in each mixed case. Controls alter temporary copies, require the intended failure reason, and leave distributed inputs unchanged. They support validation of the checkers and do not replace the proof run.

The execution logs and their machine-readable summary are provided in `validation/`, relative to the supplement root. These records identify the Python version and report the proof and control results. Their elapsed times are descriptive, not acceptance thresholds.

The inequalities of Section 6.2 (Proposition 6.6 and Remark 6.8) are not part of the proof of Theorem 1.1. They have their own checker, `Section6/verify_zarc_test.py`, which reads no data file.

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for commands, portability, and detailed acceptance criteria.
