# H15 verification supplement

Supplement to *A One-Sided Chord Trim of Wetzel's Sector Covers Unit Arcs*, by Chatchawan Panraksa and Wacharin Wichiramala, submitted to Experimental Mathematics in November 2026.

The archive holds the coefficient data and exact checkers for the computer-assisted steps of the proof, a checker for the inequalities of Section 6.2, control tests, an independent re-check, execution records, the two programs for the experiments of Section 1.1, and the program that draws the figures. The checkers need only Python 3 and its standard library. They use no optimizer, coefficient search or network access.

## Checking the proof

From this folder, run

```sh
python3 verification/verify_proof_essential.py
```

A complete run passes five stages, exits with status zero and prints `H15 PROOF CHECKS: PASS`. Do not use `python3 -O` or set `PYTHONOPTIMIZE`, because the checkers rely on assertions. On a laptop the run takes about 15 seconds with Python 3.12 and about a minute with Python 3.9.

| Stage | Program (in `verification/`) | Part of the paper |
| --- | --- | --- |
| 1 | `Geometric/verify_capcap_circles.py` | Proposition 5.3 and Appendix A.2: Cases 3 and 4, both sectors contain the arc |
| 2 | `Geometric/verify_case2_sine.py` | Proposition 5.2 and Appendix A.1: Case 2, both sectors contain the arc |
| 3 | `Case2/verify_case2_actual_mixed_prefix.py` with `Case2/case2_actual_mixed_prefix_data.json` | Lemmas 5.4 and 5.5, Table 3 and Proposition 5.6, Case 2 |
| 4 | `Case3/verify_case3_mixed.py` with `Case3/case3_mixed_certificate_data.json` | the same, Case 3 |
| 5 | `Case4/verify_case4_mixed.py` with `Case4/case4_mixed_certificate_data.json` | the same, Case 4 |

Each row of the three data files gives two angle intervals and the coefficients lambda, a and c of one rectangle. [PROOF_ESSENTIAL.md](verification/PROOF_ESSENTIAL.md) describes the certificates. [REPRODUCIBILITY.md](verification/REPRODUCIBILITY.md) gives the individual commands, the parameter domains, the arithmetic and the acceptance criteria. [PROOF_ESSENTIAL_FILES.json](verification/PROOF_ESSENTIAL_FILES.json) lists the nine files the proof run uses, with their sizes and SHA-256 hashes.

## Checking Section 6.2

```sh
python3 verification/Section6/verify_zarc_test.py
```

This program is separate from the proof of Theorem 1.1. It checks, with rational interval arithmetic, the twelve inequalities behind Table 4 and the inequality for the point `Q_3` (Proposition 6.6), and the nine inequalities for the pentagon (Remark 6.8). It prints `SECTION 6 CHECKS: PASS` and takes about a second.

## Controls and records

```sh
python3 verification/test_geometric_certificates.py
python3 verification/test_mixed_certificates.py
python3 verification/test_section6.py
python3 verify_snapshot.py
```

The control tests confirm that the checkers accept the supplied certificates and reject altered copies for the intended reason: a missing rectangle, a gap in the cover, a negative weight, and a coefficient with the wrong sign, among others. `validation/` holds the logs of complete runs, and [Environment.json](Environment.json) records the environments. `verify_snapshot.py` compares the files with `SHA256SUMS.json`; this is an integrity check, not a proof.

## A second check

```sh
python3 verification/Independent/recheck_with_mpmath.py
```

This program was written from the statements of the paper and shares no code with the checkers. It repeats the interval checks of the three mixed cases with the interval arithmetic of mpmath and must arrive at the counts of Table 3 (914 pairs settled by interval bounds, 48 by Lemma 5.5). With floating-point numbers it also samples the prefix norms, tests the value of `Phi` on random chains whose contacts lie on the rays of Table 2, and recomputes the numbers of Appendix A.1, Figure 6 and Table 5. It needs NumPy and mpmath, takes about half a minute and prints `INDEPENDENT RE-CHECK: PASS`. It is not part of the proof route; its output is in `validation/independent_recheck.log`.

## Experiments

The folder `exploration/` holds the screening program and the coefficient search described in Section 1.1, with a short guide and the log of the search that produced the Case 2 coefficients. They need NumPy, SciPy, CVXPY and Clarabel. No statement of the paper depends on them.

## Figures

`figures/make_figures.py` draws Figures 1–8 of the paper from the formulas in the text, as PDF and EPS. It needs NumPy and Matplotlib; labels are set with Matplotlib's own Computer Modern fonts, so no LaTeX installation is needed. Figure 7(b) reads the three coefficient files. The program checks the configurations it draws: in Figure 5 that `X` and `Y` are the points of the chain farthest from `V` and from `U`, and in Figure 7(a) that the chain lies in both wedges and in the unit sector at `V`, that `Y` is its point farthest from `U`, and that the two witnesses maximize the cap functionals. The figures are illustrations, and no step of the proof depends on them.

Colors have fixed meanings: blue for the wedge at `U` and what is measured from `U`, orange for the wedge at `V`, green for what is cut away or singled out (the caps, the arc `J` of Figure 4, the part of `H_15` outside `K_0`). Every colored mark also carries a label, a line style or a hatching, so the figures can be read in grayscale. The three hues (`#2a78d6`, `#eb6834`, `#1baf7a`) are separated under simulated protanopia and deuteranopia (smallest pairwise distance 9.2 in OKLab x 100, Machado–Oliveira–Fernandes model).

## Scope

The programs check the finite identities, the coverage of the parameter domains, the prefix-norm bounds and the inequalities of Section 6.2 stated in the paper. The geometric reductions are proved in the paper: the supporting wedges (the Wedge Lemma of Section 4), the localization of the witnesses, the reflections and the symmetries, and Lemma 6.5.
