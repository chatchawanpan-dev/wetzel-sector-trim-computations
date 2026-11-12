# Exploratory computations

The two programs in this folder are the numerical experiments described in Section 1.1 of the paper. No statement of the paper depends on them. They need NumPy and SciPy (the first) or NumPy, CVXPY and Clarabel (the second); the exact checkers in `verification/` need none of these.

## 1. Screening the trims `H_t`

`screen_cap_branches.py` works with a simplified model. A chain visits three points alternately on the lower, upper and lower radius of one unit sector, three more on the upper, lower and upper radius, and one point of each closed cap of `H_t`. The program runs through the 1120 orders of these eight points that keep the order within each triple, and for each order it minimizes the length of the chain by L-BFGS-B from one or more fixed starting points. The values it reports are local minima, not certified bounds, and the solver often stops with a warning, because the length is not differentiable where points coincide.

```sh
python3 exploration/screen_cap_branches.py --candidate H16 --restarts 1 --maxiter 300
```

A candidate is written `H<t>`, for example `H15` or `H15.5`. One candidate takes a minute or two. The complete table is written to `exploration/output/`.

`screening_results.txt` records the smallest length found for the values of `t` quoted in the paper, from a run on 11 November 2026 with Python 3.13.18, NumPy 2.5.3 and SciPy 1.18.1:

| `t` | smallest length found | chain |
| ---: | ---: | --- |
| 12° | 0.993129 | `O`, then a common point of the two caps |
| 14° | 0.997876 | `O`, then a common point of the two caps |
| 15° | 1.000000 | `O`, `p_15`, `p_15` |
| 15.5° | 1.017453 | `O`, `p_14.5`, `p_15.5` |
| 16° | 1.034905 | `O`, `p_14`, `p_16` |
| 17° | 1.069799 | `O`, `p_13`, `p_17` |
| 20° | 1.174311 | `O`, `p_10`, `p_20` |

For `t ≥ 15°` these are the lengths `1 + 2 sin(t − 15°)`. The first run of this experiment, in the summer of 2026 with SciPy 1.17.1, was for `t = 15°` and `16°` only and gave the same two minima.

## 2. Searching for the Case 2 coefficients

`find_case2_coefficients.py` looks for the coefficients `lambda`, `a`, `c` of Section 5.3 in Case 2. On a rectangle of the `(theta, phi)` domain it solves a second-order cone program that minimizes the largest prefix norm at the nine points of a 3 × 3 grid. If the optimal value is below 0.9995, it rounds the solution to nine decimal places and passes it to the exact checker `verification/Case2/verify_case2_actual_mixed_prefix.py`; the rectangle is kept only if every prefix passes on the whole rectangle. Otherwise the rectangle is bisected. The search starts from sixteen rectangles.

```sh
python3 exploration/find_case2_coefficients.py
python3 verification/Case2/verify_case2_actual_mixed_prefix.py exploration/output/case2_coefficients_rerun.json
```

`case2_search_2026-09-06.log` is the log of the run that produced the coefficients distributed in `verification/Case2/case2_actual_mixed_prefix_data.json`: 33 rectangles after 17 bisections, with Python 3.12, NumPy 2.5.2, CVXPY 1.9.2 and Clarabel 0.11.1. A run on 11 November 2026 with CVXPY 1.9.3 gave the same 33 rectangles, with coefficients that differ from the distributed ones by less than 0.00003, and the exact checker accepted them. The coefficients of a new run depend on the versions of the numerical software, so they need not agree with the distributed ones digit for digit. The program never writes to the distributed coefficient file.

The coefficients for Cases 3 and 4 were found by numerical optimization in the same spirit; no search program for them is included.
