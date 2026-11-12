#!/usr/bin/env python3
"""An independent re-check of the computations behind Theorem 1.1.

This program was written from the statements of the paper (Lemma 5.4, Lemma 5.5,
Proposition 5.2, Proposition 5.3 and Appendix A), not from the checkers in
verification/Case2, Case3, Case4 and Geometric, and it shares no code with them.
It is not part of the proof route; it is a second opinion.

Part A repeats the interval checks of subcase R_U C_V with the interval arithmetic of
mpmath instead of the rational intervals of the checkers.  It reads the three coefficient
files, checks that the rectangles cover [150,210] x [0,30], and bounds the squared norm of
every proper prefix on every rectangle, bisecting as the checkers do.  It must reach the
counts of Table 3 of the paper: 914 pairs settled by interval bounds and 48 by Lemma 5.5.

Part B uses floating-point numbers and proves nothing.  It samples the same prefix norms,
tests the value of Phi in Lemma 5.4, Proposition 5.2 and Appendix A.2 on random chains
whose contacts lie on the rays of Table 2, and recomputes the numbers of Appendix A.1,
Figure 6 and Table 5.

Run from the supplement root:   python3 verification/Independent/recheck_with_mpmath.py
Needs NumPy and mpmath.  Takes about half a minute.
"""
import datetime, itertools, json, math, platform, random, sys
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import mpmath
from mpmath import iv

SUPP = sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).resolve().parent.parent.parent)
iv.prec = 120
print("Independent re-check with mpmath intervals and floating-point sampling")
print("UTC:", datetime.datetime.now(datetime.timezone.utc).isoformat())
print("Python:", sys.version.split()[0], "| NumPy", np.__version__, "| mpmath", mpmath.__version__)
print("Platform:", platform.platform())


def ivF(x):
    """rigorous interval for a rational number"""
    x = F(x)
    return iv.mpf(x.numerator) / iv.mpf(x.denominator)

def P(lo, hi=None):
    """interval enclosure of (cos, sin) over an angle interval in degrees"""
    lo_i = ivF(lo); hi_i = ivF(hi if hi is not None else lo)
    t = iv.mpf([lo_i.a, hi_i.b]) * iv.pi / 180
    return (iv.cos(t), iv.sin(t))

def add(u, v): return (u[0] + v[0], u[1] + v[1])
def sc(s, u): return (s * u[0], s * u[1])
def n2hi(u):
    return float(((u[0] ** 2) + (u[1] ** 2)).b)

# orders (10) of the paper
ORDERS = []
for pre in itertools.permutations(["A", "Xm", "Xp"]):
    ORDERS.append(list(pre) + ["B", "C", "D", "E", "Y", "F"])
    ORDERS.append(list(pre) + ["B", "C", "D", "E", "F", "Y"])
PREFIXES = sorted({frozenset(o[:k]) for o in ORDERS for k in range(1, 9)},
                  key=lambda s: (len(s), sorted(s)))
assert len(ORDERS) == 12 and len(PREFIXES) == 13

def coeffs_iv(case, th, ph, lam, a, c):
    """Lemma 5.4 coefficient table, interval version. th, ph = (lo, hi) degrees."""
    eps = 1 if case == 3 else 0
    T, Pp = ("D", "C") if case == 2 else ("C", "D")
    lam = ivF(lam); c = ivF(c)
    s3 = iv.sqrt(3)
    a = -s3 * lam if a is None else ivF(a)
    cph, sph = P(*ph)
    mu0 = -lam * (sph + s3 * cph)
    mu1 = -2 * lam * cph
    g = -2 * (1 - lam) * P(15)[0]
    kap = 1 / (2 * P(F(15, 2))[0])
    e0 = P(90); e1 = P(-60)
    d = [P(th[0] + 90, th[1] + 90), P(th[0] - 60, th[1] - 60)]
    co = {"A": sc(a, e0), "B": sc(mu1, e1), T: sc(mu0 - a, e0),
          "Y": sc(lam, (cph, sph)),
          Pp: sc(c, d[eps]), "E": sc(g, d[1 - eps]), "F": sc(g - c, d[eps]),
          "Xm": sc((1 - lam) * kap, P(th[0] + F(15, 2), th[1] + F(15, 2))),
          "Xp": sc((1 - lam) * kap, P(th[0] + F(45, 2), th[1] + F(45, 2)))}
    return co

def certify(case, labels, th, ph, lam, a, c, depth=0, stat=None):
    co = coeffs_iv(case, th, ph, lam, a, c)
    v = (iv.mpf(0), iv.mpf(0))
    for l in labels: v = add(v, co[l])
    stat["nodes"] += 1; stat["depth"] = max(stat["depth"], depth)
    if n2hi(v) <= 1: return
    assert depth < 34, ("undecided", case, sorted(labels), th, ph)
    if th[1] - th[0] >= ph[1] - ph[0]:
        m = (th[0] + th[1]) / 2
        certify(case, labels, (th[0], m), ph, lam, a, c, depth + 1, stat)
        certify(case, labels, (m, th[1]), ph, lam, a, c, depth + 1, stat)
    else:
        m = (ph[0] + ph[1]) / 2
        certify(case, labels, th, (ph[0], m), lam, a, c, depth + 1, stat)
        certify(case, labels, th, (m, ph[1]), lam, a, c, depth + 1, stat)

def load(case):
    f = {2: "Case2/case2_actual_mixed_prefix_data.json",
         3: "Case3/case3_mixed_certificate_data.json",
         4: "Case4/case4_mixed_certificate_data.json"}[case]
    d = json.load(open(f"{SUPP}/verification/{f}"))
    rows = d["boxes"] if case == 2 else d["mixed_boxes"]
    out = []
    for r in rows:
        q = lambda s: F(str(s))
        a = None if r["a"] == "-sqrt(3)*lambda" else q(r["a"])
        out.append((r["name"], tuple(map(q, r["theta"])), tuple(map(q, r["phi"])),
                    q(r["lambda"]), a, q(r["c"])))
    return out

def coverage(rows):
    """exact: union of rectangles = [150,210] x [0,30] (independent method: inclusion of
    every cell of the grid cut by all end points, tested at the cell's midpoint and corners)"""
    tx = sorted({F(150), F(210)} | {x for r in rows for x in r[1]})
    px = sorted({F(0), F(30)} | {x for r in rows for x in r[2]})
    assert tx[0] == 150 and tx[-1] == 210 and px[0] == 0 and px[-1] == 30
    area = F(0)
    for t0, t1 in zip(tx, tx[1:]):
        for p0, p1 in zip(px, px[1:]):
            ok = any(r[1][0] <= t0 and t1 <= r[1][1] and r[2][0] <= p0 and p1 <= r[2][1] for r in rows)
            assert ok, ("gap", t0, t1, p0, p1)
            area += (t1 - t0) * (p1 - p0)
    assert area == 60 * 30
    overlap = sum((r[1][1] - r[1][0]) * (r[2][1] - r[2][0]) for r in rows) - area
    return len(tx), len(px), overlap

LEMMA55_i = frozenset(["A", "B", "Xm", "Xp"])
LEMMA55_ii = frozenset(["A", "B", "C", "Xm", "Xp"])

print("PART A: rigorous re-check with mpmath intervals (prec %d bits)" % iv.prec)
tot_iv = tot_an = 0
for case in (2, 3, 4):
    rows = load(case)
    nt, npx, ov = coverage(rows)
    stat = {"nodes": 0, "depth": 0}
    niv = nan = 0
    for name, th, ph, lam, a, c in rows:
        assert 150 <= th[0] < th[1] <= 210 and 0 <= ph[0] < ph[1] <= 30 and 0 < lam < 1
        for labels in PREFIXES:
            if a is None and labels == LEMMA55_i:
                assert ph[1] <= 60; nan += 1; continue      # Lemma 5.5(i)
            if case in (3, 4) and labels == LEMMA55_ii:
                nan += 1; continue                           # Lemma 5.5(ii)
            certify(case, labels, th, ph, lam, a, c, 0, stat); niv += 1
    tot_iv += niv; tot_an += nan
    print(f"  Case {case}: {len(rows)} rectangles cover the domain exactly "
          f"(grid {nt-1}x{npx-1}, overlap area {float(ov):.3g}); "
          f"interval {niv}, Lemma 5.5 {nan}, nodes {stat['nodes']}, depth {stat['depth']}")
print(f"  total: interval {tot_iv}, Lemma 5.5 {tot_an}, pairs {tot_iv + tot_an}")
assert (tot_iv, tot_an) == (914, 48)

# ---------------------------------------------------------------- Part B
print("PART B: floating-point checks")
rad = math.pi / 180
p = lambda t: np.array([math.cos(t * rad), math.sin(t * rad)])
R = lambda t: np.array([[math.cos(t * rad), -math.sin(t * rad)], [math.sin(t * rad), math.cos(t * rad)]])
e0, e1 = p(90), p(-60)
rho = math.cos(7.5 * rad); kap = 1 / (2 * rho); b = 2 * math.cos(15 * rad) - 1
assert np.allclose(kap * p(22.5), e0 + b * e1) and np.allclose(kap * p(7.5), b * e0 + e1)
print("  (19): kappa p_22.5 = e0 + b e1, kappa p_7.5 = b e0 + e1  ok;  b = %.6f, kappa = %.6f" % (b, kap))

def coeffs(case, th, ph, lam, a, c):
    eps = 1 if case == 3 else 0
    T, Pp = ("D", "C") if case == 2 else ("C", "D")
    if a is None: a = -math.sqrt(3) * lam
    mu0 = -lam * (math.sin(ph * rad) + math.sqrt(3) * math.cos(ph * rad))
    mu1 = -2 * lam * math.cos(ph * rad)
    g = -2 * (1 - lam) * math.cos(15 * rad)
    d = [R(th) @ e0, R(th) @ e1]
    return {"A": a * e0, "B": mu1 * e1, T: (mu0 - a) * e0, "Y": lam * p(ph),
            Pp: c * d[eps], "E": g * d[1 - eps], "F": (g - c) * d[eps],
            "Xm": (1 - lam) * kap * (R(th) @ p(7.5)), "Xp": (1 - lam) * kap * (R(th) @ p(22.5))}

for case in (2, 3, 4):
    rows = load(case); worst = (0, None); wnon = (0, None); zs = 0
    for name, th, ph, lam, a, c in rows:
        lam, c = float(lam), float(c); a = None if a is None else float(a)
        for t in np.linspace(float(th[0]), float(th[1]), 41):
            for f in np.linspace(float(ph[0]), float(ph[1]), 41):
                co = coeffs(case, t, f, lam, a, c)
                zs = max(zs, np.linalg.norm(sum(co.values())))
                for labels in PREFIXES:
                    n = np.linalg.norm(sum(co[l] for l in labels))
                    if n > worst[0]: worst = (n, (name, sorted(labels), round(float(t), 3), round(float(f), 3)))
                    lem = (a is None and labels == LEMMA55_i) or (case != 2 and labels == LEMMA55_ii)
                    if not lem and n > wnon[0]: wnon = (n, (name, sorted(labels), round(float(t), 3), round(float(f), 3)))
    assert worst[0] <= 1 + 1e-12 and zs < 1e-12
    print(f"  Case {case}: sampled max prefix norm {worst[0]:.9f} at {worst[1]}")
    print(f"           max over interval-checked prefixes {wnon[0]:.9f} at {wnon[1]}; |sum of coefficients| <= {zs:.1e}")

# Phi identity of Lemma 5.4 on random configurations
random.seed(1)
for case in (2, 3, 4):
    eps = 1 if case == 3 else 0
    T, Pp = ("D", "C") if case == 2 else ("C", "D")
    err = 0
    for _ in range(2000):
        th = random.uniform(150, 210); ph = random.uniform(0, 30); lam = random.uniform(.05, .95)
        a = random.uniform(-2, 2); c = random.uniform(-2, 2)
        V = np.array([random.uniform(-2, 2), random.uniform(-2, 2)])
        ray = {0: R(th) @ p(0), 1: R(th) @ p(30)}       # V_0, V_1 directions
        pts = {"A": random.uniform(0, 2) * p(0), T: random.uniform(0, 2) * p(0), "B": random.uniform(0, 2) * p(30),
               Pp: V + random.uniform(0, 2) * ray[eps], "F": V + random.uniform(0, 2) * ray[eps],
               "E": V + random.uniform(0, 2) * ray[1 - eps]}
        r = random.uniform(1, 1.5); pts["Y"] = r * p(ph)
        pts["Xm"] = V + np.array([random.uniform(-1, 1), random.uniform(-1, 1)])
        pts["Xp"] = V + np.array([random.uniform(-1, 1), random.uniform(-1, 1)])
        co = coeffs(case, th, ph, lam, a, c)
        Phi = sum(co[l] @ pts[l] for l in co)
        sm = (pts["Xm"] - V) @ (R(th) @ p(7.5)) - rho; sp = (pts["Xp"] - V) @ (R(th) @ p(22.5)) - rho
        err = max(err, abs(Phi - (1 + lam * (r - 1) + (1 - lam) * kap * (sm + sp))))
    assert err < 1e-12
    print(f"  Case {case}: Phi = 1 + lam(r-1) + (1-lam) kappa (s- + s+) on 2000 random chains, error {err:.1e}")
print("  incidences used: Case 2 U0:A,D U1:B V0:C,F V1:E; Case 3 U0:A,C U1:B V0:E V1:D,F; Case 4 U0:A,C U1:B V0:D,F V1:E (Table 2)")

# Lemma 5.5
mx = 0
for _ in range(20000):
    ph = random.uniform(0, 30); th = random.uniform(150, 210); lam = random.uniform(0, 1)
    w = np.array([-math.cos(ph * rad), math.sqrt(3) * (math.cos(ph * rad) - 1)])
    co = coeffs(3, th, ph, lam, None, 0.3)
    assert np.allclose(co["A"] + co["B"] + co["Xm"] + co["Xp"], lam * w + (1 - lam) * p(th + 15))
    assert np.allclose(co["A"] + co["B"] + co["C"] + co["Xm"] + co["Xp"], -lam * p(ph) + (1 - lam) * p(th + 15))
    mx = max(mx, np.linalg.norm(w))
print(f"  Lemma 5.5: both prefix formulas hold; max |w| on [0,30] = {mx:.6f} <= 1")

# Proposition 5.2 / Appendix A.1
Q = e0 + b * e1
def capcap2(th):
    x = math.sin(th * rad); a = -(2 + x) / 3; c = -(1 + x) / 3; Rt = R(th)
    return {"A": a * e0, "B": -b * e1, "D": (-1 - a) * e0, "XU": Q,
            "C": c * (Rt @ e0), "E": -b * (Rt @ e1), "F": (-1 - c) * (Rt @ e0), "XV": Rt @ Q}
O4 = [list(l) + list("BCDE") + list(r) for l in (["A", "XV"], ["XV", "A"]) for r in (["XU", "F"], ["F", "XU"])]
PF4 = sorted({frozenset(o[:k]) for o in O4 for k in range(1, 8)}, key=lambda s: (len(s), sorted(s)))
assert len(PF4) == 9
mx = {}
for th in np.linspace(150, 210, 6001):
    co = capcap2(th)
    assert np.linalg.norm(sum(co.values())) < 1e-12
    for s in PF4:
        n = np.linalg.norm(sum(co[l] for l in s)); k = "".join(sorted(s))
        if n > mx.get(k, (0, 0))[0]: mx[k] = (n, th)
print("  Prop 5.2 (Case 2, both sectors): max norm of each of the 9 prefixes over [150,210]:")
for k, (n, th) in sorted(mx.items(), key=lambda kv: (len(kv[0]), kv[0])):
    print(f"      {k:10s} {n:.5f} at {th:.2f}")
assert max(n for n, _ in mx.values()) < 49 / 50
err = 0
for _ in range(2000):
    th = random.uniform(150, 210); V = np.array([random.uniform(-2, 2), random.uniform(-2, 2)]); Rt = R(th)
    pts = {"A": random.uniform(0, 2) * p(0), "D": random.uniform(0, 2) * p(0), "B": random.uniform(0, 2) * p(30),
           "C": V + random.uniform(0, 2) * (Rt @ p(0)), "F": V + random.uniform(0, 2) * (Rt @ p(0)),
           "E": V + random.uniform(0, 2) * (Rt @ p(30)),
           "XU": np.array([random.uniform(-1, 1), random.uniform(-1, 1)]),
           "XV": V + np.array([random.uniform(-1, 1), random.uniform(-1, 1)])}
    co = capcap2(th); Phi = sum(co[l] @ pts[l] for l in co)
    sU = pts["XU"] @ p(22.5) - rho; sV = (pts["XV"] - V) @ (Rt @ p(22.5)) - rho
    err = max(err, abs(Phi - (1 + kap * (sU + sV))))
assert err < 1e-12
print(f"      Phi = 1 + kappa (s_U + s_V) on 2000 random chains, error {err:.1e}")
q1 = b / 2; q2 = 1 - math.sqrt(3) * q1; k = q1 * q1 + (q2 - 1 / 3) ** 2
print(f"      q1 = {q1:.6f}, q2 = {q2:.6f}, 4k = {4*k:.6f} (< 19/20), |Q|^2 = {q1*q1+q2*q2:.6f} (< 13/50)")
A1 = 1/9 - q2/2 - 2*q1/3; B1 = 2*q2/3 - 4*q1/3 + 4/9; C1 = q1*q1 + (q2 + 2/3)**2
A2 = -q2*q2 + q2/2 - q1*q1 - 2*q1/3 + 1/9; B2 = 4*q2/3 + 2*q1/3 - 2/9; C2 = 4*q2*q2 - 4*q2/3 + 4*q1*q1 + 1/9
print(f"      A1,B1,C1 = {A1:.5f}, {B1:.5f}, {C1:.5f}  (table: -A1>=0.29, |B1|<=0.05, C1<=0.957)")
print(f"      A2,B2,C2 = {A2:.5f}, {B2:.5f}, {C2:.5f}  (table: -A2>=0.35, |B2|<=0.35, C2<=0.872)")
assert -A1 >= .29 and abs(B1) <= .05 and C1 <= .957 and -A2 >= .35 and abs(B2) <= .35 and C2 <= .872
assert F(957, 1000) + F(1, 464) == F(55631, 58000) and F(109, 125) + F(7, 80) == F(1919, 2000)
assert F(55631, 58000) < F(49, 50) ** 2 and F(1919, 2000) < F(49, 50) ** 2

# Appendix A.2 / Table 5
print("  Table 5 (Cases 3 and 4, both sectors): largest norms and where")
h = {22.5: (-1, -b), 7.5: (-b, -1)}
tab = {(3, 0): (-13/20, -5/16), (3, 1): (-1/2, -1/2), (4, 0): (-3/4, -1/2), (4, 1): (-1/2, -1/4)}
E = [e0, e1]
for case in (3, 4):
    eps = 1 if case == 3 else 0
    for half, (lo, hi) in enumerate([(150, 180), (180, 210)]):
        fU, fV = (22.5, 7.5) if half == 0 else (7.5, 22.5)
        a, c = tab[(case, half)]
        def co34(th):
            Rt = R(th)
            return {"A": a * e0, "B": h[fU][1] * e1, "C": (h[fU][0] - a) * e0, "XU": kap * p(fU),
                    "D": c * (Rt @ E[eps]), "E": h[fV][1 - eps] * (Rt @ E[1 - eps]),
                    "F": (h[fV][eps] - c) * (Rt @ E[eps]), "XV": kap * p(th + fV)}
        res = {}
        for th in np.linspace(lo, hi, 30001):
            co = co34(th); assert np.linalg.norm(sum(co.values())) < 1e-12
            for s in PF4:
                n = np.linalg.norm(sum(co[l] for l in s)); kk = "".join(sorted(s))
                if n > res.get(kk, (0, 0))[0]: res[kk] = (n, th)
        row = [res[kk] for kk in ("AXV", "ABXV", "ABCDXV", "ABCDEXV")]
        crit = res["ABCXV"]; const = max(res[kk][0] for kk in res if kk not in ("AXV", "ABXV", "ABCXV", "ABCDXV", "ABCDEXV"))
        print(f"      Case {case} [{lo},{hi}]: " + "  ".join(f"{n:.6f} ({th:.2f})" for n, th in row)
              + f" | critical {crit[0]:.6f} ({crit[1]:.1f}) | constant <= {const:.4f}")
        assert max(n for n, _ in row) < .999 and crit[0] <= 1 + 1e-12 and const <= .75 + 1e-12
        # Phi identity
        err = 0
        for _ in range(500):
            th = random.uniform(lo, hi); V = np.array([random.uniform(-2, 2), random.uniform(-2, 2)]); Rt = R(th)
            ray = {0: Rt @ p(0), 1: Rt @ p(30)}
            pts = {"A": random.uniform(0, 2) * p(0), "C": random.uniform(0, 2) * p(0), "B": random.uniform(0, 2) * p(30),
                   "D": V + random.uniform(0, 2) * ray[eps], "F": V + random.uniform(0, 2) * ray[eps],
                   "E": V + random.uniform(0, 2) * ray[1 - eps],
                   "XU": np.array([random.uniform(-1, 1), random.uniform(-1, 1)]),
                   "XV": V + np.array([random.uniform(-1, 1), random.uniform(-1, 1)])}
            co = co34(th); Phi = sum(co[l] @ pts[l] for l in co)
            sU = pts["XU"] @ p(fU) - rho; sV = (pts["XV"] - V) @ p(th + fV) - rho
            err = max(err, abs(Phi - (1 + kap * (sU + sV))))
        assert err < 1e-12

# Lemma 3.3 and scalars quoted in the text
for Th in (60, 90, 120):
    w = math.sin(Th / 2 * rad) / math.cos(15 * rad); w2 = math.cos((Th / 2 - 15) * rad) / math.cos(15 * rad)
    cA, cX, cY, cF = -w2 * p(90), w * p(15), w * p(Th - 15), w2 * p(90 + Th)
    assert np.linalg.norm(cA + cX + cY + cF) < 1e-12 and abs(np.linalg.norm(cA + cX) - 1) < 1e-12
    assert w2 <= 1 + 1e-12 and w < 1
print("  Lemma 3.3: coefficients sum to zero and |c_A + c_X| = 1 for Theta = 60, 90, 120")
area = math.pi / 24 + (math.sqrt(6) - math.sqrt(2)) / 8
print(f"  area(H15) = {area:.10f}; sector - H15 = {math.pi/12 - area:.6f}; triangle - H15 = {0.2609562398 - area:.6f}; "
      f"d/dt = {math.pi/360*(1-math.cos(15*rad)):.4e}; 1 - rho = {1-rho:.6f}; 1 + 2 sin 1 = {1+2*math.sin(rad):.6f}")
print("INDEPENDENT RE-CHECK: PASS")
