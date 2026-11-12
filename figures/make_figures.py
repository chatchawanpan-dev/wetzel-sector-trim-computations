"""Draw the eight figures of the manuscript.

Run:  python3 make_figures.py          (needs NumPy and Matplotlib; no LaTeX installation)
Writes fig1_trim ... fig8_example as PDF and EPS next to this script, and PNG previews in
preview/.  All coordinates come from the formulas in the paper; only the curves of Figures
1(b,c) and 6 are sampled.  Figure 7(b) reads the three coefficient files in ../verification.

Colors have fixed meanings throughout:
    blue    the wedge at U and what is measured from U
    orange  the wedge at V and what is measured from V
    green   what is cut away or singled out (the caps, the arc J, the set H_15 minus K_0)
Every colored mark also has a label, a line style or a hatching, so the figures can be read
in grayscale.  The three hues were checked for color-vision deficiencies (see README).
"""
import json
import logging
import os
from fractions import Fraction

import numpy as np
import matplotlib
matplotlib.use("pdf")
logging.getLogger("fontTools").setLevel(logging.ERROR)
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb
from matplotlib.patches import Polygon, Arc, Rectangle, FancyArrowPatch

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.serif": ["cmr10"],
    "mathtext.fontset": "cm",
    "axes.formatter.use_mathtext": True,
    "axes.unicode_minus": False,
    "font.size": 9,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "lines.solid_capstyle": "round",
    "lines.dash_capstyle": "butt",
    "hatch.linewidth": 0.5,
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})
TW = 5.6            # width of a full-width figure in inches
LW = 0.8            # boundaries
LT = 0.5            # thin construction lines
LC = 1.3            # chains (pieces of an arc)
HERE = os.path.dirname(os.path.abspath(__file__))


def tint(color, weight):
    """mix a color with white; weight = share of the color"""
    r, g, b = to_rgb(color)
    return (1 - weight + weight * r, 1 - weight + weight * g, 1 - weight + weight * b)


INK = "#0b0b0b"
BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#1baf7a"
BLUE_T, ORANGE_T, GREEN_T = tint(BLUE, 0.22), tint(ORANGE, 0.22), tint(GREEN, 0.30)
GREY = "#e4e3df"     # a set that is kept (H_15, K_0, conv gamma)
GREY2 = "#f1f0ed"
MUTED = "#6b6a66"


def deg(x):
    """a number of degrees as a math-text label"""
    return r"$%s\!^\circ$" % x


def p(t):
    t = np.deg2rad(t)
    return np.array([np.cos(t), np.sin(t)])


def arcpts(t0, t1, r=1.0, c=(0.0, 0.0), n=200):
    return np.array([np.asarray(c) + r * p(t) for t in np.linspace(t0, t1, n)])


def seg(ax, P, Q, **kw):
    kw.setdefault("color", INK)
    kw.setdefault("lw", LW)
    ax.plot([P[0], Q[0]], [P[1], Q[1]], **kw)


def path(ax, pts, **kw):
    pts = np.asarray(pts)
    kw.setdefault("color", INK)
    kw.setdefault("lw", LC)
    kw.setdefault("solid_joinstyle", "round")
    ax.plot(pts[:, 0], pts[:, 1], **kw)


def dot(ax, P, s=3.0, fc=INK, hollow=False):
    ax.plot([P[0]], [P[1]], "o", ms=s, mfc="white" if hollow else fc, mec=INK, mew=0.5, zorder=6)


def lab(ax, P, text, dx=0, dy=0, ha="center", va="center", size=None, **kw):
    if size is not None:
        kw["fontsize"] = size
    ax.annotate(text, xy=P, xytext=(dx, dy), textcoords="offset points", ha=ha, va=va, zorder=7, **kw)


def arrow(ax, P, Q, color=INK, lw=LT, head=3.0):
    ax.add_patch(FancyArrowPatch(P, Q, arrowstyle="-|>", mutation_scale=head, lw=lw, color=color,
                                 shrinkA=0, shrinkB=0, zorder=5))


def plabel(ax, text, x=0.0, y=1.0):
    ax.text(x, y, text, transform=ax.transAxes, ha="left", va="top")


def geo(ax, xlim, ylim):
    ax.set_aspect("equal")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.axis("off")


def save(fig, name):
    os.makedirs(os.path.join(HERE, "preview"), exist_ok=True)
    fig.savefig(os.path.join(HERE, name + ".pdf"))
    fig.savefig(os.path.join(HERE, name + ".eps"))
    fig.savefig(os.path.join(HERE, "preview", name + ".png"), dpi=220)
    plt.close(fig)


O = np.zeros(2)
RHO = np.cos(np.deg2rad(7.5))
DASH = (0, (3, 2))
DOTS = (0, (1, 1.5))


# ---------------------------------------------------------------- Figure 1
def fig1():
    fig = plt.figure(figsize=(TW, 2.45))
    ax = fig.add_axes([0.0, 0.0, 0.56, 1.0])
    H = np.vstack([O, arcpts(0, 15), p(30)])
    capm = arcpts(0, 15)            # lower cap: between the chord [p0, p15] and the arc
    capp = arcpts(15, 30)           # upper cap: between the chord [p15, p30] and the arc
    ax.add_patch(Polygon(H, closed=True, fc=GREY, ec="none"))
    ax.add_patch(Polygon(capm, closed=True, fc=GREEN, ec=GREEN, lw=0.4))
    ax.add_patch(Polygon(capp, closed=True, fc=GREEN, ec=GREEN, lw=0.4))
    seg(ax, O, p(0)); path(ax, arcpts(0, 15), lw=LW); seg(ax, p(15), p(30)); seg(ax, p(30), O)
    seg(ax, p(0), p(15), lw=LT, ls=DASH)
    path(ax, arcpts(15, 30), lw=LT, ls=DASH)
    seg(ax, O, p(15), lw=LT, ls=DOTS)
    for t, s, dx, dy, ha in [(0, "$p_0$", 2, -8, "center"), (15, "$p_{15}$", 6, 0, "left"),
                             (30, "$p_{30}$", 3, 7, "center")]:
        dot(ax, p(t)); lab(ax, p(t), s, dx, dy, ha=ha)
    dot(ax, O); lab(ax, O, "$O$", -7, -6)
    lab(ax, 0.60 * p(7.5), r"$\Sigma_{15}$")
    lab(ax, np.array([0.70, 0.265]), r"$\Delta_{15}$")
    ax.add_patch(Arc(O, 0.5, 0.5, theta1=0, theta2=15, lw=LT))
    lab(ax, 0.30 * p(7.5), deg(15), size=8)
    # outward unit normals of the two chords, drawn at the midpoints of the chords
    for t, s, cap in [(7.5, r"$p_{7.5}$", r"$C^-_{15}$"), (22.5, r"$p_{22.5}$", r"$C^+_{15}$")]:
        arrow(ax, RHO * p(t), (RHO + 0.12) * p(t), head=5)
        lab(ax, (RHO + 0.12) * p(t), s, 4, 0, size=8, ha="left")
        lab(ax, (RHO + 0.03) * p(t + (3.4 if t > 15 else -3.6)), cap, 3, 0, size=8, ha="left")
    geo(ax, (-0.10, 1.30), (-0.10, 0.62))
    fig.text(0.0, 0.985, "(a)", ha="left", va="top")

    def caps(ax, t, tag):
        delta = (30 - t) / 2.0
        mp, mm = (t + 30) / 2.0, (30 - t) / 2.0
        phm = np.linspace(0, 30 - t, 300)
        php = np.linspace(t, 30, 300)
        rm = np.cos(np.deg2rad(delta)) / np.cos(np.deg2rad(phm - mm))
        rp = np.cos(np.deg2rad(delta)) / np.cos(np.deg2rad(php - mp))
        ax.fill_between(phm, rm, 1, fc=GREEN_T, ec="none", lw=0)
        ax.fill_between(php, rp, 1, fc="none", ec=GREEN, hatch="/////", lw=0)
        ax.plot(phm, rm, color=INK, lw=LT)
        ax.plot(php, rp, color=INK, lw=LT)
        ax.plot([0, 30], [1, 1], color=INK, lw=LW)
        ax.plot([15], [1], "o", ms=2.8, mfc=INK, mec=INK, zorder=5, clip_on=False)
        ax.annotate("$p_{15}$", xy=(15, 1), xytext=(0, 4), textcoords="offset points", ha="center", va="bottom", fontsize=8)
        ax.set_xlim(0, 30); ax.set_ylim(0.982, 1.0035)
        ax.set_xticks([0, 15, 30]); ax.set_xticklabels([deg(0), deg(15), deg(30)])
        ax.set_yticks([0.985, 0.99, 0.995, 1.0]); ax.set_yticklabels(["0.985", "0.990", "0.995", "1.000"])
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.text(29.5, 0.9835, tag, ha="right", va="bottom", fontsize=8)

    axb = fig.add_axes([0.655, 0.59, 0.325, 0.36])
    caps(axb, 15, r"$t=15\!^\circ$")
    box = dict(fc="white", ec="none", pad=0.6)
    axb.text(7.5, 0.9958, "$C^-_{15}$", ha="center", va="center", fontsize=8)
    axb.text(22.5, 0.9958, "$C^+_{15}$", ha="center", va="center", fontsize=8, bbox=box)
    axb.set_xticklabels([])
    axb.text(-0.5, 1.0035, "$r$", ha="right", va="bottom")
    plabel(axb, "(b)", -0.24, 1.10)
    axc = fig.add_axes([0.655, 0.14, 0.325, 0.36])
    caps(axc, 10, r"$t=10\!^\circ$")
    axc.text(4.0, 0.9945, "$C^-_{10}$", ha="center", va="center", fontsize=8)
    axc.text(26.0, 0.9945, "$C^+_{10}$", ha="center", va="center", fontsize=8, bbox=box)
    axc.set_xlabel(r"$\phi$", labelpad=1)
    axc.text(-0.5, 1.0035, "$r$", ha="right", va="bottom")
    plabel(axc, "(c)", -0.24, 1.10)
    save(fig, "fig1_trim")


# ---------------------------------------------------------------- Figure 2
def fig2():
    Theta = 90.0
    half = Theta / 2
    w = np.sin(np.deg2rad(half)) / np.cos(np.deg2rad(15))
    w2 = np.cos(np.deg2rad(half - 15)) / np.cos(np.deg2rad(15))
    fig = plt.figure(figsize=(TW, 2.55))
    ax = fig.add_axes([0.0, 0.0, 0.5, 1.0])
    for (a, b), col in [((0, 30), ORANGE_T), ((60, 90), BLUE_T)]:
        ax.add_patch(Polygon([O, p(a), p(b)], closed=True, fc=col, ec="none"))
        seg(ax, p(a), p(b), lw=LT)
    seg(ax, O, 1.12 * p(0)); seg(ax, O, 1.12 * p(90))
    seg(ax, O, p(30), lw=LT); seg(ax, O, p(60), lw=LT)
    path(ax, arcpts(0, 90), lw=LT, ls=DOTS)
    A, X, Y, F = 0.46 * p(0), 0.996 * p(7.5), 0.996 * p(82.5), 0.5 * p(90)
    path(ax, [A, X, Y, F])
    for P_, s, dx, dy, col in [(A, "$A$", 0, -8, INK), (X, "$X$", 8, 0, ORANGE), (Y, "$Y$", 0, 8, BLUE),
                               (F, "$F$", -8, 0, INK), (O, "$O$", -6, -6, INK)]:
        dot(ax, P_, fc=col); lab(ax, P_, s, dx, dy)
    ax.add_patch(Arc(O, 0.26, 0.26, theta1=0, theta2=90, lw=LT))
    lab(ax, 0.2 * p(45), r"$\Theta$", size=8)
    geo(ax, (-0.12, 1.18), (-0.12, 1.15))
    plabel(ax, "(a)", 0.0, 0.98)

    bx = fig.add_axes([0.5, 0.0, 0.5, 1.0])
    path(bx, arcpts(0, 360), lw=LT, color=MUTED)
    cA = -w2 * p(90); cX = w * p(15); cY = w * p(Theta - 15); cF = w2 * p(90 + Theta)
    S = cA + cX
    # each edge is one coefficient: black c_A or c_F, orange c_X, blue c_Y
    seg(bx, O, cA, lw=LC); seg(bx, cX, S, lw=LC)
    seg(bx, O, cX, lw=LC, color=ORANGE); seg(bx, cA, S, lw=LC, color=ORANGE)
    seg(bx, S, -cY, lw=LC, ls=DASH); seg(bx, -cF, O, lw=LC, ls=DASH)
    seg(bx, S, -cF, lw=LC, ls=DASH, color=BLUE); seg(bx, -cY, O, lw=LC, ls=DASH, color=BLUE)
    for P_, s, dx, dy in [(O, "$0$", -6, 5), (cA, "$c_A$", -9, -3), (cX, "$c_X$", 4, 7), (S, "$c_A+c_X$", 18, -7),
                          (-cY, "$-c_Y$", -12, -6), (-cF, "$-c_F$", -4, 7)]:
        dot(bx, P_); lab(bx, P_, s, dx, dy)
    geo(bx, (-1.12, 1.12), (-1.12, 1.12))
    plabel(bx, "(b)", 0.02, 0.98)
    save(fig, "fig2_endwedge")


# ---------------------------------------------------------------- Figure 3
def fig3():
    V = np.array([2.0, 0.55])
    dV1 = p(210)
    U_RAYS, V_RAYS = "U", "V"
    cases = {
        "Case 2": dict(pts={"A": (0.25, 0), "B": tuple(0.75 * p(30)), "C": (1.02, 0.55), "D": (0.72, 0),
                             "E": tuple(V + 0.55 * dV1), "F": (1.72, 0.55)},
                       hull="ADEFCB", own=dict(A="U", B="U", D="U", C="V", E="V", F="V"),
                       off={"A": (0, -8), "B": (-6, 5), "C": (0, 7), "D": (0, -8), "E": (7, -3), "F": (2, 7)}),
        "Case 3": dict(pts={"A": (0.25, 0), "B": tuple(0.75 * p(30)), "C": (0.62, 0), "D": tuple(V + 0.80 * dV1),
                             "E": (1.20, 0.55), "F": tuple(V + 0.30 * dV1)},
                       hull="ACDFEB", own=dict(A="U", B="U", C="U", D="V", E="V", F="V"),
                       off={"A": (0, -8), "B": (-6, 5), "C": (0, -8), "D": (2, -8), "E": (0, 7), "F": (7, -2)}),
        "Case 4": dict(pts={"A": (0.25, 0), "B": tuple(0.75 * p(30)), "C": (0.65, 0), "D": (1.05, 0.55),
                             "E": tuple(V + 0.60 * dV1), "F": (1.72, 0.55)},
                       hull="ACEFDB", own=dict(A="U", B="U", C="U", D="V", E="V", F="V"),
                       off={"A": (0, -8), "B": (-6, 5), "C": (0, -8), "D": (0, 7), "E": (7, -3), "F": (2, 7)}),
    }
    colour = {"U": BLUE, "V": ORANGE}
    fig = plt.figure(figsize=(TW, 2.3))
    for i, (name, cfg) in enumerate(cases.items()):
        ax = fig.add_axes([0.005 + i * 0.333, 0.33, 0.325, 0.62])
        pts = {k: np.array(v) for k, v in cfg["pts"].items()}
        ax.add_patch(Polygon([pts[c] for c in cfg["hull"]], closed=True, fc=GREY, ec="none"))
        seg(ax, O, np.array([2.02, 0]), color=BLUE); seg(ax, O, 0.92 * p(30), color=BLUE)
        seg(ax, V, np.array([0.45, 0.55]), ls=DASH, color=ORANGE); seg(ax, V, V + 0.95 * dV1, ls=DASH, color=ORANGE)
        path(ax, [pts[c] for c in "ABCDEF"], lw=1.1)
        for k_, P_ in pts.items():
            dot(ax, P_, s=3.2, fc=colour[cfg["own"][k_]]); lab(ax, P_, f"${k_}$", *cfg["off"][k_], size=8)
        dot(ax, O, s=3.2, fc=BLUE); lab(ax, O, "$U$", -6, -6, size=8)
        dot(ax, V, s=3.2, fc=ORANGE); lab(ax, V, "$V$", 6, 4, size=8)
        if i == 0:
            lab(ax, np.array([2.02, 0]), "$U_0$", -2, -7, size=7)
            lab(ax, 0.36 * p(30), "$U_1$", -6, 6, size=7)
            lab(ax, np.array([0.45, 0.55]), "$V_0$", -8, 0, size=7)
            lab(ax, V + 0.22 * dV1, "$V_1$", 6, -7, size=7)
        geo(ax, (-0.18, 2.28), (-0.22, 0.70))
        ax.text(0.5, -0.02, name, transform=ax.transAxes, ha="center", va="top")
    fig.text(0.005, 0.985, "(a)", ha="left", va="top")
    ax = fig.add_axes([0.03, 0.0, 0.94, 0.27])
    xs = dict(zip("ABCDEF", [0.10, 0.27, 0.42, 0.57, 0.72, 0.88]))
    own = dict(A=BLUE, B=BLUE, C=INK, D=INK, E=ORANGE, F=ORANGE)
    ax.plot([0, 1], [0, 0], color=INK, lw=LT)
    for k_, x in xs.items():
        ax.plot([x], [0], "o", ms=3.2, mfc=own[k_], mec=INK, mew=0.5); ax.text(x, -0.18, f"${k_}$", ha="center", va="top")
    for (a, b, s, col) in [(0.0, xs["B"], r"$X_V^{\pm},\ Y_V$", ORANGE), (xs["E"], 1.0, r"$X_U^{\pm},\ Y_U$", BLUE)]:
        ax.plot([a, a, b, b], [0.16, 0.30, 0.30, 0.16], color=col, lw=LC)
        ax.text((a + b) / 2, 0.36, s, ha="center", va="bottom")
    ax.text(1.0, -0.05, r"$\tau$", ha="left", va="center")
    ax.set_xlim(-0.01, 1.03); ax.set_ylim(-0.55, 0.75); ax.axis("off")
    fig.text(0.005, 0.29, "(b)", ha="left", va="top")
    save(fig, "fig3_configurations")


# ---------------------------------------------------------------- Figure 4
def fig4():
    fig = plt.figure(figsize=(2.6, 2.1))
    ax = fig.add_axes([0, 0, 1, 1])
    a_, b_ = 1.0, 0.78
    E = lambda t: np.array([a_ * np.cos(np.deg2rad(t)), b_ * np.sin(np.deg2rad(t))])
    bd = np.array([E(t) for t in np.linspace(0, 360, 400)])
    ax.add_patch(Polygon(bd, closed=True, fc=GREY2, ec=INK, lw=LW))
    J = np.array([E(t) for t in np.linspace(60, 215, 200)])
    path(ax, J, lw=2.6, color=GREEN)
    A, X, B, T, G = E(215), E(150), E(60), E(320), E(268)
    path(ax, [A, np.array([-0.38, -0.30]), G, np.array([0.32, 0.02]), B], lw=LC, color=BLUE)
    path(ax, [X, np.array([-0.15, 0.34]), np.array([0.12, -0.18]), T], lw=LC, ls=DASH, color=ORANGE)
    for P_, s, dx, dy in [(A, "$A$", -7, -5), (X, "$X$", -8, 2), (B, "$B$", 4, 6), (T, "$T$", 6, -5)]:
        dot(ax, P_); lab(ax, P_, s, dx, dy)
    dot(ax, G, s=2.4, hollow=True)
    lab(ax, E(118), "$J$", -6, 7)
    lab(ax, np.array([0.62, 0.22]), "$K$")
    geo(ax, (-1.2, 1.2), (-0.95, 0.95))
    save(fig, "fig4_crossing")


# ---------------------------------------------------------------- Figure 5
def refl(x, c, ang):
    c2, s2 = np.cos(np.deg2rad(2 * ang)), np.sin(np.deg2rad(2 * ang))
    d = np.asarray(x) - c
    return c + np.array([c2 * d[0] + s2 * d[1], s2 * d[0] - c2 * d[1]])


def fig5():
    d = 1.6
    U, V = O, np.array([d, 0.0])
    A, D = np.array([0.45, 0.0]), np.array([d - 0.45, 0.0])
    B, C = 0.6 * p(30), V + 0.6 * p(150)
    X, Y = np.array([0.22, 0.10]), np.array([d - 0.22, 0.10])
    chain = [A, X, B, C, Y, D]
    # X is the point of the chain farthest from V, and Y the one farthest from U
    assert max(np.linalg.norm(P_ - V) for P_ in chain) == np.linalg.norm(X - V) > 1
    assert max(np.linalg.norm(P_ - U) for P_ in chain) == np.linalg.norm(Y - U) > 1
    fig = plt.figure(figsize=(TW, 2.3))
    ax = fig.add_axes([0.0, 0.0, 0.5, 1.0])
    seg(ax, U - np.array([0.08, 0]), V + np.array([0.08, 0]))
    seg(ax, U, 1.15 * p(30), color=BLUE); seg(ax, V, V + 1.15 * p(150), color=ORANGE)
    seg(ax, 1.15 * p(30), 1.62 * p(30), lw=LT, ls=(0, (4, 1.5, 1, 1.5)), color=BLUE)
    seg(ax, V + 1.15 * p(150), V + 1.62 * p(150), lw=LT, ls=(0, (4, 1.5, 1, 1.5)), color=ORANGE)
    path(ax, arcpts(0, 30, 1.0, U), lw=LT, ls=DOTS, color=BLUE)       # unit circle about U
    path(ax, arcpts(150, 180, 1.0, V), lw=LT, ls=DOTS, color=ORANGE)  # unit circle about V
    path(ax, chain)
    for P_, s, dx, dy, col in [(U, "$U$", -5, -8, BLUE), (V, "$V$", 5, -8, ORANGE), (A, "$A$", 0, -8, INK),
                               (D, "$D$", 0, -8, INK), (B, "$B$", -4, 7, INK), (C, "$C$", 4, 7, INK),
                               (X, "$X$", -3, 8, ORANGE), (Y, "$Y$", 3, 8, BLUE)]:
        dot(ax, P_, fc=col); lab(ax, P_, s, dx, dy)
    lab(ax, 1.62 * p(30), r"$\ell_U$", 7, 1, size=8); lab(ax, V + 1.62 * p(150), r"$\ell_V$", -7, 1, size=8)
    geo(ax, (-0.15, 1.75), (-0.22, 1.00))
    fig.text(0.0, 0.985, "(a)", ha="left", va="top")

    W = d * p(60)
    Ap, Xp = refl(A, U, 30), refl(X, U, 30)
    Yp, Dp = refl(Y, V, 150), refl(D, V, 150)
    bx = fig.add_axes([0.5, 0.0, 0.5, 1.0])
    for (a, b), col in [((240, 270), ORANGE_T), ((270, 300), BLUE_T)]:
        bx.add_patch(Polygon([W, W + p(a), W + p(b)], closed=True, fc=col, ec="none"))
        seg(bx, W + p(a), W + p(b), lw=LT)
    seg(bx, W, W + 1.62 * p(240)); seg(bx, W, W + 1.62 * p(300)); seg(bx, W, W + p(270), lw=LT)
    path(bx, arcpts(240, 300, 1.0, W), lw=LT, ls=DOTS)
    path(bx, [Ap, Xp, B, C, Yp, Dp])
    seg(bx, Xp, Yp, lw=LW, ls=DASH)
    for P_, s, dx, dy, col in [(W, "$W$", 0, 7, INK), (Ap, "$A'$", -9, 3, INK), (Xp, "$X'$", -3, -9, ORANGE),
                               (B, "$B$", 0, -8, INK), (C, "$C$", 0, -8, INK), (Yp, "$Y'$", 3, -9, BLUE),
                               (Dp, "$D'$", 9, 3, INK)]:
        dot(bx, P_, fc=col); lab(bx, P_, s, dx, dy)
    bx.add_patch(Arc(W, 0.34, 0.34, theta1=240, theta2=300, lw=LT)); lab(bx, W + 0.27 * p(270), deg(60), 0, -3, size=7)
    geo(bx, (-0.22, 1.82), (-0.10, 1.50))
    fig.text(0.5, 0.985, "(b)", ha="left", va="top")
    save(fig, "fig5_case1")


# ---------------------------------------------------------------- Figure 6
def fig6():
    e0, e1 = p(90), p(-60)
    b = 2 * np.cos(np.deg2rad(15)) - 1
    Q = e0 + b * e1
    R = lambda th, v: np.array([np.cos(np.deg2rad(th)) * v[0] - np.sin(np.deg2rad(th)) * v[1],
                                np.sin(np.deg2rad(th)) * v[0] + np.cos(np.deg2rad(th)) * v[1]])
    th = np.linspace(150, 210, 1201)
    names = ["A,X_V", "A,B,X_V", "A,B,C,X_V", r"A,\ldots,D,X_V", r"A,\ldots,E,X_V"]
    curves = {k_: [] for k_ in names}
    for t in th:
        x = np.sin(np.deg2rad(t)); a = -(2 + x) / 3; c = -(1 + x) / 3
        curves[names[0]].append(a * e0 + R(t, Q))
        curves[names[1]].append(a * e0 - b * e1 + R(t, Q))
        curves[names[2]].append(a * e0 - b * e1 + R(t, Q + c * e0))
        curves[names[3]].append(-Q + R(t, Q + c * e0))
        curves[names[4]].append(-Q + (1 + c) * R(t, e0))
    norms = {k_: np.linalg.norm(np.array(v), axis=1) for k_, v in curves.items()}
    assert np.allclose(norms[names[4]], norms[names[0]][::-1]) and np.allclose(norms[names[3]], norms[names[1]][::-1])
    ramp = ["#084c34", "#14926a", "#4fbf96"]          # one hue, dark to light, in the order of the prefixes
    styles = ["-", (0, (4, 2)), (0, (1, 1.2))]
    fig = plt.figure(figsize=(3.5, 2.2))
    bx = fig.add_axes([0.15, 0.17, 0.82, 0.80])
    bx.axhspan(49 / 50, 1.0, fc=GREY2, ec="none")
    bx.axhline(49 / 50, color=MUTED, lw=LT)
    bx.axhline(1.0, color=INK, lw=LT)
    for k_, mirror, col in [(names[0], names[4], ramp[0]), (names[1], names[3], ramp[1])]:
        bx.plot(th, norms[mirror], color=col, lw=0.6)          # mirror images in theta = 180
    for k_, col, st in zip(names[:3], ramp, styles):
        bx.plot(th, norms[k_], color=col, lw=1.5, ls=st, label="$\\{" + k_ + "\\}$")
    bx.text(209.5, 0.99, r"$49/50$", fontsize=7, va="center", ha="right")
    bx.set_xlim(150, 210); bx.set_xticks([150, 165, 180, 195, 210])
    bx.set_xticklabels([deg(150), deg(165), deg(180), deg(195), deg(210)])
    bx.set_xlabel(r"$\vartheta$", labelpad=1)
    bx.set_ylabel(r"$|S|$", labelpad=2)
    bx.set_ylim(0.70, 1.012); bx.set_yticks([0.7, 0.8, 0.9, 1.0])
    for s in ("top", "right"):
        bx.spines[s].set_visible(False)
    bx.legend(frameon=False, fontsize=8, loc="lower left", handlelength=2.6, ncol=1)
    save(fig, "fig6_case2")


# ---------------------------------------------------------------- Figure 7
def load_rectangles():
    base = os.path.join(HERE, "..", "verification")
    files = {2: ("Case2/case2_actual_mixed_prefix_data.json", "boxes"),
             3: ("Case3/case3_mixed_certificate_data.json", "mixed_boxes"),
             4: ("Case4/case4_mixed_certificate_data.json", "mixed_boxes")}
    out = {}
    for case, (name, key) in files.items():
        with open(os.path.join(base, name), encoding="utf-8") as handle:
            rows = json.load(handle)[key]
        out[case] = [(tuple(float(Fraction(x)) for x in r["theta"]), tuple(float(Fraction(x)) for x in r["phi"]),
                      float(Fraction(r["lambda"])), r["a"] == "-sqrt(3)*lambda") for r in rows]
    return out


def fig7():
    fig = plt.figure(figsize=(TW, 3.75))
    # (a) a chain of subcase R_U C_V in Case 2, drawn with theta = 180 and phi = 11
    ax = fig.add_axes([0.05, 0.40, 0.90, 0.60])
    U, V = O, np.array([1.5, 0.4])
    th = 180.0
    r0, r1 = p(th), p(th + 30)
    A, B, C, D = np.array([0.62, 0.0]), 0.66 * p(30), np.array([0.78, 0.4]), np.array([0.74, 0.0])
    E, F = V + 0.52 * r1, np.array([1.0, 0.4])
    Xm, Xp = V + 0.997 * p(th + 9), V + 0.997 * p(th + 20)
    phi, r = 11.0, 1.12
    Y = r * p(phi)
    chain = [A, Xp, Xm, B, C, D, E, Y, F]
    # checks: the chain lies in both wedges and in the unit sector at V; Y is its point farthest
    # from U; the two witnesses lie in the caps of that sector and maximize the two cap functionals
    for P_ in chain:
        ang = np.degrees(np.arctan2(P_[1], P_[0])); w = P_ - V
        angv = np.degrees(np.arctan2(w[1], w[0])) % 360
        assert -1e-9 <= ang <= 30 + 1e-9 and th - 1e-9 <= angv <= th + 30 + 1e-9 and np.linalg.norm(w) <= 1
    assert max(np.linalg.norm(P_) for P_ in chain) == np.linalg.norm(Y) > 1
    assert max((P_ - V) @ p(th + 7.5) for P_ in chain) == (Xm - V) @ p(th + 7.5) > RHO
    assert max((P_ - V) @ p(th + 22.5) for P_ in chain) == (Xp - V) @ p(th + 22.5) > RHO
    hull = [A, D, E, Y, F, C, B, Xm, Xp]
    ax.add_patch(Polygon(hull, closed=True, fc=GREY2, ec="none"))
    # unit sector at U (blue) and at V (orange, dashed rays)
    seg(ax, U, np.array([1.66, 0]), color=BLUE); seg(ax, U, 1.0 * p(30), color=BLUE)
    path(ax, arcpts(0, 30), lw=LW, color=BLUE)
    seg(ax, V, V + r0, color=ORANGE, ls=DASH); seg(ax, V, V + r1, color=ORANGE, ls=DASH)
    path(ax, arcpts(th, th + 30, 1.0, V), lw=LT, color=ORANGE)
    # the two caps of the unit sector at V, with their chords
    for t0 in (th, th + 15):
        ax.add_patch(Polygon(arcpts(t0, t0 + 15, 1.0, V), closed=True, fc=GREEN, ec=GREEN, lw=0.3))
        seg(ax, V + p(t0), V + p(t0 + 15), color=GREEN, lw=LW)
    seg(ax, U, Y, lw=LT, ls=DOTS, color=BLUE)
    ax.add_patch(Arc(U, 0.44, 0.44, theta1=0, theta2=phi, lw=LT, color=BLUE))
    lab(ax, 0.275 * p(phi / 2), r"$\phi$", size=8)
    lab(ax, 0.47 * p(phi), r"$r$", 0, 5, size=8)
    path(ax, chain)
    box = dict(fc="white", ec="none", pad=0.4)
    marks = [(A, "$A$", 4, -8, INK, None), (B, "$B$", 8, -5, INK, None), (C, "$C$", 7, -7, INK, None),
             (D, "$D$", 3, -8, INK, None), (E, "$E$", 7, -5, INK, None), (F, "$F$", 0, 7, INK, None),
             (Xm, "$X_V^-$", -13, 4, GREEN, box), (Xp, "$X_V^+$", -14, -3, GREEN, box),
             (Y, "$Y$", 8, 1, BLUE, None), (U, "$U$", -6, -6, BLUE, None), (V, "$V$", 7, 3, ORANGE, None)]
    for P_, s, dx, dy, col, bb in marks:
        dot(ax, P_, s=3.2, fc=col)
        lab(ax, P_, s, dx, dy, **({} if bb is None else {"bbox": bb}))
    lab(ax, np.array([1.66, 0]), "$U_0$", -2, -7, size=7); lab(ax, 0.93 * p(30), "$U_1$", -8, 3, size=7)
    lab(ax, V + 0.26 * r0, "$V_0$", 0, 6, size=7); lab(ax, V + 0.22 * r1, "$V_1$", 6, -7, size=7)
    lab(ax, 0.36 * p(22), r"$\Pi_U$", size=8); lab(ax, V + 0.36 * p(th + 13), r"$\Pi_V$", size=8)
    geo(ax, (-0.10, 1.72), (-0.16, 0.56))
    fig.text(0.005, 0.99, "(a)", ha="left", va="top")

    # (b) the rectangles of the three certificates, colored by the weight lambda
    cmap = LinearSegmentedColormap.from_list("weight", [(0.0, "#b8481c"), (0.25, ORANGE), (0.5, "#e9e8e4"),
                                                        (0.75, BLUE), (1.0, "#184f95")])
    norm = Normalize(0.30, 0.70)
    data = load_rectangles()
    for i, case in enumerate((2, 3, 4)):
        bx = fig.add_axes([0.10 + i * 0.275, 0.085, 0.225, 0.215])
        for (t0, t1), (f0, f1), lam, special in data[case]:
            assert norm.vmin <= lam <= norm.vmax
            bx.add_patch(Rectangle((t0, f0), t1 - t0, f1 - f0, fc=cmap(norm(lam)), ec="white", lw=0.8))
            if special:        # Lemma 5.5(i) settles one prefix here
                bx.add_patch(Rectangle((t0, f0), t1 - t0, f1 - f0, fc="none", ec=INK, lw=0, hatch="//////"))
        if case != 2:          # the prefix of Lemma 5.5(ii) has norm one on this line
            bx.plot([165, 195], [0, 30], color=INK, lw=LW, ls=DASH)
        bx.set_xlim(150, 210); bx.set_ylim(0, 30)
        bx.set_xticks([150, 180, 210]); bx.set_xticklabels([deg(150), deg(180), deg(210)], fontsize=7)
        bx.set_yticks([0, 15, 30])
        bx.set_yticklabels([deg(0), deg(15), deg(30)] if i == 0 else [], fontsize=7)
        bx.set_xlabel(r"$\vartheta$", labelpad=0)
        if i == 0:
            bx.set_ylabel(r"$\phi$", labelpad=1)
        bx.set_title("Case %d" % case, fontsize=9, pad=3)
        bx.tick_params(length=2)
    cax = fig.add_axes([0.935, 0.085, 0.014, 0.215])
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, ticks=[0.3, 0.4, 0.5, 0.6, 0.7])
    cb.ax.tick_params(labelsize=7, length=2, width=0.5)
    cb.outline.set_linewidth(0.5)
    cb.ax.set_title(r"$\lambda$", fontsize=9, pad=3)
    fig.text(0.005, 0.375, "(b)", ha="left", va="top")
    save(fig, "fig7_mixed")


# ---------------------------------------------------------------- Figure 8
def fig8():
    s111 = np.sqrt(111.0)
    Q0 = np.array([3 / (4 * s111), 0]); Q1 = np.array([10 / s111, 0])
    Q2 = np.array([5 / (4 * s111), 5 * np.sqrt(3) / (12 * s111)]); Q3 = np.array([21 / (2 * s111), 5 * np.sqrt(3) / (12 * s111)])
    P = np.array([0.75, np.sqrt(3) / 4])
    fig = plt.figure(figsize=(4.4, 2.5))
    ax = fig.add_axes([0, 0, 1, 1])
    H = np.vstack([O, arcpts(0, 15), p(30)])
    ax.add_patch(Polygon(H, closed=True, fc=GREEN_T, ec=INK, lw=LT))                 # H_15, of which
    ax.add_patch(Polygon([O, p(0), Q3, P], closed=True, fc=GREY, ec=INK, lw=LT))     # K_0 is kept
    d_ = 0.7
    St = [np.array([1 - d_, 0]), 2 * (1 - d_) * P, 2 * (1 - d_) * P + (2 * d_ - 1) * p(0), p(0)]
    path(ax, St, lw=LW, ls=DASH)
    path(ax, [Q0, Q2, Q1, Q3], lw=LC)
    for P_, s, dx, dy in [(O, "$O$", -7, -5), (p(0), "$p_0$", 7, -5), (p(15), "$p_{15}$", 9, 0), (p(30), "$p_{30}$", 3, 6),
                          (P, "$P_*$", -7, 5), (Q0, "$Q_0$", 0, -8), (Q1, "$Q_1$", -3, -8), (Q2, "$Q_2$", -3, 7),
                          (Q3, "$Q_3$", 9, 2)]:
        dot(ax, P_); lab(ax, P_, s, dx, dy)
    lab(ax, np.array([0.55, 0.10]), "$K_0$")
    geo(ax, (-0.08, 1.10), (-0.09, 0.57))
    save(fig, "fig8_example")


if __name__ == "__main__":
    for f in (fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8):
        f()
    print("done")
