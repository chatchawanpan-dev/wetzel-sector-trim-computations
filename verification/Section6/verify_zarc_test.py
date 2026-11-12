#!/usr/bin/env python3
"""Exact check of the finite inequalities in Section 6.2 of the paper.

The z-arc Z of Example 6.4 has end segments of length 1/12 and a middle
segment of length 5/6; its convex hull is a parallelogram with sides 1/12
and sqrt(111)/12 that make the angle beta, where cos(beta) = 6/sqrt(111)
and sin(beta) = 5*sqrt(3)/sqrt(111).

For an angle chi put

    f_chi(phi) = sin(120 - chi)|sin(phi)| + cos(chi)|cos(phi - 120)|
                 + (1/2)|cos(phi - chi)|        (angles in degrees),

    F(chi, eps, theta) = f_chi(theta)/12 + sqrt(111) f_chi(theta + eps*beta)/12,

with eps = +1 for direct copies of Z and eps = -1 for reflected copies
(the upper and the lower sign in the paper).

Lemma 6.5 of the paper shows that a copy of Z in the triangle
x_2 >= 0, <x, p_120> <= 0, <x, p_chi> <= h satisfies
F(chi, eps, theta) <= h, and that F(chi, eps, .) is concave between
consecutive points of

    {0, 30, 90 + chi} and {0, 30, 90 + chi} - eps*beta   (mod 180).

This program checks, with rational interval arithmetic:

  (1) Proposition 6.6: F(16.8, eps, theta) > cos(13.2) at these twelve
      points; so Z has no congruent copy in H_t for t <= 3.6 degrees.
  (2) Proposition 6.6: the point Q_3 of Example 6.4 has norm less than one
      and satisfies <Q_3, p_16.85> <= cos(13.15); so the copy of Z in
      Example 6.4 lies in H_t for t >= 3.7 degrees.
  (3) Remark 6.8: F(15, +1, theta) > cos(5) at the breakpoints in [0, 140]
      and at 140, and F(25, +1, theta) > cos(5) at the breakpoints in
      [140, 180] and at 140 and 180; so Z has no congruent copy in the
      pentagon conv{O, p_0, p_10, p_20, p_30}.

Only the Python standard library is used.  Run without -O.
"""
from __future__ import annotations

import datetime
import platform
import sys
from fractions import Fraction as F

if sys.flags.optimize:
    raise SystemExit("ERROR: run without -O/-OO or PYTHONOPTIMIZE; the checks are assertions.")


class I:
    """Closed interval with rational end points."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = F(lo)
        self.hi = self.lo if hi is None else F(hi)
        if self.lo > self.hi:
            raise ValueError((self.lo, self.hi))

    def __add__(self, other):
        other = as_interval(other)
        return I(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-as_interval(other))

    def __rsub__(self, other):
        return as_interval(other) - self

    def __mul__(self, other):
        other = as_interval(other)
        products = (self.lo * other.lo, self.lo * other.hi, self.hi * other.lo, self.hi * other.hi)
        return I(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = as_interval(other)
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError("interval contains zero")
        return self * I(1 / other.hi, 1 / other.lo)

    def abs(self):
        if self.lo >= 0:
            return self
        if self.hi <= 0:
            return -self
        return I(0, max(-self.lo, self.hi))

    def __repr__(self):
        return f"[{float(self.lo):.10f}, {float(self.hi):.10f}]"


def as_interval(x):
    return x if isinstance(x, I) else I(x)


def outward(x: I, digits: int = 40) -> I:
    scale = 10 ** digits
    lo = F((x.lo * scale).numerator // (x.lo * scale).denominator, scale)
    hi = F(-((-(x.hi * scale).numerator) // (x.hi * scale).denominator), scale)
    return I(lo, hi)


def sqrt_interval(n: int, digits: int = 40) -> I:
    """Rational enclosure of sqrt(n), from the integer square root."""
    from math import isqrt

    scale = 10 ** digits
    root = isqrt(n * scale * scale)
    lo, hi = F(root, scale), F(root + 1, scale)
    assert lo * lo <= n <= hi * hi
    return I(lo, hi)


def atan_inverse(q: int, terms: int) -> I:
    """Enclosure of arctan(1/q) by the alternating series."""
    x = F(1, q)
    total, power = F(0), x
    for k in range(terms):
        term = power / (2 * k + 1)
        total += term if k % 2 == 0 else -term
        power *= x * x
    nxt = power / (2 * terms + 1)
    return I(total - nxt, total) if terms % 2 == 1 else I(total, total + nxt)


PI = outward(16 * atan_inverse(5, 40) - 4 * atan_inverse(239, 14), 50)  # Machin


def _sin_series(x: F, pairs: int = 14) -> I:
    lo = up = F(0)
    term, x2 = x, x * x
    for k in range(2 * pairs + 2):
        if k:
            term = term * x2 / ((2 * k) * (2 * k + 1))
        signed = term if k % 2 == 0 else -term
        if k <= 2 * pairs:
            up += signed
        lo += signed
    return I(lo, up)


def _cos_series(x: F, pairs: int = 14) -> I:
    lo = up = F(0)
    term, x2 = F(1), x * x
    for k in range(2 * pairs + 2):
        if k:
            term = term * x2 / ((2 * k - 1) * (2 * k))
        signed = term if k % 2 == 0 else -term
        if k <= 2 * pairs:
            up += signed
        lo += signed
    return I(lo, up)


_CACHE: dict[F, tuple[I, I]] = {}


def sin_cos(deg) -> tuple[I, I]:
    """Enclosures of sin and cos of a rational number of degrees."""
    a = F(deg) % 360
    if a in _CACHE:
        return _CACHE[a]
    quadrant = int(a // 90)
    r = a - 90 * quadrant
    swap = r > 45
    base = 90 - r if swap else r
    # 0 <= base <= 45 degrees, so the alternating series bounds apply.
    x = I(base * PI.lo / 180, base * PI.hi / 180)
    s = outward(I(_sin_series(x.lo).lo, _sin_series(x.hi).hi))
    c = outward(I(_cos_series(x.hi).lo, _cos_series(x.lo).hi))
    if swap:
        s, c = c, s
    result = {0: (s, c), 1: (c, -s), 2: (-s, -c), 3: (-c, s)}[quadrant]
    _CACHE[a] = result
    return result


# Proposition 6.6, first statement: t = 3.6, so m_t = 16.8 and delta_t = 13.2.
CHI_SMALL_TRIM = F(84, 5)
DELTA_SMALL_TRIM = F(66, 5)

SQRT3 = sqrt_interval(3)
SQRT111 = sqrt_interval(111)
SIN_BETA = 5 * SQRT3 / SQRT111
COS_BETA = I(6) / SQRT111
# sin^2 + cos^2 = 75/111 + 36/111 = 1 exactly.


class Angle:
    """The angle  deg + k*beta  (deg rational, k integer)."""

    def __init__(self, deg, k: int = 0):
        self.deg, self.k = F(deg), int(k)

    def shift(self, deg=0, k=0):
        return Angle(self.deg + F(deg), self.k + k)

    def sin_cos(self) -> tuple[I, I]:
        s, c = sin_cos(self.deg)
        sb, cb = I(0), I(1)
        for _ in range(abs(self.k)):  # angle addition for k*beta
            sb, cb = sb * COS_BETA + cb * SIN_BETA, cb * COS_BETA - sb * SIN_BETA
        if self.k < 0:
            sb = -sb
        return s * cb + c * sb, c * cb - s * sb

    def describe(self) -> str:
        if self.k == 0:
            return f"{float(self.deg):g}"
        if self.deg == 0 and self.k == 1:
            return "beta"
        beta = " + beta" if self.k == 1 else (" - beta" if self.k == -1 else f" + {self.k} beta")
        return f"{float(self.deg):g}{beta}"


def f(chi: F, phi: Angle) -> I:
    """Enclosure of f_chi(phi)."""
    s120, _ = sin_cos(120 - chi)
    _, cchi = sin_cos(chi)
    s, _ = phi.sin_cos()
    _, c120 = phi.shift(deg=-120).sin_cos()
    _, cchi_shift = phi.shift(deg=-chi).sin_cos()
    return s120 * s.abs() + cchi * c120.abs() + F(1, 2) * cchi_shift.abs()


def big_f(chi: F, eps: int, theta: Angle) -> I:
    return (f(chi, theta) + SQRT111 * f(chi, theta.shift(k=eps))) / 12


def breakpoints(chi: F, eps: int) -> list[Angle]:
    base = [F(0), F(30), 90 + chi]
    pts = [Angle(d) for d in base]
    for d in base:
        # d - eps*beta, moved into [0, 180) by a multiple of 180 degrees
        # (55 < beta < 56, checked below).
        value_lo = d - eps * F(56) if eps > 0 else d + F(55)
        shift = 0
        while value_lo + shift < 0:
            shift += 180
        while value_lo + shift >= 180:
            shift -= 180
        pts.append(Angle(d + shift, -eps))
    return pts


def check_beta_range():
    # 55 degrees < beta < 56 degrees, by comparing tangents (increasing on (0, 90)).
    tan_beta = SIN_BETA / COS_BETA
    s55, c55 = sin_cos(55)
    s56, c56 = sin_cos(56)
    assert (s55 / c55).hi < tan_beta.lo and tan_beta.hi < (s56 / c56).lo


def angle_enclosure(a: Angle) -> tuple[F, F]:
    """Crude enclosure in degrees, using 55 < beta < 56."""
    lo = a.deg + (55 * a.k if a.k >= 0 else 56 * a.k)
    hi = a.deg + (56 * a.k if a.k >= 0 else 55 * a.k)
    return lo, hi


def check_family(name: str, chi, eps: int, bound_deg, window=None, extra=()):
    """Assert F(chi, eps, theta) > cos(bound_deg) at the breakpoints in the window."""
    chi = F(chi)
    _, cos_bound = sin_cos(F(bound_deg))
    points = []
    for a in breakpoints(chi, eps):
        lo, hi = angle_enclosure(a)
        if window is None:
            points.append(a)
        elif window[0] <= lo and hi <= window[1]:
            points.append(a)
        else:
            # the breakpoint must lie strictly outside the window
            assert hi < window[0] or lo > window[1], (name, a.describe())
    points.extend(Angle(d) for d in extra)
    worst = None
    for a in points:
        value = big_f(chi, eps, a)
        assert value.lo > cos_bound.hi, (name, a.describe(), value, cos_bound)
        margin = value.lo - cos_bound.hi
        if worst is None or margin < worst[0]:
            worst = (margin, a)
        print(f"  {name}  eps={eps:+d}  theta = {a.describe():<16} F in {value}  > {float(cos_bound.hi):.10f}")
    return len(points), worst


def main() -> int:
    print("Section 6.2 checks: the arc Z in small trims and in the pentagon")
    print("UTC:", datetime.datetime.now(datetime.timezone.utc).isoformat())
    print("Python:", sys.version.replace("\n", " "))
    print("Platform:", platform.platform())
    check_beta_range()

    print("(1) Proposition 6.6: chi = 16.8, h = cos(13.2)  [H_t for t = 3.6]")
    total = 0
    worst = None
    for eps in (+1, -1):
        n, w = check_family("H_3.6", CHI_SMALL_TRIM, eps, DELTA_SMALL_TRIM)
        total += n
        worst = w if worst is None or w[0] < worst[0] else worst
    assert total == 12
    print(f"  smallest margin {float(worst[0]):.3e} at theta = {worst[1].describe()}")

    print("(2) Proposition 6.6: the copy of Example 6.4 lies in H_t for t = 3.7")
    x3 = I(21) / (2 * SQRT111)
    y3 = 5 * SQRT3 / (12 * SQRT111)
    assert F(5317, 5328) < 1  # |Q_3|^2 = 5317/5328
    norm2 = x3 * x3 + y3 * y3
    assert norm2.lo <= F(5317, 5328) <= norm2.hi
    assert norm2.hi < 1
    s, c = sin_cos(F(337, 20))  # m_t = 16.85 for t = 3.7
    _, c_delta = sin_cos(F(263, 20))  # delta_t = 13.15
    value = x3 * c + y3 * s
    assert value.hi < c_delta.lo, (value, c_delta)
    print(f"  <Q_3, p_16.85> in {value}  <  cos(13.15), which is at least {float(c_delta.lo):.10f}")

    print("(3) Remark 6.8: the pentagon conv{O, p_0, p_10, p_20, p_30}")
    n1, w1 = check_family("K10", 15, +1, 5, window=(F(0), F(140)), extra=(F(140),))
    n2, w2 = check_family("K10", 25, +1, 5, window=(F(140), F(180)), extra=(F(140), F(180)))
    assert (n1, n2) == (6, 3), (n1, n2)
    print(f"  smallest margins {float(w1[0]):.3e} (chi = 15) and {float(w2[0]):.3e} (chi = 25)")

    print("SECTION 6 CHECKS: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
