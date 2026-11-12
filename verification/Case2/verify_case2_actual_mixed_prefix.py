#!/usr/bin/env python3
"""Exact-rational verifier for actual Case 2 mixed-state prefix certificates.

Reconstructs the 12 actual mixed-state words, checks exact coverage of
[150,210] degrees x [0,30] degrees by stored coefficient rectangles, and
verifies their proper-prefix norms with outward rational intervals.
The special boundary prefix uses the analytic identity stated in the proof.
No optimizer or third-party package is used by this verifier.
"""
from __future__ import annotations
from fractions import Fraction as F
import itertools, math, sys, time

if sys.flags.optimize:
    raise SystemExit("ERROR: Python optimization disables required proof checks; run without -O/-OO or PYTHONOPTIMIZE.")

class I:
    __slots__=('lo','hi')
    def __init__(self,lo,hi=None):
        self.lo=lo if isinstance(lo,F) else F(lo)
        self.hi=self.lo if hi is None else (hi if isinstance(hi,F) else F(hi))
        if self.lo>self.hi: raise ValueError((self.lo,self.hi))
    def __add__(self,o):o=asI(o);return I(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,o):return self+(-asI(o))
    def __rsub__(self,o):return asI(o)-self
    def __mul__(self,o):
        o=asI(o);v=(self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi);return I(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=asI(o)
        if o.lo<=0<=o.hi:raise ZeroDivisionError
        return self*I(1/o.hi,1/o.lo)
    def __repr__(self):return f'[{float(self.lo):.12g},{float(self.hi):.12g}]'
def asI(x):return x if isinstance(x,I) else I(x)
def Q(x):return F(str(x))
def floor_fraction(x):return x.numerator//x.denominator
def ceil_fraction(x):return -((-x.numerator)//x.denominator)
def outward_decimal(x,d=28):
    s=10**d;return I(F(floor_fraction(x.lo*s),s),F(ceil_fraction(x.hi*s),s))
def atan_iv(q,n):
    x=F(1,q);tot=F(0);power=x
    for k in range(n):
        term=power/(2*k+1);tot += term if k%2==0 else -term;power*=x*x
    nxt=power/(2*n+1)
    return I(tot-nxt,tot) if n%2==1 else I(tot,tot+nxt)
def pi_iv():return outward_decimal(16*atan_iv(5,32)-4*atan_iv(239,10),55)
PI=pi_iv()
def sin_point(x,pairs=12):
    up=F(0);lo=F(0);x2=x*x;term=x
    for k in range(2*pairs+2):
        if k:term=term*x2/((2*k)*(2*k+1))
        signed=term if k%2==0 else -term
        if k<=2*pairs:up+=signed
        lo+=signed
    return I(lo,up)
def cos_point(x,pairs=12):
    up=F(0);lo=F(0);x2=x*x;term=F(1)
    for k in range(2*pairs+2):
        if k:term=term*x2/((2*k-1)*(2*k))
        signed=term if k%2==0 else -term
        if k<=2*pairs:up+=signed
        lo+=signed
    return I(lo,up)
PC={}
def sin_cos_deg(a):
    a%=360
    if a in PC:return PC[a]
    q=floor_fraction(a/90);r=a-90*q;swap=r>45;base=90-r if swap else r
    x=I(base*PI.lo/180,base*PI.hi/180)
    sb=I(sin_point(x.lo).lo,sin_point(x.hi).hi);cb=I(cos_point(x.hi).lo,cos_point(x.lo).hi)
    sb=outward_decimal(sb);cb=outward_decimal(cb);s,c=(cb,sb) if swap else (sb,cb)
    if q==0:z=(s,c)
    elif q==1:z=(c,-s)
    elif q==2:z=(-s,-c)
    elif q==3:z=(-c,s)
    else:raise AssertionError
    PC[a]=z;return z

def has_lattice(lo,hi,start,period=360):
    # exists start+k period in [lo,hi]
    k=ceil_fraction((lo-start)/period)
    return start+period*k<=hi
IC={}
def trig_interval(lo,hi):
    key=(lo,hi)
    if key in IC:return IC[key]
    if hi-lo>=360:return I(-1,1),I(-1,1)
    slo,clo=sin_cos_deg(lo);shi,chi=sin_cos_deg(hi)
    s=I(min(slo.lo,shi.lo),max(slo.hi,shi.hi));c=I(min(clo.lo,chi.lo),max(clo.hi,chi.hi))
    if has_lattice(lo,hi,F(90)):s=I(s.lo,1)
    if has_lattice(lo,hi,F(270)):s=I(-1,s.hi)
    if has_lattice(lo,hi,F(0)):c=I(c.lo,1)
    if has_lattice(lo,hi,F(180)):c=I(-1,c.hi)
    z=(outward_decimal(s),outward_decimal(c));IC[key]=z;return z

def ppoint(a):s,c=sin_cos_deg(a);return [c,s]
def pinterval(lo,hi):s,c=trig_interval(lo,hi);return [c,s]
def vadd(a,b):return [a[0]+b[0],a[1]+b[1]]
def vscale(s,v):return [s*v[0],s*v[1]]
def sq_upper(x):
    if x.lo<=0<=x.hi:return max(x.lo*x.lo,x.hi*x.hi)
    return max(x.lo*x.lo,x.hi*x.hi)
def norm2up(v):return sq_upper(v[0])+sq_upper(v[1])
ZERO=[I(0),I(0)]
P90=ppoint(F(90));PM60=ppoint(F(-60));SQRT3=2*sin_cos_deg(F(60))[0]
RHO=sin_cos_deg(F(15,2))[1]

SUP = tuple("ABCDEF")
WORDS = []
for pre in itertools.permutations(("A", "XVm", "XVp")):
    for tail in [("E", "XR", "F"), ("E", "F", "XR")]:
        WORDS.append(tuple(list(pre) + ["B", "C", "D"] + list(tail)))
PREFIXES = {frozenset(w[:k]) for w in WORDS for k in range(1, len(w))}
EXCEPTION = frozenset(("A", "B", "XVm", "XVp"))


def coeff_mix(labels, theta, phi, lam, a_spec, c):
    lam = I(lam)
    a = -SQRT3*lam if a_spec == "special" else I(a_spec)
    c = I(c)
    sp, cp = trig_interval(phi[0], phi[1])
    alpha = -lam*(sp + SQRT3*cp)
    beta = -2*lam*cp
    mu = (I(1)-lam)/(2*RHO)
    gamma = -2*(I(1)-lam)*sin_cos_deg(F(15))[1]
    out = [I(0), I(0)]
    vectors = {
        "A": vscale(a, P90),
        "B": vscale(beta, PM60),
        "D": vscale(alpha-a, P90),
        "XR": vscale(lam, pinterval(*phi)),
        "C": vscale(c, pinterval(theta[0]+90, theta[1]+90)),
        "E": vscale(gamma, pinterval(theta[0]-60, theta[1]-60)),
        "F": vscale(gamma-c, pinterval(theta[0]+90, theta[1]+90)),
        "XVm": vscale(mu, pinterval(theta[0]+F(15,2), theta[1]+F(15,2))),
        "XVp": vscale(mu, pinterval(theta[0]+F(45,2), theta[1]+F(45,2))),
    }
    for label in labels:
        out = vadd(out, vectors[label])
    return out


stats = {"nodes": 0, "depth": 0, "worst_leaf_squared": F(0)}

def certify(labels, theta, phi, lam, a, c, depth=0):
    stats["nodes"] += 1
    stats["depth"] = max(stats["depth"], depth)
    value = norm2up(coeff_mix(labels, theta, phi, lam, a, c))
    if value <= 1:
        stats["worst_leaf_squared"] = max(stats["worst_leaf_squared"], value)
        return
    if depth >= 28:
        raise AssertionError(("uncertified prefix", sorted(labels), theta, phi, float(value)))
    if theta[1]-theta[0] >= phi[1]-phi[0]:
        mid = sum(theta)/2
        certify(labels, (theta[0], mid), phi, lam, a, c, depth+1)
        certify(labels, (mid, theta[1]), phi, lam, a, c, depth+1)
    else:
        mid = sum(phi)/2
        certify(labels, theta, (phi[0], mid), lam, a, c, depth+1)
        certify(labels, theta, (mid, phi[1]), lam, a, c, depth+1)


def parse_data(path):
    import json
    data = json.loads(open(path, encoding="utf-8").read())
    assert data["schema"] == "case2-actual-mixed-prefix-v1"
    assert list(map(Q, data["theta_domain"])) == [F(150), F(210)]
    assert list(map(Q, data["phi_domain"])) == [F(0), F(30)]
    rows, names = [], set()
    for row in data["boxes"]:
        name = row["name"]
        assert name not in names
        names.add(name)
        theta = tuple(map(Q, row["theta"]))
        phi = tuple(map(Q, row["phi"]))
        lam, c = Q(row["lambda"]), Q(row["c"])
        a = "special" if row["a"] == "-sqrt(3)*lambda" else Q(row["a"])
        assert F(150) <= theta[0] < theta[1] <= F(210)
        assert F(0) <= phi[0] < phi[1] <= F(30)
        assert F(0) < lam < F(1)
        rows.append((name, theta, phi, lam, a, c))
    tx = sorted({F(150), F(210)} | {x for _,t,_,_,_,_ in rows for x in t})
    px = sorted({F(0), F(30)} | {x for _,_,p,_,_,_ in rows for x in p})
    for t0,t1 in zip(tx,tx[1:]):
        for p0,p1 in zip(px,px[1:]):
            assert any(t[0] <= t0 and t1 <= t[1] and p[0] <= p0 and p1 <= p[1]
                       for _,t,p,_,_,_ in rows), ("coverage gap", t0,t1,p0,p1)
    return data, rows


def main(path):
    data, rows = parse_data(path)
    assert len(WORDS) == 12 and len(set(WORDS)) == 12 and len(PREFIXES) == 13
    checked = exceptions = 0
    for name, theta, phi, lam, a, c in rows:
        for labels in sorted(PREFIXES, key=lambda x: (len(x), sorted(x))):
            if a == "special" and labels == EXCEPTION:
                # Here S = W+(1-lambda)p_(theta+15), with
                # W=lambda(-cos(phi),sqrt(3)(cos(phi)-1)).
                # ||W||^2-lambda^2 =
                # 2 lambda^2 (2 cos(phi)-1)(cos(phi)-1) <= 0 on[0,60].
                # Thus ||S|| <= lambda+(1-lambda) = 1 exactly.
                assert F(0) <= phi[0] < phi[1] <= F(60)
                exceptions += 1
                continue
            certify(labels, theta, phi, lam, a, c)
            checked += 1
    print("schema =", data["schema"])
    print("actual mixed words =", len(WORDS))
    print("unique proper prefixes =", len(PREFIXES))
    print("coefficient boxes =", len(rows))
    print("interval prefix checks =", checked)
    print("analytic prefix checks =", exceptions)
    print("adaptive interval nodes =", stats["nodes"])
    print("maximum recursion depth =", stats["depth"])
    print("worst accepted leaf squared norm <=", float(stats["worst_leaf_squared"]))
    print("PASS: all actual Case-2 mixed-state prefix norms are at most 1")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("data", nargs="?", default="case2_actual_mixed_prefix_data.json")
    main(parser.parse_args().data)
