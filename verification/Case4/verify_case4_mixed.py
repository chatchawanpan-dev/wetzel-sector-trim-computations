#!/usr/bin/env python3
"""Independent exact-rational verifier for the H15 Case 4 mixed-state proof.

The script reconstructs the mixed-state branch grammar,
loads only finite-decimal coefficient data, and proves every required
proper-prefix norm bound with outward rational interval arithmetic.
It uses no optimizer and no non-standard Python package.
"""
from __future__ import annotations
from fractions import Fraction as F
import itertools, math, sys, time

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

MIX_PREFIXES=set();MIX_WORDS=[]
for pre in itertools.permutations(('A','XVm','XVp')):
 for tail in [('E','XR','F'),('E','F','XR')]:
  w=tuple(list(pre)+['B','C','D']+list(tail));MIX_WORDS.append(w)
  s=[]
  for x in w[:-1]:s.append(x);MIX_PREFIXES.add(frozenset(s))

def parse_coefficient_data(path):
    import json
    with open(path, encoding="utf-8") as handle:
        payload = json.load(handle)
    if payload.get("schema") != "case4-mixed-prefix-certificates-v1":
        raise AssertionError("unexpected certificate-data schema")
    if list(map(Q, payload["orientation_domain_degrees"])) != [F(150), F(210)]:
        raise AssertionError("unexpected orientation domain")

    if payload.get("mixed_state") != "R_U C_V":
        raise AssertionError("unexpected mixed state")
    if payload.get("mixed_weight_rule") != {
        "radial_weight": "lambda",
        "each_cap_weight": "(1-lambda)/(2*cos(7.5 degrees))",
    }:
        raise AssertionError("unexpected mixed weight rule")

    boxes = []
    seen = set()
    for row in payload["mixed_boxes"]:
        name = row["name"]
        if name in seen:
            raise AssertionError(("duplicate box", name))
        seen.add(name)
        tlo, thi = map(Q, row["theta"])
        plo, phi = map(Q, row["phi"])
        lam = Q(row["lambda"])
        if not (F(0) < lam < F(1)):
            raise AssertionError((name, "lambda outside (0,1)"))
        a_text = row["a"]
        a = "special" if a_text == "-sqrt(3)*lambda" else Q(a_text)
        c = Q(row["c"])
        if not (F(150) <= tlo < thi <= F(210)):
            raise AssertionError((name, "invalid theta interval"))
        if not (F(0) <= plo < phi <= F(30)):
            raise AssertionError((name, "invalid phi interval"))
        boxes.append((name, (tlo, thi), (plo, phi), lam, a, c))

    if len(boxes) != 21:
        raise AssertionError(("unexpected mixed-box count", len(boxes)))

    # Exact rectangular coverage check.  Every open cell and every grid
    # vertex determined by all endpoints must be covered by a stored closed
    # rectangle.
    tx = sorted({F(150), F(210)} | {x for _, t, _, _, _, _ in boxes for x in t})
    px = sorted({F(0), F(30)} | {x for _, _, p, _, _, _ in boxes for x in p})
    for i in range(len(tx) - 1):
        for j in range(len(px) - 1):
            tm = (tx[i] + tx[i + 1]) / 2
            pm = (px[j] + px[j + 1]) / 2
            if not any(t[0] <= tm <= t[1] and p[0] <= pm <= p[1]
                       for _, t, p, _, _, _ in boxes):
                raise AssertionError(("mixed-box coverage gap", tx[i], tx[i+1], px[j], px[j+1]))
    for t0 in tx:
        for p0 in px:
            if F(150) <= t0 <= F(210) and F(0) <= p0 <= F(30):
                if not any(t[0] <= t0 <= t[1] and p[0] <= p0 <= p[1]
                           for _, t, p, _, _, _ in boxes):
                    raise AssertionError(("mixed-box boundary gap", t0, p0))
    return payload, boxes


def coeff_mix(labels,theta,phi,lam,a_spec,c):
    lam=I(lam);one=I(1)
    a=-SQRT3*lam if a_spec=='special' else I(a_spec)
    c=I(c)
    sphi,cphi=trig_interval(phi[0],phi[1])
    alpha=-lam*(sphi+SQRT3*cphi);beta=-2*lam*cphi
    mu=(one-lam)/(2*RHO)
    # gamma=delta=-2(1-lam) cos15
    c15=sin_cos_deg(F(15))[1]
    gamma=delta=-2*(one-lam)*c15
    out=[I(0),I(0)]
    if 'A'in labels:out=vadd(out,vscale(a,P90))
    if 'B'in labels:out=vadd(out,vscale(beta,PM60))
    if 'C'in labels:out=vadd(out,vscale(alpha-a,P90))
    if 'XR'in labels:out=vadd(out,vscale(lam,pinterval(phi[0],phi[1])))
    if 'D'in labels:out=vadd(out,vscale(c,pinterval(theta[0]+90,theta[1]+90)))
    if 'E'in labels:out=vadd(out,vscale(delta,pinterval(theta[0]-60,theta[1]-60)))
    if 'F'in labels:out=vadd(out,vscale(gamma-c,pinterval(theta[0]+90,theta[1]+90)))
    if 'XVm'in labels:out=vadd(out,vscale(mu,pinterval(theta[0]+F(15,2),theta[1]+F(15,2))))
    if 'XVp'in labels:out=vadd(out,vscale(mu,pinterval(theta[0]+F(45,2),theta[1]+F(45,2))))
    return out

stats={'nodes':0,'maxdepth':0}
def certify_mix(labels,tlo,thi,plo,phi,name,lam,a,c,depth=0):
    stats['nodes']+=1;stats['maxdepth']=max(stats['maxdepth'],depth)
    v=coeff_mix(labels,(tlo,thi),(plo,phi),lam,a,c)
    if norm2up(v)<=1:return
    if depth>=30:raise AssertionError(('MIX',name,sorted(labels),float(tlo),float(thi),float(plo),float(phi),v,float(norm2up(v))))
    # split larger angular side
    if (thi-tlo) >= (phi-plo):
        m=(tlo+thi)/2;certify_mix(labels,tlo,m,plo,phi,name,lam,a,c,depth+1);certify_mix(labels,m,thi,plo,phi,name,lam,a,c,depth+1)
    else:
        m=(plo+phi)/2;certify_mix(labels,tlo,thi,plo,m,name,lam,a,c,depth+1);certify_mix(labels,tlo,thi,m,phi,name,lam,a,c,depth+1)


def verify_analytic_exceptions(boxes):
    # Mixed critical prefix:
    #   S = -lambda p_phi + (1-lambda) p_(theta+15),
    # so ||S|| <= lambda + (1-lambda) = 1.
    for name, _, _, lam, _, _ in boxes:
        if not (F(0) < lam < F(1)):
            raise AssertionError((name, "critical-prefix weight failure"))

    # In special boxes a = -sqrt(3) lambda.  Before B,
    # W/lambda = (-cos phi, sqrt(3)(cos phi-1)); for 0<=phi<=60,
    # ||W/lambda||^2 - 1 = 2(2 cos phi-1)(cos phi-1) <= 0.
    special = [row for row in boxes if row[4] == "special"]
    if {row[0] for row in special} != {"T1P1", "T2P1"}:
        raise AssertionError("unexpected special-box set")
    for name, _, (plo, phi), lam, _, _ in special:
        if not (F(0) <= plo < phi <= F(60)):
            raise AssertionError((name, "special-prefix phi range"))
        if not (F(0) < lam < F(1)):
            raise AssertionError((name, "special-prefix lambda"))


def main(data_path):
    import time
    payload, boxes = parse_coefficient_data(data_path)
    verify_analytic_exceptions(boxes)

    start = time.time()

    # Sanity-check the exact zero-sum block construction and positivity of
    # all strict-escape multipliers.  Interval dependence can widen the sum,
    # but zero must remain enclosed in both coordinates.
    full_mix_labels = frozenset(("A", "B", "C", "D", "E", "F", "XR", "XVm", "XVp"))
    for name, theta, phi, lam, a, c in boxes:
        total = coeff_mix(full_mix_labels, theta, phi, lam, a, c)
        if not (total[0].lo <= 0 <= total[0].hi and total[1].lo <= 0 <= total[1].hi):
            raise AssertionError(("mixed block does not enclose zero", name, total))
        mu = (I(1) - I(lam)) / (2 * RHO)
        if mu.lo <= 0:
            raise AssertionError((name, "nonpositive cap multiplier"))

    critical = frozenset(("A", "B", "C", "XVm", "XVp"))
    special_prefix = frozenset(("A", "B", "XVm", "XVp"))
    mixed_checks = 0
    for name, (tlo, thi), (plo, phi), lam, a, c in boxes:
        for labels in MIX_PREFIXES:
            if labels == critical:
                continue
            if a == "special" and labels == special_prefix:
                continue
            certify_mix(labels, tlo, thi, plo, phi, name, lam, a, c)
            mixed_checks += 1

    # Independent grammar counts.
    if len(MIX_WORDS) != 12 or len(MIX_PREFIXES) != 13:
        raise AssertionError("mixed grammar count mismatch")

    elapsed = time.time() - start
    print("schema                     =", payload["schema"])
    print("pi interval                =", PI)
    print("mixed orders               =", len(MIX_WORDS))
    print("mixed unique prefixes      =", len(MIX_PREFIXES))
    print("mixed coefficient boxes    =", len(boxes))
    print("mixed interval checks      =", mixed_checks)
    print("mixed analytic checks      =", len(boxes) + sum(a == "special" for _, _, _, _, a, _ in boxes))
    print("adaptive interval nodes    =", stats["nodes"])
    print("maximum recursion depth    =", stats["maxdepth"])
    print("point trig cache entries   =", len(PC))
    print("interval trig cache entries=", len(IC))
    print("elapsed seconds            =", round(elapsed, 3))
    print("PASS: all Case-4 mixed-state prefix norms are at most 1")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate_data", nargs="?", default="case4_mixed_certificate_data.json")
    args = parser.parse_args()
    main(args.certificate_data)
