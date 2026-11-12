#!/usr/bin/env python3
"""Exact certificate: one sine coefficient rule for Case 2 cap--cap.

Standard library only. Algebra uses formal polynomials in x,y,r,h and
rational intervals for sqrt(3) and b. The analytic domain is
-1/2 <= x <= 1/2, y >= 0, x*x+y*y=1. No sampled angles are used.
The canonical geometric incidences and localized four words remain inputs.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json


def require(test, message):
    if not test:
        raise ValueError(message)


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F
    def __init__(self, lo, hi=None):
        lo = F(lo)
        hi = lo if hi is None else F(hi)
        require(lo <= hi, 'valid interval')
        object.__setattr__(self, 'lo', lo)
        object.__setattr__(self, 'hi', hi)
    def __add__(self, z):
        z = iv(z)
        return Interval(self.lo+z.lo, self.hi+z.hi)
    __radd__ = __add__
    def __neg__(self):
        return Interval(-self.hi, -self.lo)
    def __sub__(self, z):
        return self + -iv(z)
    def __rsub__(self, z):
        return iv(z) + -self
    def __mul__(self, z):
        z = iv(z)
        ps = [a*b for a in (self.lo, self.hi) for b in (z.lo, z.hi)]
        return Interval(min(ps), max(ps))
    __rmul__ = __mul__
    def square(self):
        return Interval(0 if self.lo <= 0 <= self.hi else min(self.lo**2, self.hi**2),
                        max(self.lo**2, self.hi**2))


def iv(z):
    return z if isinstance(z, Interval) else Interval(z)


# Formal polynomial arithmetic in (x,y,r,h); exact Fraction coefficients.
class Poly:
    def __init__(self, terms=0):
        if isinstance(terms, Poly):
            self.terms = terms.terms.copy()
        elif isinstance(terms, dict):
            self.terms = {e:F(c) for e,c in terms.items() if c}
        else:
            self.terms = {(0,0,0,0):F(terms)} if terms else {}
    def __add__(self, z):
        z = Poly(z)
        terms = self.terms.copy()
        for e,c in z.terms.items():
            terms[e] = terms.get(e, F(0)) + c
        return Poly(terms)
    __radd__ = __add__
    def __neg__(self):
        return Poly({e:-c for e,c in self.terms.items()})
    def __sub__(self, z):
        return self + -Poly(z)
    def __rsub__(self, z):
        return Poly(z) + -self
    def __mul__(self, z):
        z = Poly(z)
        terms = {}
        for e,c in self.terms.items():
            for f,d in z.terms.items():
                g = tuple(a+b for a,b in zip(e,f))
                terms[g] = terms.get(g, F(0)) + c*d
        return Poly(terms)
    __rmul__ = __mul__
    def __pow__(self, n):
        require(isinstance(n, int) and n >= 0, 'polynomial exponent')
        result = Poly(1)
        for _ in range(n):
            result = result*self
        return result
    def __eq__(self, z):
        return self.terms == Poly(z).terms
    def substitute(self, index, value):
        result = Poly(0)
        value = Poly(value)
        for e,c in self.terms.items():
            f = list(e)
            power = f[index]
            f[index] = 0
            result = result + Poly({tuple(f):c})*value**power
        return result
    def circle_reduce(self):
        # Replace y^(2k) by (1-x^2)^k; leave at most one factor y.
        result = Poly(0)
        for e,c in self.terms.items():
            f = list(e)
            power = f[1]//2
            f[1] %= 2
            result = result + Poly({tuple(f):c})*(1-x*x)**power
        return result
    def evaluate(self, values):
        result = Interval(0)
        for e,c in self.terms.items():
            term = Interval(c)
            for value,power in zip(values,e):
                for _ in range(power):
                    term = term*value
            result = result+term
        return result


x = Poly({(1,0,0,0):1})
y = Poly({(0,1,0,0):1})
r = Poly({(0,0,1,0):1})
h = Poly({(0,0,0,1):1})
zero = Poly(0)
e0 = (zero, Poly(1))
Q = (r,h)
a = -(2+x)*F(1,3)
c = -(1+x)*F(1,3)


def add(u,v):
    return tuple(a+b for a,b in zip(u,v))


def scale(t,u):
    return tuple(t*a for a in u)


def rotation(u):
    return (-y*u[0]-x*u[1], x*u[0]-y*u[1])


def sumvec(us):
    result = (zero, zero)
    for u in us:
        result = add(result,u)
    return result


def norm2(u):
    return (u[0]*u[0]+u[1]*u[1]).circle_reduce()


# Here b e1=(r,h-1), Q=e0+b e1.
coefficients = {
    'A':scale(a,e0), 'B':(-r,1-h), 'D':scale(-1-a,e0), 'XU':Q,
    'C':scale(c,rotation(e0)), 'E':rotation((-r,1-h)),
    'F':scale(-1-c,rotation(e0)), 'XV':rotation(Q),
}
require(sumvec(coefficients[k] for k in ('A','B','D','XU')) == (zero,zero),
        'U coefficient balance')
require(sumvec(coefficients[k] for k in ('C','E','F','XV')) == (zero,zero),
        'V coefficient balance')
words = [tuple(left)+tuple('BCDE')+tuple(right)
         for left,right in product([('A','XV'),('XV','A')],[('XU','F'),('F','XU')])]
prefixes = {frozenset(w[:j]) for w in words for j in range(1,len(w))}
full = frozenset(coefficients)
constant_sets = {frozenset(['A']),frozenset(['XV']),full-{'F'},full-{'XU'}}
require(len(words)==4 and len(prefixes)==9, 'four words and nine prefix sets')
require(prefixes-constant_sets == {frozenset(list('ABCDE'[:j])+['XV']) for j in range(1,6)},
        'five moving prefix sets')
expected = [add(scale(a,e0),rotation(Q)),
            add(add(scale(a,e0),(-r,1-h)),rotation(Q)),
            add(add(scale(a,e0),(-r,1-h)),rotation(add(Q,scale(c,e0)))),
            add(scale(-1,Q),rotation(add(Q,scale(c,e0)))),
            add(scale(-1,Q),scale(1+c,rotation(e0)))]
for j,p in enumerate(expected,1):
    require(sumvec(coefficients[k] for k in list('ABCDE'[:j])+['XV']) == p,
            f'prefix {j} binding')
require(norm2(coefficients['A']) == a*a, 'A singleton prefix')
require(norm2(coefficients['XV']) == r*r+h*h, 'XV singleton prefix')
require(norm2(sumvec(coefficients[k] for k in full-{'F'})) == (1+c)**2,
        'F complement prefix')
require(norm2(sumvec(coefficients[k] for k in full-{'XU'})) == r*r+h*h,
        'XU complement prefix')
require(a.substitute(0,-x)==-1-c and c.substitute(0,-x)==-1-a,
        'coefficient reversal identity')
N = [norm2(p) for p in expected]
require(N[4]==N[0].substitute(0,-x), 'fifth norm reverses first')
require(N[3]==N[1].substitute(0,-x), 'fourth norm reverses second')
K = r*r+(h-F(1,3))**2
require(N[2] == 2*K*(1+y)+2*x*x*(1-y)*F(1,9)-4*r*x*x*F(1,3),
        'central norm identity')
A = [F(1,9)-h*F(1,2)-r*F(2,3),
     -h*h+h*F(1,2)-r*r-r*F(2,3)+F(1,9)]
B = [h*F(2,3)-r*F(4,3)+F(4,9),
     h*F(4,3)+r*F(2,3)-F(2,9)]
C = [r*r+(h+F(2,3))**2,
     4*h*h-h*F(4,3)+4*r*r+F(1,9)]
Y = [2*h*(2+x)*F(1,3),
     2*(r*r+h*h-h*F(1,3)+h*x*F(1,3))]
for j in range(2):
    cubic = N[j].substitute(1,1-x*x*F(1,2))
    require(N[j]-cubic == Y[j]*(y-1+x*x*F(1,2)),
            f'norm {j+1}: y upper bound identity')
    quadratic = A[j]*x*x+B[j]*x+C[j]
    require(cubic-quadratic == -h*x*x*(x+F(1,2))*F(1,3),
            f'norm {j+1}: cubic remainder identity')

# Fixed algebraic constants.  These are exactly the same enclosures as
# in the four-block proposal, so no additional algebraic numbers enter.
s_i = Interval('1.73205','1.73206')
b_i = Interval('.9318','.9319')
require(s_i.lo**2<3<s_i.hi**2, 'sqrt(3) isolation')
require((b_i.lo+1)**2<2+s_i.lo and (b_i.hi+1)**2>2+s_i.hi, 'b isolation')
r_i = b_i*F(1,2)
h_i = 1-s_i*r_i
values = (Interval(F(-1,2),F(1,2)), Interval(0,1), r_i,h_i)
require(F(-1,2)<=values[0].lo<=values[0].hi<=F(1,2) and
        0<=values[1].lo<=values[1].hi<=1, 'analytic domain for x and y')
require(values[0].lo == -values[0].hi, 'domain invariant under x reversal')
require((1-x*x*F(1,2))**2-(1-x*x) == x**4*F(1,4),
        'squaring identity for y <= 1-x^2/2 on the upper unit semicircle')
require((1-x*x*F(1,2)).evaluate(values).lo>0, 'positive squared upper comparison')
require((x+F(1,2)).evaluate(values).lo>=0, 'sign of cubic remainder')
require(r_i.lo>0 and h_i.lo>0, 'positive r and h')
require((r*r+h*h-h*F(1,2)).evaluate(values).lo>0,
        'nonnegative y coefficient for second norm on whole x interval')
require(K.evaluate(values).lo>F(1,36),
        'nonnegative y coefficient for central norm on whole x interval')
require((4*K).evaluate(values).hi<F(19,20), 'central norm square < 19/20')
require((r*r+h*h).evaluate(values).hi<F(13,50), 'Q norm square < 13/50')
ai = a.evaluate(values)
ci = (1+c).evaluate(values)
require(F(-5,6)<=ai.lo<=ai.hi<=F(-1,2), 'constant prefix a range')
require(F(1,2)<=ci.lo<=ci.hi<=F(5,6), 'constant prefix 1+c range')
M = F(49,50)
require(0<M<1, 'positive norm bound strictly below one')
require(F(25,36)<M*M and F(13,50)<M*M and F(19,20)<M*M,
        'constant and central prefixes below M')
# Exact upper bounds for the two concave quadratics.  The entries mean
# A<=-alpha, |B|<=beta, C<=gamma.  Completing the square gives the result.
table = [(F(29,100),F(1,20),F(957,1000)),
         (F(7,20),F(7,20),F(109,125))]
require(len(table)==len(A)==len(B)==len(C)==2,
        'complete two-quadratic inventory')
records = []
for j,(alpha,beta,gamma) in enumerate(table):
    require(alpha>0 and beta>=0, f'quadratic {j+1}: valid completion parameters')
    ai,bi,ci = [p.evaluate(values) for p in (A[j],B[j],C[j])]
    require(ai.hi<=-alpha, f'quadratic {j+1}: negative leading coefficient')
    require(-beta<=bi.lo and bi.hi<=beta, f'quadratic {j+1}: linear coefficient')
    require(ci.hi<=gamma, f'quadratic {j+1}: constant coefficient')
    upper = gamma+beta*beta/(4*alpha)
    require(upper<M*M, f'quadratic {j+1}: global completed-square upper bound')
    records.append({'prefix':j+1,'alpha':str(alpha),'beta':str(beta),'gamma':str(gamma),
                    'squared_norm_upper_bound':str(upper)})
output = {'status':'PASS','construction':'one sine coefficient rule over the full orientation interval',
          'coefficients':{'a':'-(2+sin(theta))/3','c':'-(1+sin(theta))/3'},
          'orientation_degrees':[150,210],'words':4,'prefix_sets':9,
          'prefix_norm_bound':str(M),'closed_chain_length_lower_bound':str(1/M),
          'angular_sampling':False,'orientation_subdivision':False,
          'verification':'formal rational polynomial identities and exact fixed-algebraic bounds',
          'quadratic_bounds':records,'central_squared_norm_upper_bound':'19/20',
          'dependencies':['canonical Case 2 incidences','four-word upper-cap witness localization',
                          'balanced-prefix length inequality',
                          'y <= 1-x^2/2 when y>=0 and x^2+y^2=1',
                          'concave quadratic maximum by completing the square']}
Path(__file__).with_name('case2_sine_check.json').write_text(json.dumps(output,indent=2)+'\n')
print('PASS: formal prefix, reversal, and quadratic identities; exact fixed-algebraic bounds.')
print('PASS: one coefficient rule on [150,210]; all nine prefix norms < 49/50.')
print('PASS: closed-model chain length >= 50/49 > 1; no angular subdivision.')
