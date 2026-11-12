#!/usr/bin/env python3
"""Fixed rational arithmetic for the 16 noncritical localized cap-cap arcs.
Only the constants sqrt(3) and b=2*cos(15 deg)-1 are enclosed. There is
no interval of the geometric parameter and no angular subdivision.
"""
from fractions import Fraction as F
from pathlib import Path
import json,itertools,sys

if sys.flags.optimize:
    raise SystemExit("ERROR: Python optimization disables required proof checks; run without -O/-OO or PYTHONOPTIMIZE.")
class I:
 def __init__(self,a,b=None):self.lo=F(a);self.hi=F(a if b is None else b)
 def __add__(self,o):o=asI(o);return I(self.lo+o.lo,self.hi+o.hi)
 __radd__=__add__
 def __neg__(self):return I(-self.hi,-self.lo)
 def __sub__(self,o):return self+-asI(o)
 def __rsub__(self,o):return asI(o)+-self
 def __mul__(self,o):
  o=asI(o);p=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi];return I(min(p),max(p))
 __rmul__=__mul__
 def __truediv__(self,o):o=asI(o);assert o.lo>0;return self*I(1/o.hi,1/o.lo)
 def sq(self):return I(0 if self.lo<=0<=self.hi else min(self.lo**2,self.hi**2),max(self.lo**2,self.hi**2))
def asI(a):return a if isinstance(a,I) else I(a)
def add(u,v):return[u[0]+v[0],u[1]+v[1]]
def mul(c,u):return[c*u[0],c*u[1]]
def cross(u,v):return u[0]*v[1]-u[1]*v[0]
def norm2(u):return u[0].sq()+u[1].sq()
s=I('1.73205','1.73206');b=I('.9318','.9319')
# Elementary proof of both constant bounds.
assert s.lo**2<3<s.hi**2
# b > 0 and (b+1)^2 = 2+sqrt(3).
assert (b.lo+1)**2<2+s.lo and (b.hi+1)**2>2+s.hi
E=[ [I(0),I(1)], [I(F(1,2)),-s/2] ]
ql=add(E[0],mul(b,E[1])) # kappa p22.5
qr=add(mul(b,E[0]),E[1]) # kappa p7.5
rots={150:(-s/2,I(F(1,2))),180:(I(-1),I(0)),210:(-s/2,I(F(-1,2)))}
def rot(v,t):
 c,ss=rots[t];return[c*v[0]-ss*v[1],ss*v[0]+c*v[1]]
def sign(a):
 assert a.lo>0 or a.hi<0,('unresolved fixed determinant',a.lo,a.hi)
 return 1 if a.lo>0 else -1
words=[tuple(pre)+tuple('BCDE')+tuple(tail) for pre in [('A','XV'),('XV','A')] for tail in [('XU','F'),('F','XU')]]
prefixes={frozenset(w[:j]) for w in words for j in range(1,len(w))}
assert len(words)==4 and len(prefixes)==9
constant_sets={frozenset(['A']),frozenset(['XV']),frozenset(['A','B','C','D','E','XU','XV']),frozenset(['A','B','C','D','E','F','XV'])}
assert prefixes-constant_sets=={frozenset(list('ABCDE'[:j])+['XV']) for j in range(1,6)}
# The four coefficient pairs of Appendix A.2, in exact rational form.
splits={3:[('-13/20','-5/16'),('-1/2','-1/2')],
        4:[('-3/4','-1/2'),('-1/2','-1/4')]}
records=[]
for case,eps in [(3,1),(4,0)]:
 for half_index,(name,lo,hi) in enumerate([('left_half',150,180),('right_half',180,210)]):
  half={'name':name}
  left=half_index==0
  a,c=map(I,splits[case][half_index])
  qu,qv=(ql,qr) if left else (qr,ql)
  hu,hv=([-I(1),-b],[-b,-I(1)]) if left else ([-b,-I(1)],[-I(1),-b])
  # The four remaining localized prefixes have norms |a|, |qv|,
  # |hv[eps]-c| and |qu|, bounded here by 3/4.
  assert a.sq().hi<=F(9,16) and (hv[eps]-c).sq().hi<=F(9,16)
  assert norm2(qu).hi<F(9,16) and norm2(qv).hi<F(9,16)
  tests=[('AXV',mul(a,E[0]),qv),
         ('ABXV',add(mul(a,E[0]),mul(hu[1],E[1])),qv),
         ('ABCDXV',mul(-1,qu),add(qv,mul(c,E[eps]))),
         ('ABCDEXV',mul(-1,qu),mul(-(hv[eps]-c),E[eps]))]
  for name,u,v in tests:
   sg=[sign(cross(rot(v,lo),u)),sign(cross(u,rot(v,hi)))]
   m=F(999,1000)
   assert F(0)<m<F(1), 'positive circle bound below one'
   if sg==[1,1]:
    method='alignment';U=norm2(u);V=norm2(v);A=I(m*m)-U-V
    assert A.lo>0 and (A.sq()-4*U*V).lo>0,(case,half['name'],name,'triangle')
   else:
    assert sg in [[1,-1],[-1,1]]
    method='upper endpoint' if sg==[1,-1] else 'lower endpoint'
    assert norm2(add(u,rot(v,lo))).hi<m*m,(case,name,'lo')
    assert norm2(add(u,rot(v,hi))).hi<m*m,(case,name,'hi')
   records.append({'case':case,'half':half['name'],'prefix':name,'maximum_location':method,'determinant_signs':sg,'norm_bound':'999/1000'})
assert len(records)==16, 'complete sixteen-comparison inventory'
result={'status':'PASS','comparisons':len(records),'parameter_interval_arithmetic':False,
        'fixed_constant_enclosures':{'sqrt3':['1.73205','1.73206'],'b':['0.9318','0.9319']},
        'all_noncritical_localized_prefixes_bound':'999/1000',
        'constant_prefixes_bound':'3/4','critical_prefix':'analytic chord bound 1',
        'records':records}
(Path(__file__).resolve().parent/'capcap_circle_check.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS: all 16 noncritical circular arcs have norm < 999/1000; four constant-prefix forms <= 3/4. Critical chord bound remains exactly 1.')
