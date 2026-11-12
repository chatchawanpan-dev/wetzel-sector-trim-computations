#!/usr/bin/env python3
"""Numerical search for the Case 2 coefficients of Section 5.3.

This is the search described in Section 1.1 of the paper.  It is not part
of the proof: the proof uses only the coefficient file
verification/Case2/case2_actual_mixed_prefix_data.json, which the exact
checker verifies without any optimizer.

Method.  The domain [150, 210] x [0, 30] (degrees) of (theta, phi) is cut
into sixteen rectangles.  On each rectangle a second-order cone program
chooses lambda in [0.02, 0.98], a and c so as to minimize the largest
prefix norm at the nine points of a 3 x 3 grid (solver Clarabel, absolute
and relative gap and feasibility tolerances 1e-9).  If the optimal value
is below 0.9995, the solution is rounded to nine decimal places and
passed to the exact checker, which must accept every prefix on the whole
rectangle; the rectangle is then written to the output file.  Otherwise
the rectangle is bisected along its longer side, to depth at most eight.
On the rectangles with phi = 0 that contain theta = 165 the relation
a = -sqrt(3) lambda of Lemma 5.5(i) is imposed.

The run of 6 September 2026 that produced the distributed coefficients is
recorded in case2_search_2026-09-06.log (33 rectangles after 17
bisections).  A new run need not reproduce the same numbers, because they
depend on the versions of the numerical software.

Usage (NumPy, CVXPY and Clarabel are needed):

    python3 find_case2_coefficients.py
    python3 ../verification/Case2/verify_case2_actual_mixed_prefix.py output/case2_coefficients_rerun.json
"""
import sys,math,json,itertools,time,argparse
from pathlib import Path
from fractions import Fraction as F
HERE=Path(__file__).resolve().parent
BASE=HERE.parent/'verification'/'Case2'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path,
                    default=HERE/'output'/'case2_coefficients_rerun.json',
                    help='output JSON path (default: output/case2_coefficients_rerun.json beside this file)')
args=parser.parse_args()
OUTPUT=args.output.expanduser().resolve()
if OUTPUT==(BASE/'case2_actual_mixed_prefix_data.json').resolve():
 parser.error('choose a separate output path; the distributed proof certificate must be preserved')
import numpy as np
import cvxpy as cp
sys.path.insert(0,str(BASE));import verify_case2_actual_mixed_prefix as rig

def p(a):return np.array([math.cos(math.radians(float(a))),math.sin(math.radians(float(a)))])
rho=p(7.5)[0];rt3=math.sqrt(3)
rows=[];start=time.time()

def solve_box(tlo,thi,plo,phi,depth=0):
 special=plo==0 and tlo<=165<=thi
 lam=cp.Variable();c=cp.Variable();a=-rt3*lam if special else cp.Variable();t=cp.Variable()
 mu=(1-lam)/(2*rho)
 cons=[lam>=.02,lam<=.98]
 for th in [tlo,(tlo+thi)/2,thi]:
  for ph in [plo,(plo+phi)/2,phi]:
   alpha=-lam*(p(ph)[1]+rt3*p(ph)[0]);beta=-2*lam*p(ph)[0];gamma=-2*(1-lam)*p(15)[0]
   coeff={'A':a*p(90),'B':beta*p(-60),'D':(alpha-a)*p(90),'XR':lam*p(ph),'C':c*p(th+90),'E':gamma*p(th-60),'F':(gamma-c)*p(th+90),'XVm':mu*p(th+F(15,2)),'XVp':mu*p(th+F(45,2))}
   for labels in rig.PREFIXES:
    if special and labels==rig.EXCEPTION:continue
    cons.append(cp.norm(sum(coeff[x] for x in labels))<=t)
 problem=cp.Problem(cp.Minimize(t),cons);problem.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9)
 if t.value is not None and t.value<.9995:
  lv=F(f'{lam.value:.9f}');av='special' if special else F(f'{a.value:.9f}');cv=F(f'{c.value:.9f}')
  try:
   for labels in rig.PREFIXES:
    if special and labels==rig.EXCEPTION:continue
    rig.certify(labels,(tlo,thi),(plo,phi),lv,av,cv)
  except AssertionError as e:
   print('exact failed',tlo,thi,plo,phi,str(e)[:120],flush=True)
  else:
   row={'name':f'B{len(rows)+1:03}','theta':[str(tlo),str(thi)],'phi':[str(plo),str(phi)],'lambda':str(lv),'a':'-sqrt(3)*lambda' if special else str(av),'c':str(cv)}
   rows.append(row);print('CERTIFIED',row,'sampled t',t.value,'elapsed',round(time.time()-start,1),flush=True)
   data={'schema':'case2-actual-mixed-prefix-v1','theta_domain':['150','210'],'phi_domain':['0','30'],'boxes':rows}
   OUTPUT.parent.mkdir(parents=True,exist_ok=True)
   OUTPUT.write_text(json.dumps(data,indent=2)+'\n')
   return
 if depth>=8:raise RuntimeError(('deep synthesis',tlo,thi,plo,phi,t.value))
 print('SPLIT',tlo,thi,plo,phi,'sampled t',t.value,flush=True)
 if thi-tlo >= phi-plo:
  mid=(tlo+thi)/2;solve_box(tlo,mid,plo,phi,depth+1);solve_box(mid,thi,plo,phi,depth+1)
 else:
  mid=(plo+phi)/2;solve_box(tlo,thi,plo,mid,depth+1);solve_box(tlo,thi,mid,phi,depth+1)

for th in range(150,210,15):
 for ph in [F(0),F(15,2),F(15),F(45,2)]:
  solve_box(F(th),F(th+15),ph,ph+F(15,2))
print('DONE',len(rows),'boxes','elapsed',time.time()-start)
