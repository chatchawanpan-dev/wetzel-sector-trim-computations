#!/usr/bin/env python3
"""Positive and rejection controls for the current Case 2 sine certificate.

Mutations apply to temporary copies only. Each textual change must match
exactly once, and the damaged checker must reject for the intended reason
under both normal and optimized Python. No production artifact is edited.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

BASE = Path(__file__).resolve().parent


def invoke(script, optimized):
    environment = dict(os.environ)
    environment.pop('PYTHONOPTIMIZE', None)
    command = [sys.executable] + (['-O'] if optimized else []) + [str(script)]
    return subprocess.run(command, cwd=script.parent, env=environment,
                          text=True, capture_output=True)


def main():
    sources = {name:(BASE/'Geometric'/filename).read_text(encoding='utf-8')
               for name,filename in [('sine','verify_case2_sine.py')]}
    controls = [
        ('sine','bad sine coefficient','a = -(2+x)*F(1,3)',
         'a = -(2+x)*F(1,4)','coefficient reversal identity'),
        ('sine','bad support coefficient',"'B':(-r,1-h)",
         "'B':(-r,2-h)",'U coefficient balance'),
        ('sine','invalid root isolation',"Interval('1.73205','1.73206')",
         "Interval('1.73207','1.73208')",'sqrt(3) isolation'),
        ('sine','insufficient constant precision',"Interval('.9318','.9319')",
         "Interval('.9','1')",'central norm square < 19/20'),
        ('sine','missing localized words',"product([('A','XV'),('XV','A')],[('XU','F'),('F','XU')])",
         "product([('A','XV')],[('XU','F'),('F','XU')])",'four words and nine prefix sets'),
        ('sine','target above one','M = F(49,50)','M = F(2)','positive norm bound strictly below one'),
        ('sine','negative global target','M = F(49,50)','M = F(-49,50)','positive norm bound strictly below one'),
        ('sine','unsupported stronger target','M = F(49,50)','M = F(9,10)','constant and central prefixes below M'),
        ('sine','negative quadratic curvature bound','(F(29,100),F(1,20),F(957,1000))',
         '(F(-29,100),F(1,20),F(957,1000))','valid completion parameters'),
        ('sine','missing quadratic bound row',
         'table = [(F(29,100),F(1,20),F(957,1000)),\n         (F(7,20),F(7,20),F(109,125))]',
         'table = [(F(29,100),F(1,20),F(957,1000))]','complete two-quadratic inventory'),
    ]
    records=[]
    with tempfile.TemporaryDirectory(prefix='h15-geometric-controls-') as directory:
        root=Path(directory)
        for name,source in sources.items():
            script=root/(name+'.py')
            script.write_text(source,encoding='utf-8')
            for optimized in [False,True]:
                result=invoke(script,optimized)
                if result.returncode or 'PASS:' not in result.stdout:
                    raise RuntimeError(('positive control failed',name,optimized,result.stdout,result.stderr))
                records.append({'kind':'positive','certificate':name,'optimized':optimized,'status':'PASS'})
                print(f'PASS positive: {name}, optimized={optimized}',flush=True)
        for name,title,old,new,expected in controls:
            source=sources[name]
            if source.count(old)!=1 or old==new:
                raise RuntimeError(('mutation did not match exactly once',title,source.count(old)))
            damaged=source.replace(old,new,1)
            if damaged==source:
                raise RuntimeError(('mutation was not applied',title))
            script=root/(name+'_damaged.py')
            script.write_text(damaged,encoding='utf-8')
            for optimized in [False,True]:
                result=invoke(script,optimized)
                if result.returncode==0 or expected not in result.stderr:
                    raise RuntimeError(('damaged certificate not rejected as expected',title,optimized,
                                        result.returncode,expected,result.stdout,result.stderr))
                records.append({'kind':'negative','certificate':name,'test':title,
                                'optimized':optimized,'status':'PASS: rejected','reason':expected})
                print(f'PASS rejection: {title}, optimized={optimized}: {expected}',flush=True)
        circle=(BASE/'Geometric'/'verify_capcap_circles.py').read_text(encoding='utf-8')
        script=root/'circle.py'
        script.write_text(circle,encoding='utf-8')
        result=invoke(script,False)
        if result.returncode or 'PASS:' not in result.stdout:
            raise RuntimeError(('circle positive control failed',result.stdout,result.stderr))
        records.append({'kind':'positive','certificate':'circle','optimized':False,'status':'PASS'})
        print('PASS positive: circle, optimized=False',flush=True)
        result=invoke(script,True)
        if result.returncode==0 or 'optimization disables required proof checks' not in result.stderr:
            raise RuntimeError(('circle optimization guard failed',result.stdout,result.stderr))
        records.append({'kind':'negative','certificate':'circle','test':'optimization guard','status':'PASS: rejected'})
        print('PASS rejection: circle refuses optimized execution',flush=True)
        for title,old,new,expected in [
            ('missing circle case','[(3,1),(4,0)]','[(3,1)]','complete sixteen-comparison inventory'),
            ('invalid circle target','m=F(999,1000)','m=F(2)','positive circle bound below one')]:
            if circle.count(old)!=1:
                raise RuntimeError(('circle mutation did not match exactly once',title))
            script.write_text(circle.replace(old,new,1),encoding='utf-8')
            result=invoke(script,False)
            if result.returncode==0 or expected not in result.stderr:
                raise RuntimeError(('circle rejection missing',title,result.stdout,result.stderr))
            records.append({'kind':'negative','certificate':'circle','test':title,'status':'PASS: rejected','reason':expected})
            print(f'PASS rejection: {title}: {expected}',flush=True)
    positive=sum(r['kind']=='positive' for r in records)
    negative=sum(r['kind']=='negative' for r in records)
    print(f'PASS: {positive} positive controls and {negative} rejection controls; distributed sources unchanged.',flush=True)
    return records


if __name__=='__main__':
    main()
