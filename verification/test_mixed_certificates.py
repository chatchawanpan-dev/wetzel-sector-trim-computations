#!/usr/bin/env python3
"""Portable rejection controls for mixed-state inventory, coverage, signs and norms.

Positive verification is performed by verify_proof_essential.py. These controls
mutate temporary copies only and require the intended rejection, not merely
any unsuccessful subprocess. Production proof arithmetic is unchanged.
"""
import copy
from fractions import Fraction
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent.parent
CASES = [
    (2, 'Case2/verify_case2_actual_mixed_prefix.py', 'Case2/case2_actual_mixed_prefix_data.json', 'boxes'),
    (3, 'Case3/verify_case3_mixed.py', 'Case3/case3_mixed_certificate_data.json', 'mixed_boxes'),
    (4, 'Case4/verify_case4_mixed.py', 'Case4/case4_mixed_certificate_data.json', 'mixed_boxes'),
]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    inventory = json.loads((ROOT / 'verification/PROOF_ESSENTIAL_FILES.json').read_text())
    inputs = [row['path'] for row in inventory['files'] if row['runtime_required']]
    env = dict(os.environ)
    env.pop('PYTHONOPTIMIZE', None)
    results = []
    for case, script, data_path, key in CASES:
        original = json.loads((ROOT / 'verification' / data_path).read_text())
        for mutation in ['missing_box', 'coverage_gap', 'negative_weight', 'wrong_coefficient']:
            altered = copy.deepcopy(original)
            rows = altered[key]
            if mutation == 'missing_box':
                rows.pop(0)
                expected = 'coverage gap' if case == 2 else 'unexpected mixed-box count'
            elif mutation == 'coverage_gap':
                for row in rows:
                    if Fraction(row['theta'][0]) < 151:
                        row['theta'][0] = '151'
                expected = 'coverage gap'
            elif mutation == 'negative_weight':
                rows[0]['lambda'] = '-1/10'
                expected = 'assert F(0) < lam < F(1)' if case == 2 else 'lambda'
            else:
                # Change the sign of c in the first rectangle.  A prefix norm then
                # exceeds one there, and the interval check must fail at its depth limit.
                rows[0]['c'] = rows[0]['c'].lstrip('-')
                expected = 'uncertified prefix' if case == 2 else "'MIX'"
            require(altered != original, 'Mutation did not change input')
            with tempfile.TemporaryDirectory(prefix='h15-mixed-control-') as temporary:
                base = Path(temporary)
                for name in inputs:
                    target = base / name
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(ROOT / name, target)
                (base / 'verification' / data_path).write_text(json.dumps(altered))
                command = [sys.executable, str(base / 'verification' / script),
                           str(base / 'verification' / data_path)]
                run = subprocess.run(command, cwd=temporary, env=env, capture_output=True, text=True, timeout=60)
                output = run.stdout + run.stderr
                require(run.returncode != 0 and expected in output,
                        'Intended rejection missing: Case %s %s\n%s' % (case, mutation, output))
                results.append({'case': case, 'mutation': mutation, 'expected': expected, 'status': 'PASS'})
                print('PASS: Case %s rejects %s (%s)' % (case, mutation, expected))
    require(len(results) == 12, 'Expected twelve rejection controls')
    print('PASS: all twelve mixed-state rejection controls; distributed inputs unchanged')
    return results


if __name__ == '__main__':
    main()
