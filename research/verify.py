"""Reproduce the extension with Python 3.10+ and the standard library only.

Pinned third-party source copies are checked against their SHA-256 manifest.
No network access, package installation, or external proof assistant is needed.
Passing these finite checks is not verification of the complete theorem.
"""
import importlib.util
import json
from pathlib import Path
import sys
import unittest
from fractions import Fraction as Q

ROOT = Path(__file__).resolve().parent
for sub in ('scripts','independent/complex-twostage','independent/bit-improvement','independent/guard-improvement'):
    sys.path.insert(0,str(ROOT/sub))

from certificate import certificate
from check_paired_loo import check, check_leave_one_out
from reassociate import histogram
from assembly_options import AB, AC, options
from semantic_guard import certify as guard_certificate
from audit import matrix_identity, scalar_envelope, dirty_scalar_network
from hist import build
from frames import Checker
from cert import cert
from test_extension import Extension
from test_semantic_guard import SemanticGuard


def main():
    suite = unittest.TestSuite()
    for cls in (Extension,SemanticGuard):
        suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(cls))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful(): raise SystemExit(1)
    print('Checking new bit circuit and exact moment...',flush=True)
    bit = certificate()
    assert Q(bit['saving']) == AB
    bit_checks = check(23,histogram)
    bit_checks['leave_one_out_cases'] = check_leave_one_out()
    print('Checking smaller complex circuit, every frame, and exact moment...',flush=True)
    cx = build(16,copied=True)
    producer = cx.pop('c')
    frames = Checker(producer).run()
    assert frames['bad']==0 and cx['sum_ok']
    ac, slack = cert(cx['hist'],cx['W'],cx['m'])
    assert ac==AC and slack>0
    cx.update(saving=str(ac),exact_moment_slack=str(slack))
    audit = dict(all_pairs=matrix_identity(16), scalar_prefix=scalar_envelope(16),
                 dirty_scalar_network=dirty_scalar_network())
    guard = guard_certificate(16)
    assert guard['forward_G']==audit['scalar_prefix']['bounds'][0]['peak']
    assert guard['inverse_G']==audit['scalar_prefix']['bounds'][1]['peak']
    record, balanced = options()
    manifest = json.loads((ROOT/'upstream-manifest.json').read_text())
    output = dict(status='proposed conditional result; inherited theorem interfaces assumed',
                  source_repository=manifest['repository'], source_commit=manifest['commit'],
                  bit=bit, complex=cx, guard=guard, profiles=[record,balanced],
                  checks=dict(unit_tests=result.testsRun,bit=bit_checks,complex_frames=frames,scalar=audit))
    target=ROOT.parent/'certificates'/'latest.json'
    target.parent.mkdir(exist_ok=True)
    target.write_text(json.dumps(output,indent=2)+'\n')
    (ROOT/'independent'/'bit-improvement'/'assembly-options.json').write_text(json.dumps([record,balanced],indent=2)+'\n')
    print('PASS: exact witnesses',record['kappa'],'and',balanced['kappa'])
    print('PASS: complex roles',cx['W'],'guard',guard['layer_guard'],'precision p=768*b^3')


if __name__=='__main__': main()
