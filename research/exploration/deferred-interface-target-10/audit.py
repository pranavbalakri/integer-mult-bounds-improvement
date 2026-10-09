"""Reproduce credited PR97 finite checks without hiding its missing log pins.

Original audit wrapper, written with OpenAI Codex assistance. Apache-2.0.
The invoked source is Zhihao Chen's PR97 integration of Swapnil Jain's
round-seven bit and round-six complex constructions. See README.md.
No downloaded source file or its manifest is modified.
"""
from argparse import ArgumentParser
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

if not __debug__:
    raise RuntimeError('Run without -O')
HERE = Path(__file__).resolve().parent
HEAD = 'f5f9c56e637463cac1e300d1589ccf42838f688a'
MANIFEST = 'eca8ccc565dc6145b8da703a752f7dd54873993bb5bc7ec0c23bd6649a2ac1d0'
IMPORT = '4b5f15981e846084b6def674ce0903c2d4ef0c930ea0c40093a2a8584164102d'
MISSING = sorted('research/deferred-signed/validation-prior/' + name for name in
                 ('frames-replay.log', 'lifted-replay.log', 'stair-replay.log',
                  'stair_control-replay.log'))


def read(p):
    return json.loads(p.read_text())


def pins(source):
    pkg = source / 'research/deferred-signed'
    assert sha256((pkg/'SOURCE.json').read_bytes()).hexdigest() == MANIFEST
    assert sha256((pkg/'swapnil-round7/IMPORT.json').read_bytes()).hexdigest() == IMPORT
    reports = {}
    for label, root, manifest in (
            ('integration', source, pkg/'SOURCE.json'),
            ('imported', pkg/'swapnil-round7', pkg/'swapnil-round7/IMPORT.json')):
        missing, matched = [], 0
        for rel, digest in read(manifest)['files'].items():
            file = root/rel
            if not file.exists():
                missing.append(rel)
            else:
                assert sha256(file.read_bytes()).hexdigest() == digest, ('changed', rel)
                matched += 1
        assert sorted(missing) == (MISSING if label == 'integration' else [])
        reports[label] = dict(matched=matched, missing=sorted(missing), changed=0)
    return reports


def module(name, file):
    spec = importlib.util.spec_from_file_location(name, file)
    obj = importlib.util.module_from_spec(spec)
    sys.modules[name] = obj
    spec.loader.exec_module(obj)
    return obj


def moments(source):
    exact = module('independent_moment', HERE.parent/'indexed-cycle-target-10/independent_moment.py')
    pkg = source/'research/deferred-signed'
    results = {}
    for label, filename, key, m, W, mass, saving in (
            ('bit', 'round7-literal-ledger/result.json', 'histogram', 529,
             108516254, 57403754177, Q(31987, 500000000)),
            ('complex', 'round6-complex-literal-ledger/result.json', 'full_histogram',
             576, 207387136, 119453132304, Q(36926111, 500000000000))):
        rows = {int(k): v for k, v in read(pkg/filename)[key].items()}
        assert sum(k*v for k, v in rows.items()) == mass
        assert all(0 < k < m and v > 0 for k, v in rows.items())
        p = dict(m=m, W=W, rows=rows)
        lo, hi = exact.moment(p, saving)
        assert hi < 1
        reject_lo, _ = exact.moment(p, Q(1, 1024))
        assert reject_lo > 1
        results[label] = dict(m=m, W=W, mass=mass, deficit=m*W-mass,
            saving=str(saving), exact_acceptance_gap=str(1-hi),
            target_rejection_gap=str(reject_lo-1), maxchild=max(rows))
    return results


def stable(obj):
    if isinstance(obj, dict):
        return {k: stable(v) for k, v in obj.items() if k != 'elapsed'}
    if isinstance(obj, list):
        return [stable(v) for v in obj]
    return obj


def replay(source):
    source_pkg = source/'research/deferred-signed'
    jobs = [
        ('round7_literal_frame_ledger.py', 'round7-literal-ledger/result.json', 300),
        ('round7_tensor_endpoint_controls.py', 'round7-tensor-endpoint-controls.json', 180),
        ('round6_complex_interface_controls.py', 'round6-complex-interface-controls.json', 180),
        ('round6_complex_literal_ledger.py', 'round6-complex-literal-ledger/result.json', 300),
        ('check_complex_identity.py', None, 180),
        ('round7_balanced_assembly_candidate.py', 'round7-balanced-assembly-candidate.json', 180)]
    results = []
    with tempfile.TemporaryDirectory(prefix='kappa-deferred-audit-') as tmp:
        root = Path(tmp)
        pkg = root/'research/deferred-signed'
        shutil.copytree(source_pkg, pkg)
        (root/'scripts').mkdir()
        shutil.copy2(source/'scripts/structured_bulk_assembly.py', root/'scripts')
        for script, output, limit in jobs:
            proc = subprocess.run([sys.executable, script], cwd=pkg, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=limit)
            assert proc.returncode == 0, (script, proc.stdout[-5000:])
            if output:
                assert stable(read(pkg/output)) == stable(read(source_pkg/output)), script
            results.append(dict(script=script, returncode=0,
                output_matches_pinned=bool(output),
                stdout_sha256=sha256(proc.stdout.encode()).hexdigest()))
            print('Reproduced', script, file=sys.stderr, flush=True)
    return results


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, required=True, help='Pinned PR97 repository root')
    ap.add_argument('--replay', action='store_true')
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    result = dict(source_commit=HEAD, pins=pins(args.source), moments=moments(args.source),
        upstream_top_level_verifier_passes=False,
        upstream_failure='Four pinned historical validation logs are absent from this commit.',
        replay=replay(args.source) if args.replay else [],
        scope='Finite numerical and optional literal checks only; inherited analytic and tape assumptions remain.',
        new_multiplication_bound=False)
    text = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text)


if __name__ == '__main__':
    main()
