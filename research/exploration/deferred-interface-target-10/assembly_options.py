"""Exact candidate macro inequalities for the credited deferred construction.

The changed finite frames and actual signed h24 guard have separate audits.
This is not a proof of the inherited all-size routing/compiler interfaces.
The two precision choices are alternatives and must not be combined.
"""
if not __debug__:
    raise RuntimeError('Run without -O')

from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys

from audit import HERE, HEAD, module, pins, read
from tradeoffs import halving


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, required=True)
    ap.add_argument('--check-guard', action='store_true', help='Replay the full actual signed h24 guard audit')
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    pin_report = pins(args.source)
    pkg = args.source/'research/deferred-signed'
    guard_dir = HERE.parent/'deferred-guard-target-10'
    guard = (module('actual_event_guard', guard_dir/'audit.py').run(args.source)
             if args.check_guard else read(guard_dir/'results.json'))
    assert guard == read(guard_dir/'results.json')
    assert guard['source_commit'] == HEAD and guard['C1'] == 1
    assert guard['recursive_guard_coefficient'] == 51
    assert guard['layer_guard'] == '128(d+1)'
    assert guard['precision']['coefficient'] == 768
    assert guard['scalar_prefix_product'] == '3876535381862'

    sys.path.insert(0, str(HERE.parent/'indexed-cycle-target-10'))
    macro = module('candidate_macro', HERE.parent/'indexed-cycle-target-10/combined_assembly.py')
    # Explicit parameter substitution in a newly loaded module. The older
    # h18 helper source and every old certificate stay unchanged.
    ac = Q(36926111, 500000000000)
    macro.AC = ac
    exact = module('independent_moment', HERE.parent/'indexed-cycle-target-10/independent_moment.py')
    cx_rows = {int(t): n for t, n in read(pkg/'round6-complex-literal-ledger/result.json')['full_histogram'].items()}
    _, cx_upper = exact.moment(dict(m=576, W=207387136, rows=cx_rows), ac)
    assert cx_upper < 1
    dr = module('deferred', pkg/'swapnil-round7/independent/deferred-readout/deferred.py')
    W, D = dr.load()
    S = dr.Schedule(W, D)
    original = S.f[:]
    cov, xd = S.adjoint(), S.xdata()
    old_profile = {int(t): n for t, n in read(pkg/'round7-literal-ledger/result.json')['histogram'].items()}
    cases = []
    for cap, ab in ((None, Q(31987,500000000)), (21,Q(637083,10**10)),
                    (20,Q(631911,10**10)), (19,Q(624868,10**10))):
        frame_receipt = None
        S.f = original[:] if cap is None else [min(f,cap) for f in original]
        dim = lambda k: S.f[k[1]] if k[0] == 'sigma' else S.dim(k)
        p = dr.histogram(S, dr.chain_ranks(S,dim), S.ylevels(cov), xd,
                         corner=dr.STAIRCASE_CORNER(S.h))
        if cap is None:
            assert p['hist'] == old_profile
        else:
            folder = 'deferred-truncation-target-10' if cap == 21 else 'deferred-cap-target-10'
            frame_path = HERE.parent/folder/f'cap{cap}-results.json'
            frames = read(frame_path)
            frame_receipt = dict(path=f'../{folder}/cap{cap}-results.json',
                                 sha256=sha256(frame_path.read_bytes()).hexdigest())
            assert frames['source_commit'] == '741e7aa078392553815df7926ee17ac5e25a8c38'
            for rel, digest in frames['source_sha256'].items():
                assert sha256((pkg/'swapnil-round7'/rel).read_bytes()).hexdigest() == digest
            assert {int(k):v for k,v in frames['entrance_dimension_histogram'].items()} == dict(Counter(S.f))
            assert frames['NE_failures'] == 0 and frames['all_caps_nondegenerate']
            assert frames['total_rank_mass_unchanged']
            assert frames['largest_auxiliary_block'] == max(p['hist'])
            if cap == 21:
                assert frames['all_changed_edges_nested_over_Q']
            else:
                assert frames['all_unchanged_predecessors_contained']
                assert frames['all_caps_inside_every_original_high_frame']
                assert frames['singleton_first_edge_rank_histogram_matches_unsplit_dimension_count']
        rows = p.pop('hist')
        profile = dict(m=p['m'], W=p['W'], rows=rows)
        assert sum(t*n for t,n in rows.items()) == 57403754177
        db, dc = halving(529,max(rows)), halving(576,max(cx_rows))
        coefficient = p['W'].bit_length()*db + (207387136).bit_length()*dc
        degree = 1000*((coefficient*51)//25000+1)
        bridge = dict(bit=dict(m=529,W=p['W'],maxchild=max(rows),halving_degree=db),
            complex=dict(m=576,W=207387136,maxchild=552,halving_degree=dc),
            product_row_coefficient=coefficient, rounded_row_degree=degree,
            strict_row_degree_gap=str(Q(degree)-Q(51,25)*coefficient),
            guard='128(d+1)')
        choices = []
        for name,eps,delta,gap,power in (
                ('larger saving', Q(9999,10000),Q(1,10**6),Q(1,10**11),Q(3)),
                ('lower precision and wider margins',Q(71,100),Q(9,100),Q(1,10**8),Q(43,20))):
            # At these parameters ab is below both complex/internal limits;
            # round down to a 1e-8 grid and let the full helper check
            # EVERY constraint and margin rather than just this bottleneck.
            upper = eps*(ab-2*gap)
            kap = Q((upper*10**8).__floor__(),10**8)
            result = macro.evaluate(name,ab,eps,delta,gap,kap,power,bit_profile=profile)
            assert kap < Q(1,1024)
            choices.append(result)
        cases.append(dict(entrance_cap=cap, bit_saving=str(ab),frame_receipt=frame_receipt,
                          bridge=bridge,choices=choices))
    result = dict(source_commit=HEAD,pin_report=pin_report,
        status='Exact candidate assembly inequalities; global multiplication interfaces remain inherited and unreviewed.',
        new_multiplication_bound=False,target_met=False,guard_replayed_in_this_run=args.check_guard,
        actual_guard_receipt='../deferred-guard-target-10/results.json',
        actual_guard_receipt_sha256=sha256((guard_dir/'results.json').read_bytes()).hexdigest(),
        actual_guard_snapshot=guard,
        complex_moment_gap=str(1-cx_upper),cases=cases,
        assumptions=guard['assumptions']+[
            'Changed bit-frame certificates and original deferred word/geometry checks',
            'Common rational basis, ordered residuals, compact routing and exact recovery',
            'The enlarged-precision macro stack, distinct from PR97 balanced assembly'])
    text = json.dumps(macro.serial(result),indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text)


if __name__ == '__main__':
    main()
