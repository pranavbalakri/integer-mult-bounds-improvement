"""Conditional dimension-profile accounting for common entrance truncation.

This does NOT verify new frame geometry or NE profiles. Each q>0 profile
assumes those obligations have been established separately. q=0 must exactly
reproduce the credited Jain/Chen pinned histogram. See README.md.
"""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as Q
import json
from math import floor
from pathlib import Path
import sys

from audit import HERE, HEAD, module, pins, read

if not __debug__:
    raise RuntimeError('Run without -O')


def halving(m, r):
    d = 1
    while m**d <= 2*r**d:
        d += 1
    return d


def bracket(p, exact):
    def approx(a):
        return sum(n*t/(p['m']*p['W'])*(p['m']/t)**a for t, n in p['rows'].items())
    lo, hi = 0., .002
    assert approx(lo) < 1 < approx(hi)
    for _ in range(50):
        mid = (lo+hi)/2
        if approx(mid) < 1:
            lo = mid
        else:
            hi = mid
    lower = Q(floor(lo*10**10), 10**10)
    upper = lower + Q(1, 10**10)
    _, accept = exact.moment(p, lower)
    reject, _ = exact.moment(p, upper)
    assert accept < 1 < reject
    return dict(lower=str(lower), upper=str(upper), lower_decimal=float(lower),
        acceptance_gap=str(1-accept), rejection_gap=str(reject-1))


def main():
    ap = ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, required=True)
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    pin_result = pins(args.source)
    pkg = args.source/'research/deferred-signed'
    deferred = module('deferred', pkg/'swapnil-round7/independent/deferred-readout/deferred.py')
    exact = module('independent_moment', HERE.parent/'indexed-cycle-target-10/independent_moment.py')
    witness, data = deferred.load()
    S = deferred.Schedule(witness, data)
    original_f = S.f[:]
    cov, xdata = S.adjoint(), S.xdata()
    expected = {int(t): n for t, n in read(pkg/'round7-literal-ledger/result.json')['histogram'].items()}
    answers = []
    choices = [('common_cut', q) for q in (0, 1, 2, 3, 4, 8, 12, 22)] + [('cap', cap) for cap in (21, 20, 19)]
    for mode, q in choices:
        S.f = ([max(f-q, 0) for f in original_f] if mode == 'common_cut'
               else [min(f, q) for f in original_f])
        def dimension(key):
            return S.f[key[1]] if key[0] == 'sigma' else S.dim(key)
        ranks = deferred.chain_ranks(S, dimension)
        p = deferred.histogram(S, ranks, S.ylevels(cov), xdata,
            corner=deferred.STAIRCASE_CORNER(S.h))
        if mode == 'common_cut' and q == 0:
            assert p['hist'] == expected
        assert sum(t*n for t, n in p['hist'].items()) == p['s'] == 57403754177
        p['rows'] = p.pop('hist')
        interval = bracket(p, exact)
        d = halving(p['m'], max(p['rows']))
        # This is the retained product-row bound, not a runtime estimate.
        cx_d = halving(576, 552)
        coefficient = p['W'].bit_length()*d + (207387136).bit_length()*cx_d
        degree = 1000*((coefficient*51)//25000+1)
        answers.append(dict(mode=mode, parameter=q, saving_bracket=interval,
            maxchild=max(p['rows']), halving_degree=d,
            product_row_coefficient=coefficient, rounded_row_degree=degree,
            entrance_dimensions=dict(sorted(Counter(S.f).items())),
            histogram=dict(sorted(p['rows'].items())),
            exactly_matches_pinned_original_profile=(mode == 'common_cut' and q == 0)))
        print('Computed', mode, q, file=sys.stderr, flush=True)
    result = dict(source_commit=HEAD, pin_report=pin_result, results=answers,
        status='q>0 arithmetic candidates only; realizable common intersections and fixed-basis corner profiles require separate verification.',
        new_multiplication_bound=False)
    text = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text)


if __name__ == '__main__':
    main()
