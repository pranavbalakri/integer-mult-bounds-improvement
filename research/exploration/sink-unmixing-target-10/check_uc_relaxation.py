"""Reject beating the old borrowed saving using only U_c retirement.

This is scoped to the pinned first-stage producer and the old borrowed second
stage. Each first-stage role independently chooses its best common-point
anchor or no retirement. All reference roles and all extra X path costs are
optimistically omitted. Original-last-frame retirement is outside this bound.
"""
if not __debug__:raise RuntimeError('Run without -O')
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
from argparse import ArgumentParser
import json
import indexed_profile as ip
import indexed_adapter as ia
import frame_profiles as fp


def run(source):
    ip.arithmetic.check_source(source);M=ia.model(source);h=23;m=575;N=4073300
    threshold=Q(5488,10**8);ext=Counter({529:1,23:1});merged=Counter({1:1,21:1,531:1})
    @lru_cache(None)
    def weights(r):
        a,b=ip.arithmetic.log_bounds(Q(m,r))
        return r*ip.arithmetic.exp_bounds(threshold*a)[0],r*ip.arithmetic.exp_bounds(threshold*b)[1]
    def interval(profile):
        return tuple(sum(n*weights(r)[i] for r,n in profile.items()) for i in (0,1))
    original=ip.arithmetic.profile();rows=Counter(original['rows']);I=fp.identity(h)
    counts=Counter(tuple(M['serialized']['frames'][M['C']['hold'][s][-1]-1])
                   for s in range(M['R']) if s not in M['C']['out'])
    selected=Counter();not_retired=0;alt_anchors=0;local_delta=Counter()
    for (core,cover),multiplicity in sorted(counts.items()):
        A=fp.frame(h,core,cover);old=fp.profile(fp.sub(I,A));ip.change(old,ext)
        options=[(None,old)]
        for c in ia.bits(core):
            U=fp.frame(h,1<<c,(1<<h)-1);p=fp.profile(fp.sub(U,A));ip.change(p,merged)
            options.append((c,p))
        choices=[]
        for c,p in options:
            lo,hi=interval(p);choices.append((lo,hi,c,p))
        chosen=min(choices,key=lambda item:item[1])
        for option in choices:
            if option[3]!=chosen[3]:assert chosen[1]<option[0]
        c,p=chosen[2:]
        if c is None:not_retired+=multiplicity
        elif c!=(core&-core).bit_length()-1:alt_anchors+=multiplicity
        selected[str(c)]+=multiplicity
        ip.change(local_delta,old,-multiplicity);ip.change(local_delta,p,multiplicity)
    assert ip.mass(local_delta)==0
    ip.change(rows,local_delta,2300)
    ip.change(rows,Counter([525,25,1,1,23]),-N);W=original['W']-N
    assert all(n>=0 for n in rows.values())
    lo,hi=ip.arithmetic.moment(dict(rows=rows,m=m,W=W),threshold)
    assert lo>1
    # Direct X→full has [1,h−2]. Every nontrivial intermediate frame
    # yields at least two positive children summing h−1, whose concave
    # moment is at least that of [1,h−2]. A full auxiliary reference
    # path sums m and adds a strictly positive amount to moment−W.
    return dict(status='PASS scoped U_c-only retirement cannot improve the old borrowed saving',
                threshold=str(threshold),relaxed_moment_excess_lower=str(lo-1),
                approximate_root_of_threshold_minimizing_profile=ip.float_root(rows,m,W),
                eligible_roles=sum(counts.values()),not_retired_in_relaxation=not_retired,
                nonleast_point_anchors=alt_anchors,selected_anchor_counts=dict(selected),
                extra_X_path_costs_omitted=True,all_reference_costs_omitted=True,
                original_last_frame_retirement_outside_scope=True,
                W=W,child_multiplicities={str(r):n for r,n in sorted(rows.items()) if n})


if __name__=='__main__':
    p=ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);args=p.parse_args()
    result=run(args.source)
    Path(__file__).with_name('uc-relaxation-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('child_multiplicities',)},indent=2))
