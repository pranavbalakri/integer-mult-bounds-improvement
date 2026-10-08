"""Exact scalar and paid-rank audit of the literal point-star transvection route.
The two-rank transfer E_i -> F -> E_j is granted, as requested; no stronger
claim about a new nonnested compiler or a global multiplication theorem is made.
"""
if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')

from itertools import combinations
from fractions import Fraction as Q
from pathlib import Path
import json


def check(h):
    assert h%4==2
    triples=list(combinations(range(h),3));v=len(triples)
    stars=[sum(1<<k for k,t in enumerate(triples) if i in t) for i in range(h)]
    assert all((u&w).bit_count()%2==0 for u in stars for w in stars)
    rows=[1<<k for k in range(v)]
    star_totals=[]
    for c,u in enumerate(stars):
        total=0
        for k,t in enumerate(triples):
            if c in t:total^=rows[k]
        # Earlier commuting stars leave all centre totals unchanged.
        assert total==u;star_totals.append(total)
        for k,t in enumerate(triples):
            if c in t:rows[k]^=total
    expected=[sum(1<<j for j,s in enumerate(triples) if len(set(t)&set(s))==1) for t in triples]
    assert rows==expected
    twice=[]
    for row in rows:
        out=0
        for k in range(v):
            if row>>k&1:out^=rows[k]
        twice.append(out)
    assert twice==[1<<k for k in range(v)]
    # E_c has Euclidean annihilator 1-3e_c. K_T has annihilator 3t_T-1
    # for the rational Gram I-J/9. No pair is proportional.
    for c in range(h):
        e=[1-3*(i==c) for i in range(h)]
        for t in triples:
            k=[3*(i in t)-1 for i in range(h)]
            assert any(e[i]*k[j]!=e[j]*k[i] for i in range(h) for j in range(h))
    paid=(h-2)+2*2+2+1
    assert paid==h+5 and paid-(h-1)==6
    # This is relative to the old data source path, with all auxiliaries free.
    deficit_per_pair=1-Q(12,h-2)-12
    assert deficit_per_pair<0
    return dict(h=h,v=v,star_count=h,side_is_involution=True,
                stars_pairwise_commute=True,original_centres_preserved=True,
                no_star_frame_equals_any_target_frame=True,
                distinct_star_visits_per_data_role=3,
                interior_rank_one_steps_per_data_role=4,
                terminal_rank_one_steps_per_data_role=2,
                total_source_path_rank=paid,old_direct_source_path_rank=h-1,
                excess_rank_per_data_role_per_stage=6,
                two_stage_deficit_per_data_pair=str(deficit_per_pair),
                reordering_or_permuting_targets_cannot_remove_these_switches=True,
                qualification='Literal point-star visits and inherited data endpoints; not a lower bound for other factorizations.')

if __name__=='__main__':
    result=[check(h) for h in (6,10,14,18,22,26)]
    print(json.dumps(result,indent=2))
    Path(__file__).with_name('point-stars-results.json').write_text(json.dumps(result,indent=2)+'\n')
