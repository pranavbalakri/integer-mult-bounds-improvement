"""Product-centre fitting identity and explicitly unrealized moment diagnostics.
This does not certify a new common-frame circuit or multiplication saving.
"""
if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'borrowed-ceiling-target-10'))
import ideal_moment as old
from requirements import root_bracket


def fitting():
    rows=[]
    for s in range(4):
        for t in range(4):
            central=(s%2)*(t%2)
            gram=Q((s-1)*(t-1),4)
            side=central^int(s==3 and t==3)
            first=int(s==1 and t==3)
            second=int(s==3 and t==1)
            both=int(s==1 and t==1)
            assert side==first^second^both
            assert not side or gram==0
            if s==t==3:assert central==gram==1
            elif central:assert gram==0
            rows.append(dict(first_intersection=s,second_intersection=t,
                             binary_centre_coefficient=central,rational_gram=str(gram),
                             side_coefficient=side))
    return rows


def profile(a,b):
    p=old.profile(a,b,Q(2),centre=False)
    d=(a-1)*(b-1);L=a*b*d
    p['rows'][d]+=a*b;p['deficit']-=L
    assert Q(L,p['N'])==Q(36,(a-2)*(b-2))
    oldL=b*(b-1)*comb(a,3)+a*(a-1)*comb(b,3)
    target=Q(1,1023);lo,hi=old.cert.moment(p,target)
    result=dict(a=a,b=b,m=a*b,N=p['N'],auxiliary_ratio=2,
                hypothetical_product_centre_count=a*b,
                hypothetical_product_centre_rank=d,
                hypothetical_total_centre_rank=L,old_total_centre_rank=oldL,
                hypothetical_centre_rank_fraction=str(Q(L,p['N'])),
                old_centre_rank_fraction=str(Q(oldL,p['N'])),
                relaxed_root_interval=root_bracket(p),
                target_moment_interval=[str(lo),str(hi)],
                missing_step='No routing with this centre loss and these inherited data/auxiliary profiles is supplied.')
    if (a,b) in ((12,12),(10,14)):
        _,upper=old.cert.moment(p,Q(1,800));assert upper<1
        result['hypothetical_profile_passes_1_over_800']=True
        result['exact_hypothetical_moment_gap']=str(1-upper)
    return result


def direct_product_star_path(a,b):
    # Every data row belongs to a 3 by 3 grid of product stars.
    # Two distinct product-star frames intersect in the tensor product of
    # their local intersections. A monotone through-the-sum switch costs
    # twice the rank difference; the full frame cannot cost less.
    rank=(a-1)*(b-1)
    same_first_cost=2*(a-1)
    same_second_cost=2*(b-1)
    diagonal_cost=2*(a+b-3)
    cheap,expensive=sorted((same_first_cost,same_second_cost))
    # A Hamiltonian snake attains the lower bound: at least two edges must
    # change the expensive coordinate, while every other edge costs >=cheap.
    assert diagonal_cost>=expensive
    switch_cost=6*cheap+2*expensive
    initial=rank-1;final=a*b-rank
    assert initial+final==a*b-1
    return dict(a=a,b=b,product_star_visits_per_data_role=9,
                minimum_switch_count=8,minimum_extra_rank_per_data_role=switch_cost,
                initial_plus_final_rank=a*b-1,
                direct_in_place_route_cannot_fit_available_credit=True,
                qualification='Literal visit to all nine product-star frames. Other producer factorizations are not covered.')

if __name__=='__main__':
    result=dict(fitting_identity=fitting(),
                hypothetical_moments=[profile(a,b) for a,b in ((10,14),(12,12),(23,25))],
                direct_frame_path_bounds=[direct_product_star_path(a,b) for a,b in ((6,6),(10,14),(12,12),(23,25))],
                status='Useful scalar/fitting ingredients and exact hypothetical budgets; no new certified exponent.')
    Path(__file__).with_name('tensor-relaxation-results.json').write_text(json.dumps(result,indent=2)+'\n')
    for p in result['hypothetical_moments']:
        print(p['a'],p['b'],p['relaxed_root_interval'])
