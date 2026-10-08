"""Exact target budgets and root brackets for explicitly optimistic profiles.
None of these relaxed profiles is asserted to be a realizable circuit.
"""

if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')

from fractions import Fraction as Q
from math import floor,comb
from pathlib import Path
import json
import ideal_moment as model
cert=model.cert
TARGET=Q(1,1023)


def root_bracket(p):
    grid=10**12;k=floor(model.root_float(p)*grid)
    while cert.moment(p,Q(k,grid))[1]>=1:k-=1
    while cert.moment(p,Q(k+1,grid))[0]<=1:k+=1
    assert cert.moment(p,Q(k,grid))[1]<1<cert.moment(p,Q(k+1,grid))[0]
    return dict(lower=str(Q(k,grid)),upper=str(Q(k+1,grid)))


def max_rank_one_cost(h):
    p=model.profile(h,h,Q(2),centre=False)
    v=comb(h,3)
    def moment(k):
        q=dict(p);q['rows']=p['rows'].copy();q['rows'][1]+=2*v*k
        return cert.moment(q,TARGET)
    lo,hi=0,h*(h-1)
    assert moment(lo)[1]<1<moment(hi)[0]
    while hi-lo>1:
        mid=(lo+hi)//2
        if moment(mid)[1]<1:lo=mid
        else:assert moment(mid)[0]>1;hi=mid
    return dict(h=h,v=v,extra_aux_roles_per_invocation=2*v,
                maximum_rank_one_charges_per_invocation=lo,
                next_integer_budget_rejected=hi,
                current_retained_centre_rank=h*(h-1),
                accepted_target_moment_gap=str(1-moment(lo)[1]),
                rejected_target_moment_excess=str(moment(hi)[0]-1))


def check_rho_boundary(a,b,data,lo,hi,centre):
    p=model.profile(a,b,lo,centre=centre,data=data)
    q=model.profile(a,b,hi,centre=centre,data=data)
    assert cert.moment(p,TARGET)[1]<1<cert.moment(q,TARGET)[0]
    return dict(a=a,b=b,data=data,retained_centre_cost_present=centre,
                maximum_uniform_aux_ratio_interval=[str(lo),str(hi)])


def audit():
    cases=[('fixed indexed data; no auxiliary roles',model.profile(23,25,Q(0),data='indexed')),
           ('fixed indexed data; two auxiliary roles per target',model.profile(23,25,Q(2),data='indexed')),
           ('ideal data; flag selected edges; finite-scan maximizer',model.profile(24,25,Q(2))),
           ('ideal data and selected edges; finite-scan maximizer',model.profile(24,25,Q(2),'perfect')),
           ('ideal data; no centres; two auxiliary roles per target',model.profile(23,25,Q(2),centre=False))]
    profiles=[]
    for name,p in cases:
        lo,hi=cert.moment(p,TARGET)
        assert lo>1
        profiles.append(dict(name=name,a=p['a'],b=p['b'],root_interval=root_bracket(p),
                             target_moment_excess_lower=str(lo-1)))
    p=model.profile(24,25,Q(0))
    assert cert.moment(p,Q(115,100000))[1]<1
    profiles.append(dict(name='ideal data; no auxiliary roles; retained centres',a=24,b=25,
                         root_interval=root_bracket(p),target_numerically_possible=True,
                         qualification='An unrealizable zero-role relaxation, not a producer.'))
    # A hypothetical endpoint-free architecture has enough room: this is a requirement,
    # NOT an implementation or an asserted removable charge.
    p=model.profile(13,13,Q(2));p['rows'][1]-=p['N'];p['deficit']+=p['N']
    assert cert.moment(p,Q(16,10000))[1]<1
    profiles.append(dict(name='hypothetical deletion of the paid endpoint copy',a=13,b=13,
                         root_interval=root_bracket(p),target_numerically_possible=True,
                         qualification='No circuit removing that endpoint charge is supplied.'))
    ratios=[check_rho_boundary(23,25,'ideal',Q(1740684,10**7),Q(1740685,10**7),True),
            check_rho_boundary(23,25,'indexed',Q(3976826,10**7),Q(3976827,10**7),False)]
    budgets=[max_rank_one_cost(h) for h in (6,7,8,10,12,14,16,18,20)]
    return dict(target_bit_saving=str(TARGET),relaxed_profiles=profiles,
                required_aux_ratios=ratios,rank_one_reset_budgets=budgets,
                status='Quantitative necessary changes and optimistic budgets only; no target-achieving circuit.')

if __name__=='__main__':
    result=audit()
    print(json.dumps(result,indent=2))
    Path(__file__).with_name('requirements-results.json').write_text(json.dumps(result,indent=2)+'\n')
