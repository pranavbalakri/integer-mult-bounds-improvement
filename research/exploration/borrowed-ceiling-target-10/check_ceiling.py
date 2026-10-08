"""Exact all-dimension ceilings for optimistic borrowed-source profiles.
The assumptions and profile relaxations are documented in README.md.
"""

if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')

from fractions import Fraction as Q
from pathlib import Path
import json


def log_integer_lower(n):
    k=n.bit_length()-1;x=Q(n,1<<k);u=(x-1)/(x+1)
    ln2=2*sum((Q(1,3)**(2*j+1)/Q(2*j+1) for j in range(24)),Q())
    lo=k*ln2+2*sum((u**(2*j+1)/Q(2*j+1) for j in range(24)),Q())
    return Q((lo*10**12).__floor__(),10**12)


def large_block_entropy(m,width):
    # width*log(m/width), lower bounded by a positive -log(1-t) series.
    t=Q(m-width,m);power=t;answer=Q()
    for j in range(1,9):answer+=power/j;power*=t
    return width*answer


def entropy_lower(a,b,logs,selected,rho=Q(2)):
    m=a*b;d=(a-1)*(b-1)
    # The data paths are completely lumped; even the source-line edges are ideal.
    answer=2*large_block_entropy(m,d)
    answer+=(2*a-1)*logs[b]+(2*b-1)*logs[a]+4-Q(2,a)-Q(2,b)
    # The stated uniform ratio rho. Grant all local work one full h block.
    for h,other in ((a,b),(b,a)):
        if selected=='flag':
            answer+=rho*(large_block_entropy(m,m-2*h)+2*h*logs[other])
        elif selected=='perfect':
            answer+=rho*(large_block_entropy(m,m-h)+h*logs[other])
        else:raise ValueError(selected)
    # Centre-copy entropy is positive; discarding it weakens the bound safely.
    return answer


def verify(selected,denominator,tail):
    # In the tail the displayed log terms alone suffice, since log(a)>2.
    coefficient=6 if selected=='flag' else 4
    assert 2*(coefficient*tail-1)>denominator
    logs={i:log_integer_lower(i) for i in range(9,tail)}
    assert logs[9]>2
    worst=None;count=0
    for a in range(9,tail):
        for b in range(a,tail):
            deficit=1-Q(6,a-2)-Q(6,b-2)
            if deficit<=0:continue
            E=entropy_lower(a,b,logs,selected)
            slack=E-denominator*deficit
            assert slack>0,(a,b,selected,slack)
            if worst is None or slack<worst[0]:worst=(slack,a,b)
            count+=1
    return dict(selected_profile=selected,minimum_aux_per_v=2,
                necessary_saving_bound=f'a_bit < 1/{denominator}',finite_pairs=count,
                finite_minimum_slack=str(worst[0]),finite_minimum_slack_pair=list(worst[1:]),
                tail=f'max(a,b)>={tail}',target_unattainable=Q(1,denominator)<Q(1,1023))

def verify_target_ratio(selected,rho,tail=214):
    coefficient=2+(2 if selected=='flag' else 1)*rho
    assert 2*(coefficient*tail-1)>1023
    logs={i:log_integer_lower(i) for i in range(9,tail)}
    assert logs[9]>2
    worst=None;count=0
    for a in range(9,tail):
        for b in range(a,tail):
            deficit=1-Q(6,a-2)-Q(6,b-2)
            if deficit<=0:continue
            slack=entropy_lower(a,b,logs,selected,rho)-1023*deficit
            assert slack>0,(a,b,selected,rho,slack)
            if worst is None or slack<worst[0]:worst=(slack,a,b)
            count+=1
    return dict(selected_profile=selected,target='1/1023',
                target_rejected_if_both_axis_aux_ratios_at_least=str(rho),
                finite_pairs=count,finite_minimum_slack=str(worst[0]),
                finite_minimum_slack_pair=list(worst[1:]),tail=f'max(a,b)>={tail}',
                scope='Retained centre copies and fully batched inherited data; not arbitrary circuits.')


if __name__=='__main__':
    rows=[verify('flag',2500,209),verify('perfect',1700,213)]
    print(json.dumps(rows,indent=2))
    Path(__file__).with_name('ceiling-results.json').write_text(json.dumps(rows,indent=2)+'\n')
    ratios=[verify_target_ratio('flag',Q(1,5)),verify_target_ratio('perfect',Q(2,5))]
    print(json.dumps(ratios,indent=2))
    Path(__file__).with_name('role-ratio-results.json').write_text(json.dumps(ratios,indent=2)+'\n')
