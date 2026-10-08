"""Conservative chronological histogram for the proposed in-place complex word.

Retains explicit source-line and low-frame vertices, so no merged-edge gain is
claimed. Independently reconstructs the new histogram then checks the exact
old-histogram deletion identity. This remains an exploratory new interface.
"""
if not __debug__:
    raise RuntimeError("This exact checker requires assertions; do not run Python with -O.")

from collections import Counter
from pathlib import Path
import json
import sys

RESEARCH=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(RESEARCH/'independent/complex-twostage'))
from producer import NStar3,compile_roles
from hist import build,label_dims
from cert import cert


def profile(h):
    c=NStar3(h);kk=compile_roles(c);dims=label_dims(c)
    v=len(c.triples);N=v*v;m=h*h;R=kk['size']
    sources=set(kk['src'].values())
    cur={s:int(s in sources) for s in range(R)}
    local=Counter()
    def grow(s,d):
        assert cur[s]<=d
        if cur[s]<d:local[d-cur[s]]+=1
        cur[s]=d
    for n,ins,outs in kk['gates']:
        for s in set(ins+outs):grow(s,dims[n])
    centre_loss=0
    for name,s in kk['rout'].items():
        if cur[s]:local[cur[s]]+=1
        centre_loss+=cur[s]
    for k,s in kk['pout'].items():grow(s,h-1)
    for s in range(R):grow(s,h)
    assert centre_loss==h+(h-1)**2
    assert sum(r*n for r,n in local.items())==h*R-v+centre_loss
    H=Counter({r:n*2*v for r,n in local.items()})
    H[m-h]+=2*v*(R-v)  # all Z explicitly visit D0 before local growth
    H[h-1]+=2*N        # Y's target-complement local growth in both stages
    H[(h-1)**2]+=2*N   # explicit interstage entry vertices on X and Y
    H[1]+=N            # unchanged completed endpoint correction
    W=2*v*R
    raw=sum(r*n for r,n in H.items())
    assert raw==m*W-N+2*v*centre_loss
    old=build(h,copied=True)
    deletion=Counter(old['hist'])
    for r in (m-h,h-1,1):deletion[r]-=2*N
    deletion=Counter({r:n for r,n in deletion.items() if n})
    assert all(n>0 for n in deletion.values())
    assert H==deletion
    assert W==old['W']-2*N
    saving,gap=cert(dict(H),W,m)
    return dict(h=h,m=m,N=N,old_W=old['W'],W=W,R=R,remaining_auxiliary_roles_per_invocation=R-v,
                centre_loss_per_invocation=centre_loss,rank_sum=raw,rank_deficit=m*W-raw,
                maxrank=max(H),hist=dict(sorted(H.items())),saving=str(saving),moment_gap=str(gap),
                chronological_histogram_matches_conservative_deletion=True,
                scope='New local complex profile; inherited copy, recursive width and global assembly interfaces remain assumed')


if __name__=='__main__':
    print(json.dumps([profile(h) for h in map(int,sys.argv[1:] or ['16','18'])],indent=2))
