"""Independent direct reconstruction of the borrowed PR84 rectangle profile.
Uses the pinned local frame profiles, not the old global histogram subtraction.
"""

if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from math import comb
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'indexed-cycle-target-10'
sys.path.insert(0,str(SOURCE))
import independent_moment as cert
from borrowed_indexed import profile as comparison


def check():
    a,b=23,25;m=a*b;N=comb(a,3)*comb(b,3)
    # Unchanged rectangular data entrance and endpoint correction.
    rows=Counter({1:19*N,21:2*N,17:2*N,481:2*N})
    W=2*N;centre_loss=0;removed={}
    manifest=json.loads((SOURCE/'inputs/source.json').read_text())
    for h in (a,b):
        name=f'indexed-cycle-profiles-{h}.json';raw=(SOURCE/'inputs'/name).read_bytes()
        assert sha256(raw).hexdigest()==manifest['files']['certificates/'+name]
        local=json.loads(raw);v=comb(h,3);repeat=N//v
        Rnew=local['R']-v;bank=repeat*Rnew;W+=bank
        blocks=Counter({r:n for r,n in enumerate(local['blocks']) if r and n})
        # Each borrowed source already occupies its initial rank-one frame.
        blocks[1]-=v;assert blocks[1]>=0
        assert sum(r*n for r,n in blocks.items())==h*Rnew+(h-1)*v+local['loss']
        for r,n in blocks.items():rows[r]+=repeat*n
        # All surviving aux roles start at D0, even zero columns of D.
        rows[h]+=bank;rows[m-2*h]+=bank
        # Only the target's old rank-(h-1) edge remains; source now uses M's path.
        rows[1]+=N;rows[h-2]+=N
        centre_loss+=repeat*local['loss']
        removed[h]=dict(borrowed_per_invocation=v,invocations=repeat,
                        removed_global_roles=N,removed_rank=N*m,
                        old_direct_edge=[h-2,1],old_selected_edge=[m-2*h,h],
                        old_leaf_entrance=[1])
    rows=Counter({r:n for r,n in rows.items() if n})
    mass=sum(r*n for r,n in rows.items());deficit=W*m-mass
    assert (W,mass,deficit)==(124319508,71481870200,1846900)
    assert deficit==N-centre_loss
    other=comparison()
    assert rows==Counter({int(r):n for r,n in other['child_multiplicities'].items()})
    assert (W,mass)==(other['W'],other['total_rank'])
    saving=Q(5488,10**8)
    lo,hi=cert.moment(dict(W=W,m=m,rows=rows),saving)
    assert hi<1
    target_lo,_=cert.moment(dict(W=W,m=m,rows=rows),Q(1,1023))
    assert target_lo>1
    result=dict(W=W,m=m,N=N,total_rank=mass,rank_deficit=deficit,
                direct_reconstruction_matches_subtraction=True,
                per_axis_accounting=removed,bit_saving=str(saving),
                exact_moment_gap=str(1-hi),target_10_bit_moment_excess_lower=str(target_lo-1),
                scope='Independent finite-profile arithmetic; assumes the new lifted frame word and inherited global interfaces.')
    return result

if __name__=='__main__':
    result=check();print(json.dumps(result,indent=2))
    (HERE/'borrowed-rect-independent-audit.json').write_text(json.dumps(result,indent=2)+'\n')
