"""Optimistic moment relaxations of the borrowed two-stage architecture.
Relaxed profiles are lower bounds, not proposed realizable circuits.
"""

if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')

from collections import Counter
from fractions import Fraction as Q
from math import comb,exp,log
from pathlib import Path
import json,sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'exploration/indexed-cycle-target-10'))
import independent_moment as cert


def profile(a,b,rho=Q(2),selected='flag',centre=True,data='ideal'):
    a,b=int(a),int(b);m=a*b;N=comb(a,3)*comb(b,3)
    assert a>=5 and b>=5 and rho>=0 and data in ('ideal','indexed')
    # Every data path is batched as favorably as its total edge rank permits.
    # Add the small blocks separately: they coincide for symmetric dimensions.
    rows=Counter({(a-1)*(b-1):2*N,1:N});rows[a-1]+=2*N;rows[b-1]+=2*N
    if data=='indexed':
        assert (a,b)==(23,25)
        rows=Counter({1:19*N,21:2*N,17:2*N,481:2*N})
        # Both source paths optimally match the superseded direct growth.
        for h in (a,b):rows[1]+=2*N;rows[h-2]+=2*N
    for h in (a,b):
        bank=rho*N
        if selected=='flag':rows[m-2*h]+=bank;rows[h]+=bank
        elif selected=='perfect':rows[m-h]+=bank
        else:raise ValueError(selected)
        rows[h]+=bank # optimistically put ALL internal mass into one h block
        if centre:rows[h-1]+=Q(N*h,comb(h,3))
    W=N*(2+2*rho);mass=sum(w*n for w,n in rows.items());D=W*m-mass
    assert D==N*(1-(Q(6,a-2)+Q(6,b-2) if centre else 0))
    return dict(a=a,b=b,m=m,N=N,W=W,rows=rows,deficit=D,rho=str(rho),selected=selected,centre=centre,data=data)


def root_float(p):
    if p['deficit']<=0:return 0.
    terms=[(float(Q(w*n,p['m']*p['W'])),log(p['m']/w)) for w,n in p['rows'].items()]
    lo,hi=0.,.1
    for _ in range(55):
        mid=(lo+hi)/2
        if sum(weight*exp(mid*z) for weight,z in terms)<1:lo=mid
        else:hi=mid
    return (lo+hi)/2

if __name__=='__main__':
    for rho,selected,centre in ((Q(2),'flag',True),(Q(2),'perfect',True),(Q(0),'flag',True),(Q(2),'flag',False),(Q(0),'flag',False)):
        best=max((root_float(p:=profile(a,b,rho,selected,centre)),a,b) for a in range(9,101) for b in range(a,101))
        print(rho,selected,centre,best,flush=True)
    for rho in (Q(0),Q(1),Q(2),Q(3)):
        print('indexed23x25',rho,root_float(profile(23,25,rho,data='indexed')),flush=True)
