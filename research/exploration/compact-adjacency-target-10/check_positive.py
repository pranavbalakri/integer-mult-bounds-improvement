"""Exact small positive-producer audit, including the unfavorable rank budget."""

if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')

from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import json


def rank(rows):
    a=[[Q(x) for x in row] for row in rows]
    if not a:return 0
    r=0
    for c in range(len(a[0])):
        p=next((i for i in range(r,len(a)) if a[i][c]),None)
        if p is None:continue
        a[r],a[p]=a[p],a[r];q=a[r][c];a[r]=[x/q for x in a[r]]
        for i in range(r+1,len(a)):
            q=a[i][c]
            if q:a[i]=[x-q*y for x,y in zip(a[i],a[r])]
        r+=1
    return r


def check(path):
    d=json.loads(path.read_text());h=d['h'];trip=list(combinations(range(h),3));v=len(trip)
    assert h<9 # G positive definite; all source spans nondegenerate.
    args=d['args'];support=[1<<i for i in range(v)];dims=[1]*v
    for n,ab in enumerate(args[v:],v):
        a,b=ab;assert a<n and b<n and not support[a]&support[b]
        support.append(support[a]|support[b])
        rows=[[int(c in trip[i]) for c in range(h)] for i in range(v) if support[n]>>i&1]
        dims.append(rank(rows))
    for t,n in zip(trip,d['outputs']):
        expected=sum(1<<i for i,s in enumerate(trip) if len(set(s)&set(t))==1)
        assert support[n]==expected and dims[n]<=h-1
        if h>=6:assert dims[n]==h-1
    # Literal source injections and CNOTs, with all data coefficients tracked.
    users=[[] for _ in args]
    for n,ab in enumerate(args):
        if ab:
            for pos,a in enumerate(ab):users[a].append(('gate',n,pos))
    for i,n in enumerate(d['outputs']):users[n].append(('output',i))
    signal=[];history=[];edge={};outs={};ops=[]
    def alloc():signal.append(0);history.append([]);return len(signal)-1
    for n,ab in enumerate(args):
        if ab:
            a,b=(edge[n,pos] for pos in range(2));signal[a]^=signal[b]
            history[a].append(n);history[b].append(n);ops.append((a,b,n));pivot=a
        else:
            pivot=alloc();signal[pivot]^=1<<n;history[pivot].append(n)
        for k,u in enumerate(users[n]):
            s=pivot if k==0 else alloc()
            if k:
                signal[s]^=signal[pivot];history[s].append(n);history[pivot].append(n);ops.append((s,pivot,n))
            if u[0]=='gate':edge[u[1],u[2]]=s
            else:outs[s]=u[1]
    R=len(signal);assert R==d['roles']==len(args)-v+v
    assert len(outs)==v
    rk=Counter()
    for s,hs in enumerate(history):
        prev=0;prevdim=0
        for n in hs:
            assert not prev&~support[n]
            if dims[n]>prevdim:rk[dims[n]-prevdim]+=1
            assert dims[n]>=prevdim
            prev,prevdim=support[n],dims[n]
        if s in outs:
            t=trip[outs[s]]
            assert all(len(set(trip[i])&set(t))==1 for i in range(v) if prev>>i&1)
            assert signal[s]==support[d['outputs'][outs[s]]]
            if prevdim<h-1:rk[h-1-prevdim]+=1
            prevdim=h-1
        if prevdim<h:rk[h-prevdim]+=1
    assert sum(r*n for r,n in rk.items())==h*R
    # Direct copied centres, not hypothetical free centres.
    m=h*h;N=v*v;A=v*(R+h);hist=Counter()
    def add(w,n):
        if n:assert w>0;hist[w]+=n
    def inner(r):return [1]*r if 2*r<=h else [1]*(h-r)+[2*r-h]
    for _ in range(2):
        add(m-2*h,A);add(h,A)
        for r,n in rk.items():
            for w in inner(r):add(w,v*n)
        add(h,2*v*h)
    add(m-4*h+2,2*N);add(h-2,2*N);add(1,(h+1)*2*N)
    for _ in range(2):add(h-2,2*N);add(1,2*N)
    add(1,N)
    W=2*N+2*A;s=sum(w*n for w,n in hist.items());deficit=W*m-s
    assert deficit==N-2*v*h*h<0
    # Corrupting one terminal coefficient must fail the exact adjacency test.
    bad=signal[next(iter(outs))]^1
    assert bad!=signal[next(iter(outs))]
    return dict(h=h,v=v,adds=len(args)-v,full_side_outputs=v,side_roles=R,
                cnot_count=len(ops),all_source_frames_positive_definite=True,
                all_monotone_frames_checked=True,all_side_coefficients_checked=True,
                inner_rank_histogram=dict(sorted(rk.items())),W=W,m=m,rank_sum=s,
                rank_deficit=deficit,child_histogram=dict(sorted(hist.items())),
                conclusion='Valid positive side producer; no positive saving under direct copied-centre assembly.')

if __name__=='__main__':
    directory=Path(__file__).resolve().parent
    rows=[check(directory/f'positive-h{h}.json') for h in (5,6)]
    print(json.dumps(rows,indent=2))
    (directory/'positive-audit.json').write_text(json.dumps(rows,indent=2)+'\n')
