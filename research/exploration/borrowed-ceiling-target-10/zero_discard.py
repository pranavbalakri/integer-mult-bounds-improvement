"""Honest optimistic rank budget for zero-input, discardable producer auxiliaries.

Every auxiliary can begin at its first required frame and disappear after its
last use; no entrance, full-frame exit, dirty correction, or cleanup is charged.
Only actual intervening growth and the necessary side/centre transfer remain.
The data roles retain the inherited endpoints. A negative rank deficit rejects
this relaxation before any global zero-auxiliary transposition theorem is needed.
No certified producer, profile, or external source is changed.
"""

if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')

from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
from math import comb
import gzip,json,sys
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]
sys.path.insert(0,str(RESEARCH/'independent/two-stage-bit'))
import rtgm
import ideal_moment as model


def scalar_positive(G,C):
    """Independent complete clean-input scalar replay, linked to compiler slots."""
    users={n:[] for n in G['active']}
    for n in sorted(G['active']):
        if G['args'][n]:
            for pos,k in enumerate(G['args'][n]):users[k].append(('gate',n,pos))
    for target,n in sorted(G['outputs'].items()):users[n].append(('side',target))
    for c,n in sorted(G['retained'].items()):users[n].append(('centre',c))
    values=[];edge={};output=[0]*len(G['trip']);lookup={t:i for i,t in enumerate(G['trip'])}
    for n in sorted(G['active']):
        if G['args'][n]:
            a,b=edge[n,0],edge[n,1];values[a]^=values[b];pivot=a
        else:
            pivot=len(values);values.append(1<<(n-1))
        assert values[pivot]==G['sup'][n]
        for k,u in enumerate(users[n]):
            s=pivot
            if k:s=len(values);values.append(values[pivot])
            if u[0]=='gate':edge[u[1],u[2]]=s
            elif u[0]=='side':output[lookup[u[1][1]]]^=values[s]
            else:
                for i,t in enumerate(G['trip']):
                    if u[1] in t:output[i]^=values[s]
    assert len(values)==C['roles']
    assert output==[1<<i for i in range(len(output))]


def ordinary(h, optimized=False):
    if optimized:
        sys.path.insert(0,str(RESEARCH/'independent/bit-improvement'))
        import reassociate
        old,G,C=reassociate.histogram(h)
    else:
        G=rtgm.build(h,retain='search');C=rtgm.compile_(G)
    scalar_positive(G,C)
    local=Counter();initial=final=0
    for s,hold in enumerate(C['hold']):
        ranks=[G['dn'][n] for n in hold]
        assert ranks
        if s in C['out']:ranks.append(h-1)
        if s in C['ret']:assert ranks[-1]==h-1
        assert all(a<=b for a,b in zip(ranks,ranks[1:]))
        initial+=ranks[0];final+=h-ranks[-1]
        local.update(b-a for a,b in zip(ranks,ranks[1:]) if b>a)
    q=sum(r*n for r,n in local.items())
    assert q==h*C['roles']-initial-final
    return dict(kind='paired/reassociated' if optimized else 'ordinary retained total',h=h,
                v=comb(h,3),roles=C['roles'],initial_rank_removed=initial,
                final_rank_removed=final,remaining_local_rank=q,
                local_histogram={str(r):n for r,n in sorted(local.items())},
                clean_source_to_combined_scatter_is_identity=True)


def indexed(root,h):
    model.cert.check_source(root)
    d=json.loads(gzip.decompress((root/f'certificates/indexed-cycle-word-{h}.json.gz').read_bytes()))
    triples=list(combinations(range(h),3));v=len(triples)
    rank=[1 if c==u else u.bit_count()-c.bit_count() for c,u in d['frames']]
    values=[0]*d['R'];first=[None]*d['R'];last=[None]*d['R'];local=Counter()
    def arrive(slot,fr):
        if first[slot] is None:first[slot]=fr
        if last[slot] is not None:
            delta=rank[fr]-rank[last[slot]];assert delta>=0
            if delta:local[delta]+=1
        last[slot]=fr
    lookup={tuple(f):i for i,f in enumerate(d['frames'])}
    for i,slot in d['sources'].items():
        i=int(i);mask=sum(1<<k for k in triples[i]);arrive(slot,lookup[mask,mask]);values[slot]=1<<i
    for a,b,fr in d['ops']:
        arrive(a,fr);arrive(b,fr);values[a]^=values[b]
    output=[0]*v;terminal=set()
    for slot,fr,c,t in d['outputs']:
        assert slot not in terminal;terminal.add(slot);arrive(slot,fr)
        if len(t)==3:
            assert rank[fr]==h-3
            # Two paid ranks reach t_T-perp. The old third singleton was cleanup.
            local[1]+=2;output[triples.index(tuple(t))]^=values[slot]
        else:
            assert len(t)==1 and rank[fr]==h-1
            for i,q in enumerate(triples):
                if c in q:output[i]^=values[slot]
    assert output==[1<<i for i in range(v)]
    assert all(f is not None for f in first+last)
    q=sum(r*n for r,n in local.items())
    initial=sum(rank[f] for f in first)
    final=sum(h-(h-1 if s in terminal else rank[f]) for s,f in enumerate(last))
    assert q==h*d['R']-initial-final
    return dict(kind='pinned indexed PR84',h=h,v=v,roles=d['R'],initial_rank_removed=initial,
                final_rank_removed=final,remaining_local_rank=q,
                local_histogram={str(r):n for r,n in sorted(local.items())},
                clean_source_to_combined_scatter_is_identity=True)


def global_budget(A,B,data='ideal'):
    a,b=A['h'],B['h'];p=model.profile(a,b,Q(0),data=data)
    for local,other in ((A,B),(B,A)):
        h=local['h']
        for r,n in local['local_histogram'].items():
            r=int(r)
            blocks=[1]*r if 2*r<=h else [1]*(h-r)+[2*r-h]
            assert sum(blocks)==r
            for w in blocks:p['rows'][w]+=other['v']*n
    mass=sum(r*n for r,n in p['rows'].items())
    D=p['m']*p['W']-mass
    assert D==p['N']*(1-Q(6,a-2)-Q(6,b-2)-Q(A['remaining_local_rank'],A['v'])-Q(B['remaining_local_rank'],B['v']))
    assert D<0
    return dict(a=a,b=b,data_profile=data,W=int(p['W']),N=p['N'],m=p['m'],total_rank=int(mass),
                deficit=int(D),deficit_per_data_pair=str(Q(D,p['N'])),
                normalized_moment_at_zero=str(Q(mass,p['m']*p['W'])),
                any_positive_saving_rejected=True,
                child_histogram={str(r):int(n) for r,n in sorted(p['rows'].items())})

if __name__=='__main__':
    ap=ArgumentParser();ap.add_argument('--source',type=Path);args=ap.parse_args()
    locals=[ordinary(h) for h in (5,6,7,8,10,12,16,18,23)]
    locals.append(ordinary(23,optimized=True))
    budgets=[global_budget(p,p) for p in locals]
    if args.source:
        records=[indexed(args.source,h) for h in (23,25)]
        locals.extend(records);budgets.append(global_budget(*records,data='indexed'))
    result=dict(qualification='Optimistically removes all auxiliary entrances and cleanup. The inherited data endpoints and centre transfer costs remain; negative rank deficits alone reject these producers.',
                producers=locals,global_budgets=budgets)
    (HERE/'zero-discard-results.json').write_text(json.dumps(result,indent=2)+'\n')
    for p in locals:print(p['kind'],p['h'],'q/v=',Q(p['remaining_local_rank'],p['v']))
    for p in budgets:print(p['a'],p['b'],p['data_profile'],'D/N=',p['deficit_per_data_pair'])
