"""Read-only adapter from a pinned PR84 word to sparse cleanup plan inputs.

No upstream compiler is imported. Frame supports are all triples in the
specified common-core/envelope space, so tests using these supports are
conservative for the actual serialized rational frame.
"""
if not __debug__: raise RuntimeError('Run without -O')
from itertools import combinations
from pathlib import Path
import gzip,json


def bits(x):
    while x:
        bit=x&-x;yield bit.bit_length()-1;x^=bit


def model(source,h=23):
    source=Path(source)
    raw=gzip.decompress((source/f'certificates/indexed-cycle-word-{h}.json.gz').read_bytes())
    d=json.loads(raw);triples=list(combinations(range(h),3));v=len(triples);R=d['R']
    masks=[sum(1<<i for i in t) for t in triples]
    source_by_core={}
    for i,t in enumerate(triples):
        for c in (*[(z,) for z in t],*combinations(t,2),t):
            mask=sum(1<<z for z in c);source_by_core.setdefault(mask,[]).append(i)
    supports=[0];ranks=[0]
    for core,cover in d['frames']:
        sup=sum(1<<i for i in source_by_core[core] if masks[i]&~cover==0)
        assert sup
        supports.append(sup);ranks.append(1 if core==cover else cover.bit_count()-core.bit_count())
    hold=[[] for _ in range(R)]
    for s,before,after in d['events']:
        assert (hold[s][-1]-1 if hold[s] else -1)==before
        hold[s].append(after+1)
    out={};ret={}
    for s,g,c,t in d['outputs']:
        if len(t)==1:ret[s]=c
        else:out[s]=(c,tuple(t))
    C=dict(roles=R,hold=hold,out=out,ret=ret)
    G=dict(trip=triples,sup=supports,dn=ranks)
    mid=[(a,b,g+1) for a,b,g in d['ops']]
    sources={int(i):s for i,s in d['sources'].items()};signal=[0]*R
    for i,s in sources.items():signal[s]=1<<i
    for a,b,g in mid:signal[a]^=signal[b]
    candidates=[]
    for s in range(R):
        n=hold[s][-1]
        assert signal[s] and not signal[s]&~supports[n]
        if s in out:continue
        core,cover=d['frames'][n-1]
        for c in bits(core):candidates.append((s,c,tuple(bits(signal[s]))))
    return dict(h=h,G=G,C=C,mid=mid,sources=sources,triples=triples,
                v=v,R=R,signal=signal,candidates=candidates,serialized=d)
