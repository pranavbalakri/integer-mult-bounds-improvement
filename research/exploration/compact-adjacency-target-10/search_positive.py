"""Experimental positive shared-sum producers for full intersection-one outputs.
Only scalar DAG counts; no multiplication theorem is certified by this search.
"""
from collections import Counter
from itertools import combinations
from random import Random
import argparse, json


def search(h, seed):
    rng = Random(seed)
    trip = list(combinations(range(h),3)); v=len(trip)
    masks = [sum(1<<i for i in t) for t in trip]
    rows = [{j for j,s in enumerate(masks) if (s&t).bit_count()==1} for t in masks]
    args = [None]*v; supports=[1<<i for i in range(v)]
    while True:
        counts=Counter(p for row in rows for p in combinations(sorted(row),2))
        if not counts: break
        best=max(counts.values())
        if best<2: break
        cand=[p for p,c in counts.items() if c==best]
        # Size tie break varies by seed while retaining exact sharing count.
        if seed%3:
            score=lambda p:(supports[p[0]]|supports[p[1]]).bit_count()
            opt=(max if seed%3==1 else min)(map(score,cand))
            cand=[p for p in cand if score(p)==opt]
        a,b=rng.choice(cand); k=len(args)
        assert not supports[a]&supports[b]
        args.append((a,b)); supports.append(supports[a]|supports[b])
        for row in rows:
            if a in row and b in row:
                row.difference_update((a,b));row.add(k)
    outputs=[]
    for row in rows:
        queue=sorted(row,key=lambda n:supports[n].bit_count())
        if not queue: outputs.append(None);continue
        x=queue[0]
        for y in queue[1:]:
            assert not supports[x]&supports[y]
            args.append((x,y));supports.append(supports[x]|supports[y]);x=len(args)-1
        outputs.append(x)
    adds=len(args)-v
    return dict(h=h,v=v,seed=seed,adds=adds,roles=adds+sum(x is not None for x in outputs),args=args,outputs=outputs)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--h',type=int,nargs='+',default=[5,6,7,8]);p.add_argument('--seeds',type=int,default=10);p.add_argument('--output');a=p.parse_args()
    for h in a.h:
        best=min((search(h,s) for s in range(a.seeds)),key=lambda d:d['roles'])
        print({k:v for k,v in best.items() if k not in ('args','outputs')},flush=True)
        if a.output:
            from pathlib import Path
            Path(a.output.format(h=h)).write_text(json.dumps(best,indent=2)+'\n')
