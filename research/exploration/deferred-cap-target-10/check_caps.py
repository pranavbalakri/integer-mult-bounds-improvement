"""Exact componentwise caps on the pinned deferred entrance frames.

Original frame data, scalar word, and common U,V basis are credited to the
pinned Swapnil round-seven source. This changes only selected sigma frames.
"""
if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')

from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace
import argparse
import hashlib
import json
import random
import sys
import time

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'deferred-truncation-target-10'))
from check_truncation import PINS


def run(source,cap):
    started=time.monotonic()
    for name,digest in PINS.items():
        assert hashlib.sha256((source/name).read_bytes()).hexdigest()==digest,name
    sys.path.insert(0,str(source/'independent/deferred-readout'))
    import deferred as dr
    import check_frames as cf
    from linalg import Q31,null_exact,rank_mod,rightmost_pivots,predicted
    assert cap in (19,20)
    S=dr.Schedule(*dr.load(23));h=S.h;cov=S.adjoint()
    high=[s for s in S.readout if S.f[s]>cap]
    parent={s:s for s in high};target_owner={}
    def find(s):
        while parent[s]!=s:
            parent[s]=parent[parent[s]];s=parent[s]
        return s
    for s in high:
        for t in cov[s]:
            if t in target_owner:parent[find(s)]=find(target_owner[t])
            else:target_owner[t]=s
    components=defaultdict(list)
    for s in high:components[find(s)].append(s)
    by_target=defaultdict(list)
    for s in S.readout:
        for t in cov[s]:by_target[t].append(s)
    @lru_cache(None)
    def annihilator(B):return null_exact(B,h)
    def inside(A,B):
        null=annihilator(tuple(map(tuple,B)))
        return all(sum(x*y for x,y in zip(row,z))==0 for row in A for z in null)
    def nondeg(B):
        sums=[sum(row) for row in B]
        gram=[[9*sum(x*y for x,y in zip(a,b))-sa*sb
               for b,sb in zip(B,sums)] for a,sa in zip(B,sums)]
        return rank_mod(gram,Q31)==len(B)
    caps={};component_records=[];inclusions=0;extensions=0
    rng=random.Random(2026100800+cap)
    for _,ss in sorted(components.items()):
        targets=sorted({t for s in ss for t in cov[s]})
        proper=sorted({s for t in targets for s in by_target[t] if S.f[s]<=cap})
        container=min(ss,key=lambda s:(S.f[s],s));H=S.sigma[container]
        assert nondeg(H)
        # One least-dimensional old high frame is the intersection in these
        # instances; check the actual inclusion, not only its dimension.
        for s in ss:
            assert inside(H,S.sigma[s]);inclusions+=1
        base=max(proper,key=lambda s:(S.f[s],-s)) if proper else None
        B=[row[:] for row in S.sigma[base]] if base is not None else []
        assert len(B)<=cap and nondeg(B) and inside(B,H)
        for s in proper:
            assert inside(S.sigma[s],B);inclusions+=1
        extra=[];rejections=0
        while len(B)<cap:
            for attempt in range(100):
                weights=[rng.randrange(-32,33) for row in H]
                row=[sum(a*r[j] for a,r in zip(weights,H)) for j in range(h)]
                if nondeg(B+[row]):break
            else:raise AssertionError(('extension search exhausted',targets,len(B)))
            assert inside([row],H)
            B.append(row);extra.append(row);rejections+=attempt
        assert len(B)==cap and nondeg(B)
        for t in targets:
            normal=[9*int(j in S.trip[t])-3 for j in range(h)]
            assert all(sum(x*y for x,y in zip(row,normal))==0 for row in B)
        for s in ss:caps[s]=B
        extensions+=bool(extra)
        component_records.append(dict(slots=ss,targets=targets,base_slot=base,
          original_base_dimension=0 if base is None else S.f[base],
          container_slot=container,extra_integer_rows=extra,
          rejected_extension_candidates=rejections))
    print(f'cap{cap}: {len(components)} exact component caps; {extensions} extended components',flush=True)

    bases=[];indices={}
    def frame(B):
        k=tuple(map(tuple,B))
        if k not in indices:indices[k]=len(bases);bases.append(B)
        return ('sigma',indices[k])
    oldkeys={s:frame(B) for s,B in S.sigma.items()}
    capkeys={s:frame(B) for s,B in caps.items()}
    def key(s):return capkeys.get(s,oldkeys[s])
    targetkeys={}
    for t in target_owner:
        normal=[9*int(j in S.trip[t])-3 for j in range(h)]
        targetkeys[t]=frame(null_exact([normal],h))
    extra={('sigma',i):B for i,B in enumerate(bases)}
    mq=cf.ModQ7(SimpleNamespace(h=h,pr=None),{},extra)
    def dim(k):return h if k==('F',) else 0 if k==('0',) else len(bases[k[1]])
    edges={};raw=Counter()
    def edge(A,B,kind):
        r=dim(B)-dim(A);assert r>=0
        if not r:return
        raw[kind]+=1;edges[A,B]=r
    for s in high:
        # Explicitly retain old sigma as an identity vertex. Its onward
        # original edge and common-basis certificate remain unchanged.
        edge(key(s),oldkeys[s],'inserted_first_edge')
        edge(key(s),('F',),'auxiliary_corner')
        old_first=S.dim(S.start_key(s))-S.f[s]
        assert 0<=old_first<=h-cap
        assert (S.f[s]-cap)+old_first==S.dim(S.start_key(s))-cap
        assert 2*max(S.f[s]-cap,old_first)<=h
        assert 2*(S.dim(S.start_key(s))-cap)<=h
    for ring in ('Z','F2'):
        for t in sorted(target_owner):
            slots=[s for s in by_target[t] if ring=='Z' or cov[s][t]&1]
            previous=('0',);changed_before=False
            for s in slots:
                current=key(s);changed=s in caps
                if changed or changed_before:edge(previous,current,'Y_'+ring)
                previous=current;changed_before=changed
            if changed_before:edge(previous,targetkeys[t],'Y_'+ring)
    used={k for e in edges for k in e}
    assert all(mq.nondeg(k) for k in used) and mq.dimfail==0
    want=predicted(h);checked=Counter()
    for A,B in sorted(edges,key=str):
        r=edges[A,B];qa,qb=mq.conjugated(A),mq.conjugated(B)
        for side in range(2):
            M=[[(qb[side][i*h+j]-qa[side][i*h+j])%Q31 for j in range(h)] for i in range(h)]
            assert rightmost_pivots(M,Q31)==want[r],(A,B,r,side)
            checked[r]+=1
    assert mq.dimfail==0
    print(f'cap{cap}: both fixed U,V NE profiles pass for {len(edges)} changed edges',flush=True)
    return dict(source_repository='https://github.com/Swapnil-jain/integer-mult-kappa',
      source_commit='741e7aa078392553815df7926ee17ac5e25a8c38',source_sha256=PINS,
      h=h,cap=cap,changed_slots=len(high),components=len(components),
      old_high_slots_with_multiple_integer_targets=sum(len(cov[s])>1 for s in high),
      predecessor_dimensions=dict(sorted(Counter(r['original_base_dimension'] for r in component_records).items())),
      extended_components=extensions,exact_inclusion_comparisons=inclusions,
      all_caps_integer_defined=True,all_caps_nondegenerate=True,
      all_unchanged_predecessors_contained=True,all_caps_inside_every_original_high_frame=True,
      common_axis_basis='check_lifted.point_UV(23,1)',modulus=Q31,
      explicit_old_sigma_identity_vertices=True,
      singleton_first_edge_rank_histogram_matches_unsplit_dimension_count=True,
      distinct_changed_edges=len(edges),changed_edges_by_kind=dict(raw),
      NE_checks_by_rank=dict(sorted(checked.items())),NE_failures=0,
      entrance_dimension_histogram=dict(sorted(Counter(min(f,cap) for f in S.f).items())),
      largest_auxiliary_block=h*h-2*(h-cap),total_rank_mass_unchanged=True,
      extension_seed=2026100800+cap,component_witnesses=component_records,
      scalar_word='Unchanged pinned signed word; only readout/auxiliary entrance frames are capped.',
      inherited_frame_obligations='Original sigma-to-start inclusions, unchanged scalar word, and unchanged frame edges are inherited from the pinned source checks.',
      full_multiplication_bound=False,elapsed_seconds=round(time.monotonic()-started,3))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--cap',type=int,choices=[19,20],required=True)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=run(args.source.resolve(),args.cap)
    output=args.output or HERE/f'cap{args.cap}-results.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('component_witnesses','source_sha256')},indent=2))
