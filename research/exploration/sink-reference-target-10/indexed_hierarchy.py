"""Profile a literal zero-reference hierarchy on the pinned PR84 producer.

The original word, baseline and reversed rectangular basis are credited to
the adjacent indexed-cycle package. Only first-stage paths change here.
All source detours are paid; no rank-only single-child assumption is used.
"""
if not __debug__:
    raise RuntimeError('Run without -O')

from collections import Counter
from functools import lru_cache
from pathlib import Path
import argparse
import json
import sys
import time

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'sink-unmixing-target-10'))
import indexed_adapter
import indexed_profile as ip
import frame_profiles as fp
import reference_plans as rp
import hierarchy


def run(source,orders=('large-first',),full=True):
    started=time.monotonic()
    commit=ip.arithmetic.check_source(source)
    M=indexed_adapter.model(source,23)
    P=rp.choose(M,[min(t) for t in M['triples']])
    G=M['G'];C=M['C'];data=M['serialized']
    a,b,m,N,repeat=23,25,575,4073300,2300
    baseline=ip.arithmetic.profile()
    assert M['h']==a and M['v']==N//repeat
    exterior=Counter({m-2*a:1,a:1})
    fullkey=('full',)
    def frame(key):
        return fp.identity(a) if key==fullkey else fp.frame(a,*key)
    @lru_cache(None)
    def difference(before,after):
        result=fp.profile(fp.sub(frame(after),frame(before)))
        assert all(w>0 and n>0 for w,n in result.items())
        return result
    complements={}
    def complement(key):
        if key not in complements:
            complements[key]=fp.complement_from_local(a,0,frame(key))[0]
            if len(complements)%100==0:
                print(f'exact merged-frame profiles: {len(complements)}',flush=True)
        return complements[key]
    def nodekey(n):return tuple(data['frames'][n-1])
    def rank(key):
        if key==fullkey:return a
        core,cover=key
        return 1 if core==cover else cover.bit_count()-core.bit_count()
    def star(c):return (1<<c,(1<<a)-1)
    cases=[]
    for order in orders:
        Q=hierarchy.select(M,P,order)
        audit=hierarchy.audit(M,Q,full)
        local_delta=Counter();retirement_delta=Counter();source_delta=Counter()
        actual_exit_pairs=Counter()
        for s,c in Q['anchors'].items():
            A=nodekey(C['hold'][s][-1])
            endpoint=A if s in Q['local_nodes'] else star(c)
            assert s not in Q['local_nodes'] or Q['local_nodes'][s]==C['hold'][s][-1]
            actual_exit_pairs[A,endpoint]+=1
        for (A,endpoint),multiplicity in actual_exit_pairs.items():
            old=difference(A,fullkey)
            new=difference(A,endpoint)
            exit=complement(endpoint)
            assert ip.mass(old)==a-rank(A)
            assert ip.mass(new)==rank(endpoint)-rank(A)
            assert ip.mass(exit)==m-rank(endpoint)
            ip.change(retirement_delta,old,-multiplicity)
            ip.change(retirement_delta,exterior,-multiplicity)
            ip.change(retirement_delta,new,multiplicity)
            ip.change(retirement_delta,exit,multiplicity)
        for i in Q['used_primary']:
            mask=sum(1<<j for j in M['triples'][i]);first=(mask,mask)
            old=difference(first,fullkey)
            assert old==Counter({1:1,a-2:1})
            ip.change(source_delta,old,-1)
            previous=first
            for n in Q['source_chains'][i]:
                current=nodekey(n)
                ip.change(source_delta,difference(previous,current))
                previous=current
            current=star(Q['primary'][i])
            ip.change(source_delta,difference(previous,current))
            ip.change(source_delta,difference(current,fullkey))
        assert ip.mass(retirement_delta)==0==ip.mass(source_delta)
        ip.change(local_delta,retirement_delta);ip.change(local_delta,source_delta)
        for second in ('original','borrowed'):
            rows=Counter(baseline['rows'])
            ip.change(rows,local_delta,repeat)
            W=baseline['W']
            if second=='borrowed':
                ip.change(rows,Counter([m-2*b,b,1,1,b-2]),-N)
                W-=N
            certificate=ip.certify(rows,m,W)
            assert certificate['rank_deficit']==baseline['deficit']
            certificate.update(selection_order=order,second_stage=second,
                               first_stage_reference_roles=0,
                               first_stage_early_roles=len(Q['anchors']),
                               first_stage_original_frame_cleanups=len(Q['local_nodes']))
            cases.append(certificate)
        audit.update(selection_order=order,
                     retirement_delta_per_first_stage=dict(sorted(retirement_delta.items())),
                     source_data_delta_per_first_stage=dict(sorted(source_delta.items())),
                     total_delta_per_first_stage=dict(sorted(local_delta.items())))
        yield dict(type='plan',record=audit)
    yield dict(type='summary',record=dict(
        status='Exact scalar replay and rational child profile; global integration is not certified.',
        source_commit=commit,dimensions=[a,b],m=m,N=N,
        source='Pinned PR84 indexed-cycle word and baseline; see ../indexed-cycle-target-10/UPSTREAM-NOTICE.',
        exact_merged_frames=len(complements),exact_local_differences=difference.cache_info().currsize,
        every_new_path_rank_accounted=True,cases=cases,elapsed_seconds=time.monotonic()-started))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--order',action='append',choices=['small-first','large-first','multiplicity'])
    parser.add_argument('--skip-full-dirty',action='store_true')
    parser.add_argument('--output',type=Path,default=HERE/'indexed-hierarchy-results.json')
    args=parser.parse_args()
    result=list(run(args.source,args.order or ['large-first'],not args.skip_full_dirty))
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    for item in result:
        if item['type']=='summary':
            print(json.dumps(item['record'],indent=2))
