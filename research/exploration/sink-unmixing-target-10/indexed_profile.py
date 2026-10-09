"""Exact first-stage sink-unmixing profile on the pinned PR84 producer.

The original producer and its baseline histogram are credited to PR84. This
file changes only complete-word exits, adds every reference path, and leaves
the second tensor stage unchanged (or uses our separately audited borrowed
wrapper there). No claimed second-stage retirement merge is included.
"""
if not __debug__: raise RuntimeError('Run without -O')
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from math import floor, log
from pathlib import Path
import gzip,json,sys,time
import frame_profiles as fp
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'indexed-cycle-target-10'))
import independent_moment as arithmetic


def change(counter, profile, factor=1):
    for r,n in profile.items():counter[r]+=factor*n


def mass(counter):return sum(r*n for r,n in counter.items())


def float_root(rows,m,W):
    def f(a):return sum(n*(r/m)**(1-a) for r,n in rows.items())/W
    lo=0.;hi=.01
    assert f(lo)<1<f(hi)
    for _ in range(70):
        mid=(lo+hi)/2
        if f(mid)<1:lo=mid
        else:hi=mid
    return lo


def certify(rows,m,W):
    rows=Counter({r:n for r,n in rows.items() if n})
    assert all(0<r<m and n>0 for r,n in rows.items())
    approximate=float_root(rows,m,W)
    saving=Q(floor(approximate*10**12)-1,10**12)
    p=dict(rows=rows,m=m,W=W)
    _,hi=arithmetic.moment(p,saving)
    targetlo,_=arithmetic.moment(p,Q(1,1023))
    assert hi<1<targetlo
    return dict(W=W,total_rank=mass(rows),rank_deficit=m*W-mass(rows),
                bit_saving=str(saving),approximate_root=approximate,
                exact_moment_gap=str(1-hi),target_10_excess_lower=str(targetlo-1),
                child_multiplicities={str(r):n for r,n in sorted(rows.items())})


def run(source):
    started=time.monotonic();commit=arithmetic.check_source(source)
    a,b,m,n=23,25,575,4073300;repeat=2300;v=1771
    d=json.loads(gzip.decompress((source/'certificates/indexed-cycle-word-23.json.gz').read_bytes()))
    baseline=arithmetic.profile();last=[None]*d['R']
    for s,before,after in d['events']:
        assert last[s]==(None if before==-1 else before)
        last[s]=after
    outputs={row[0]:row for row in d['outputs']}
    side={s for s,g,c,t in d['outputs'] if len(t)==3}
    eligible=Counter()
    for s,g in enumerate(last):
        if s in side:continue
        core,cover=d['frames'][g]
        assert core and core&~cover==0
        c=(core&-core).bit_length()-1
        eligible[core,cover,c]+=1
    assert sum(eligible.values())==d['R']-3*v
    ext=Counter({m-2*a:1,a:1})
    complements={c:fp.complement_profile(a,0,c)[0] for c in range(a)}
    assert all(p==Counter({1:1,a-2:1,m-2*a+2:1}) for p in complements.values())
    I=fp.identity(a);delta=Counter();local_records=[]
    for index,((core,cover,c),multiplicity) in enumerate(sorted(eligible.items())):
        A=fp.frame(a,core,cover);U=fp.frame(a,1<<c,(1<<a)-1)
        rank=1 if core==cover else cover.bit_count()-core.bit_count()
        old=fp.profile(fp.sub(I,A));new=fp.profile(fp.sub(U,A))
        assert mass(old)==a-rank and mass(new)==a-1-rank
        change(delta,old,-multiplicity);change(delta,new,multiplicity)
        change(delta,ext,-multiplicity);change(delta,complements[c],multiplicity)
        local_records.append([core,cover,c,multiplicity,dict(old),dict(new)])
        if index and index%2000==0: print('exact local exit pairs',index,flush=True)
    assert mass(delta)==0
    triples=list(combinations(range(a),3));ref_profiles={};detours=Counter()
    for S in triples:
        P=fp.frame(a,sum(1<<i for i in S),sum(1<<i for i in S))
        direct=fp.profile(fp.sub(I,P));assert direct==Counter({1:1,a-2:1})
        for c in S:
            U=fp.frame(a,1<<c,(1<<a)-1)
            split=fp.profile(fp.sub(U,P));end=fp.profile(fp.sub(I,U))
            assert mass(split)==a-2 and end==Counter({1:1})
            ref=Counter({1:1});change(ref,split);change(ref,end);change(ref,ext)
            assert mass(ref)==m
            ref_profiles[c,S]=ref
            if c==min(S):
                change(detours,split);change(detours,end);change(detours,direct,-1)
    assert mass(detours)==0
    cases=[]
    # Original same-producer baseline and separately checked borrowed baseline.
    original=certify(baseline['rows'],m,baseline['W']);original['name']='original_unborrowed'
    cases.append(original)
    borrowed=Counter(baseline['rows'])
    for h in (a,b):change(borrowed,Counter([m-2*h,h,1,1,h-2]),-n)
    q=certify(borrowed,m,baseline['W']-2*n);q['name']='old_borrowed_both_stages';cases.append(q)
    for refs_per_source in (3,2):
        local_delta=Counter(delta)
        for S in triples:
            for c in S:
                if refs_per_source==2 and c==min(S):continue
                change(local_delta,ref_profiles[c,S])
        if refs_per_source==2:change(local_delta,detours)
        assert mass(local_delta)==refs_per_source*v*m
        for second in ('original','borrowed'):
            rows=Counter(baseline['rows']);change(rows,local_delta,repeat)
            W=baseline['W']+refs_per_source*n
            if second=='borrowed':
                change(rows,Counter([m-2*b,b,1,1,b-2]),-n);W-=n
            out=certify(rows,m,W)
            assert out['rank_deficit']==baseline['deficit']
            out.update(name=f'first_stage_{refs_per_source}v_references_second_{second}',
                       first_stage_references=refs_per_source*v,
                       first_stage_retired_roles=sum(eligible.values()))
            cases.append(out)
    return dict(status='Exact Q profiles for a proposed first-stage scalar-word replacement; global interfaces and full literal integration remain separate.',
                source_commit=commit,dimensions=[a,b],m=m,N=n,
                baseline_producer_roles=d['R'],first_stage_early_retirement_roles=sum(eligible.values()),
                no_anchor_fallback_roles=0,exact_distinct_exit_pairs=len(eligible),
                merged_complement_profiles={str(c):dict(p) for c,p in complements.items()},
                first_stage_retirement_delta_per_invocation=dict(sorted(delta.items())),
                primary_X_detour_delta_per_invocation=dict(sorted(detours.items())),
                local_exit_profile_records=local_records,cases=cases,
                elapsed_seconds=time.monotonic()-started)


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__);parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,default=HERE/'indexed-profile-results.json')
    args=parser.parse_args();result=run(args.source)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([{k:c[k] for k in ('name','W','bit_saving','approximate_root')} for c in result['cases']],indent=2))
