"""Dirty-reference word with an input-independent auxiliary map undone at sink.
This checks the complete scalar word and common-point cleanup incidences.
It does NOT certify a new fixed-basis child histogram or exponent.
"""
if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')
from pathlib import Path
from collections import Counter
import json,sys
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]
sys.path.insert(0,str(RESEARCH/'independent/two-stage-bit'))
import rtgm


def bits(row):
    while row:
        bit=row&-row;yield bit.bit_length()-1;row^=bit


def producer(G):
    C=rtgm.compile_(G);users={n:[] for n in G['active']}
    for n in sorted(G['active']):
        if G['args'][n]:
            for pos,k in enumerate(G['args'][n]):users[k].append(('gate',n,pos))
    for t,n in sorted(G['outputs'].items()):users[n].append(('side',t))
    for c,n in sorted(G['retained'].items()):users[n].append(('centre',c))
    nextslot=0;edge={};middle=[];sources={}
    for n in sorted(G['active']):
        if G['args'][n]:
            a,b=edge[n,0],edge[n,1];middle.append((a,b,n));pivot=a
        else:
            pivot=nextslot;nextslot+=1;sources[n-1]=pivot
        for k,u in enumerate(users[n]):
            s=pivot
            if k:s=nextslot;nextslot+=1;middle.append((s,pivot,n))
            if u[0]=='gate':edge[u[1],u[2]]=s
    assert nextslot==C['roles']
    return C,middle,sources


def check(h,two_references=True):
    G=rtgm.build(h,retain='search');C,mid,sources=producer(G)
    triples=G['trip'];v=len(triples);R=C['roles'];primary=[min(t) for t in triples]
    refs=[(c,s) for s,t in enumerate(triples) for c in t if not two_references or c!=primary[s]]
    bref={key:2*v+R+i for i,key in enumerate(refs)};B=len(refs);size=2*v+R+B
    assert B==(2 if two_references else 3)*v
    middle=[(2*v+a,2*v+b) for a,b,n in mid]
    inject=[(2*v+s,i) for i,s in sorted(sources.items())]
    injectB=[(bref[c,s],s) for c,s in refs]
    Jcentre=[(v+i,2*v+s) for s,c in sorted(C['ret'].items()) for i,t in enumerate(triples) if c in t]
    lookup={t:i for i,t in enumerate(triples)}
    Jside=[(v+lookup[t],2*v+s) for s,(c,t) in sorted(C['out'].items())]
    J=Jcentre+Jside
    initial=[1<<i for i in range(size)]
    dirt=initial[:]
    for a,b in middle:dirt[a]^=dirt[b]
    D=[0]*v
    for a,b in J:D[a-v]^=dirt[b]
    rmask=((1<<R)-1)<<(2*v)
    assert all(row&~rmask==0 for row in D)
    early=[(v+i,s) for i,row in enumerate(D) for s in bits(row)]
    # Compute every producer row's clean-source coefficient, independently.
    signal=[0]*size
    for i in range(v):signal[i]=1<<i
    for a,b in inject+middle:signal[a]^=signal[b]
    cleanout=[0]*v
    for a,b in J:cleanout[a-v]^=signal[b]
    assert cleanout==[1<<i for i in range(v)]
    early_clean=[];late_clean=[];KB=[];anchors={};fallback=[]
    for s in range(R):
        finalnode=C['hold'][s][-1]
        label=G['sup'][finalnode]
        assert signal[2*v+s]&~label==0
        common=set(range(h))
        for i in bits(label):common.intersection_update(triples[i])
        if common:c=min(common)
        else:c=None;fallback.append(s)
        anchors[s]=c
        group=early_clean if s not in C['out'] and c is not None else late_clean
        for i in bits(signal[2*v+s]):
            # A no-anchor fallback is kept at D1 and uses any valid reference.
            chosen=c if c is not None else primary[i]
            assert chosen in triples[i]
            control=i if two_references and chosen==primary[i] else bref[chosen,i]
            group.append((2*v+s,control))
            if control>=2*v+R:KB.append((2*v+s,control))
    clearB=injectB
    sink_unmix=KB+list(reversed(middle))
    word=early+injectB+inject+middle+J+early_clean+late_clean+clearB+sink_unmix
    for transpose in (False,True):
        state=initial[:]
        for a,b in word:
            if transpose:a,b=b,a
            state[a]^=state[b]
        expected=initial[:]
        if transpose:
            for i in range(v):expected[i]^=initial[v+i]
        else:
            for i in range(v):expected[v+i]^=initial[i]
        assert state==expected
    # Check the exact intermediate block before free sink unmixing.
    state=initial[:]
    prefix=word[:-len(sink_unmix)]
    for a,b in prefix:state[a]^=state[b]
    xmask=(1<<v)-1;ymask=xmask<<v
    assert all(row&(xmask|ymask)==0 for row in state[2*v:])
    assert state[2*v+R:]==initial[2*v+R:]
    expectedR=dirt[2*v:2*v+R]
    for a,b in KB:expectedR[a-2*v]^=initial[b]
    assert state[2*v:2*v+R]==expectedR
    # The reference dirt must be unmixed; deleting that whole block fails.
    bad=initial[:]
    for a,b in prefix+list(reversed(middle)):bad[a]^=bad[b]
    assert bad[2*v:]!=initial[2*v:]
    # Frame obligations checked combinatorially: each early label is in U_c;
    # each reference (c,S) has P_S subset U_c and stays in this ONE U_c.
    assert all(c in triples[s] for c,s in refs)
    for s,c in anchors.items():
        if s not in C['out'] and c is not None:
            assert all(c in triples[i] for i in bits(G['sup'][C['hold'][s][-1]]))
    early_roles={a-2*v for a,b in early_clean}
    late_roles={a-2*v for a,b in late_clean}
    assert not early_roles&late_roles
    assert early_roles|late_roles==set(range(R))
    return dict(h=h,v=v,producer_roles=R,reference_roles=B,total_auxiliary_roles=R+B,
                complete_public_basis_dimension=size,word_cnot_count=len(word),
                early_dirty_correction_cnot_count=len(early),
                reference_dirt_unmix_cnot_count=len(KB),
                common_point_cleanup_roles=len(early_roles),
                full_stage_cleanup_roles=len(late_roles),
                no_common_anchor_fallback_roles=fallback,
                forward_all_dirty_identity=True,same_chronology_transposed_identity=True,
                intermediate_auxiliary_map_independent_of_data=True,
                omitting_reference_dirt_unmix_detected=True,
                reference_mode='2v plus read-only primary X' if two_references else '3v dirty references',
                cleanup_rank_target=h-1,
                largest_stage_one_exit_rank_without_profiling=h*h-h+1,
                status='Scalar word and common-point cleanup checked. Exact frame matrices, child partitions, moment, and global interface remain to audit.')

if __name__=='__main__':
    dims=list(map(int,sys.argv[1:])) or [5,6,7,8]
    result=[check(h,mode) for h in dims for mode in (False,True)]
    Path(__file__).with_name('scalar-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
