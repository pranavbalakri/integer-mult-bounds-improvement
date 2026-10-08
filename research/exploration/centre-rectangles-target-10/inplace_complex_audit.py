"""Exact dirty-basis audit of an in-place complex producer word.

This audits scalar identities and local binary frame histories, not the complete
lifted histogram or recursive multiplication theorem. The second orientation is
same-chronology inverse-transpose, with identical scalar-gate frame incidences.
Written with OpenAI Codex assistance.
"""
if not __debug__:
    raise RuntimeError("This exact checker requires assertions; do not run Python with -O.")

from pathlib import Path
import json
import sys

RESEARCH=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(RESEARCH/'independent/complex-twostage'))
from producer import NStar3,compile_roles
from frames import Checker,perp_in,ok_res,echelon,red,dot


def add(row,other,num=1,den=1):
    for k,v in list(other.items()):
        z=num*v
        assert z%den==0
        z=row.get(k,0)+z//den
        if z:row[k]=z
        elif k in row:del row[k]


def run(h=8,frames_only=False):
    c=NStar3(h);kk=compile_roles(c);v=len(c.triples);R=kk['size']
    tid={t:i for i,t in enumerate(c.triples)}
    source={sl:tid[t] for t,sl in kk['src'].items()}
    assert len(source)==v
    other=[s for s in range(R) if s not in source]
    physical=source|{s:2*v+i for i,s in enumerate(other)}
    count=2*v+len(other)
    checker=Checker(c);ambient=[1<<i for i in range(h)]
    history={s:([sum(1<<j for j in c.triples[physical[s]])] if s in source else []) for s in range(R)}
    middle=[]
    for n,ins,outs in kk['gates']:
        L=checker.label(n)
        for s in set(ins+outs):
            assert checker.edge(history[s],L)
            history[s]=L
        pivot=physical[ins[0]]
        middle.extend((pivot,physical[s],1,1) for s in ins[1:])
        middle.extend((physical[s],pivot,1,1) for s in outs[1:])
    Jcentre=[]
    for i,t in enumerate(c.triples):
        last=h-1 in t
        Jcentre.append((v+i,physical[kk['rout'][('*',)]],5-h if last else 2,2))
        points=[j for j in range(h-1) if j not in t] if last else t
        Jcentre.extend((v+i,physical[kk['rout'][('E',j)]],1 if last else -1,2) for j in points)
    Jside=[]
    for k,slot in kk['pout'].items():
        t,n,coef=c.pieces[k]
        target_perp=perp_in(ambient,[sum(1<<j for j in t)])
        assert checker.edge(history[slot],target_perp)
        history[slot]=target_perp
        Jside.append((v+tid[t],physical[slot],coef,2))
    for name,slot in kk['rout'].items():
        assert ok_res(history[slot]) # inverse phase on a retained copy
    assert all(checker.edge(L,ambient) for L in history.values())
    # All tensor entry/outer residuals are Kronecker products of these
    # checked one-factor spaces, or their orthogonal direct sums.
    perps=[]
    for t in c.triples:
        line=[sum(1<<j for j in t)]
        A=perp_in(ambient,line)
        assert len(A)==h-1 and ok_res(A)
        assert all(dot(a,line[0])==0 for a in A)
        perps.append(A)
    if frames_only:
        return dict(h=h,v=v,old_auxiliary_roles=R,new_auxiliary_roles=R-v,
                    common_frame_middle_gates=len(kk['gates']),
                    elementary_producer_updates=len(middle),
                    exact_distinct_local_frame_transitions=len(checker.cache),
                    centre_copy_frames_checked=len(kk['rout']),side_incidence_frames_checked=len(Jside),
                    final_producer_role_growths_checked=R,triple_perpendicular_spaces_checked=len(perps),
                    tensor_data_entry_family='a_perp tensor b_perp; all factors nondegenerate and nonalternating',
                    tensor_outer_family='a_perp tensor F and F tensor b_perp; all factors nondegenerate and nonalternating',
                    tensor_local_family='a tensor residual or residual tensor b; odd norm lines preserve Gram form',
                    explicit_source_line_and_low_frame_vertices=True,
                    all_local_frame_checks_pass=True,
                    scope='Actual local role histories plus exact tensor-factor residual proof; global recursion not certified')

    # Compute the exact doubled read matrix 2JM, including every dirty column.
    midmap=[{i:2} for i in range(count)]
    for a,b,n,d in middle:add(midmap[a],midmap[b],n,d)
    image=[{} for _ in range(v)]
    for a,b,n,d in Jcentre+Jside:add(image[a-v],midmap[b],n,d)
    for i,row in enumerate(image):
        assert {j:z for j,z in row.items() if j<2*v}=={i:2}
    early=[(v+i,j,-z,2) for i,row in enumerate(image) for j,z in row.items() if j>=2*v]
    word=early+middle+Jcentre+Jside+[(a,b,-n,d) for a,b,n,d in reversed(middle)]
    for inverse_transpose in (False,True):
        values=[{i:2} for i in range(count)]
        for a,b,n,d in word:
            if inverse_transpose:a,b,n=b,a,-n
            add(values[a],values[b],n,d)
        for i in range(count):
            wanted={i:2}
            if not inverse_transpose and v<=i<2*v:wanted[i-v]=2
            if inverse_transpose and i<v:wanted[v+i]=-2
            assert values[i]==wanted,(inverse_transpose,i,values[i],wanted)
    assert early
    # Omitting one actual cancellation leaves a nonzero dirty coefficient.
    a,b,n,d=early[0]
    assert n!=0 and b>=2*v
    return dict(h=h,v=v,old_auxiliary_roles=R,new_auxiliary_roles=R-v,
                total_scalar_roles=count,borrowed_source_roles=v,
                elementary_producer_updates=len(middle),early_cancellation_updates=len(early),
                centre_reads=len(Jcentre),side_reads=len(Jside),word_length=len(word),
                full_dirty_basis_forward=True,full_dirty_basis_same_order_inverse_transpose=True,
                local_phase_histories_valid=True,omitted_cancellation_leaves_nonzero_dirty_coefficient=True,
                scope='Exact scalar and local frame audit; lifted histogram and global recursion not certified')


if __name__=='__main__':
    frames_only='--frames-only' in sys.argv
    hs=[int(a) for a in sys.argv[1:] if a!='--frames-only'] or [8]
    print(json.dumps([run(h,frames_only) for h in hs],indent=2))
