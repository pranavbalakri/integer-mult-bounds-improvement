"""Independent defining-space and direct full-tensor checks over Q."""
if not __debug__:raise RuntimeError('Run without -O')
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import json
import frame_profiles as fp


def transpose(A):return [list(row) for row in zip(*A)]
def multiply(A,B):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]


def inverse(A):
    n=len(A);M=[list(row)+[Q(i==j) for j in range(n)] for i,row in enumerate(A)]
    for j in range(n):
        k=next(k for k in range(j,n) if M[k][j]);M[j],M[k]=M[k],M[j]
        p=M[j][j];M[j]=[x/p for x in M[j]]
        for i in range(n):
            if i!=j:
                p=M[i][j];M[i]=[x-p*y for x,y in zip(M[i],M[j])]
    return [row[n:] for row in M]


def gram_projector(h,triples):
    # Select an independent rational source basis without trusting envelope rank.
    vectors=[[Q(i in t) for i in range(h)] for t in triples]
    basis=[]
    for v in vectors:
        if len(fp.pivots_exact(basis+[v]))>len(basis):basis.append(v)
    L=transpose(basis);G=[[Q(i==j)-Q(1,9) for j in range(h)] for i in range(h)]
    T=[[Q(i==j)+1 for j in range(h)] for i in range(h)]
    Ti=[[Q(i==j)-Q(1,h+1) for j in range(h)] for i in range(h)]
    gram=multiply(multiply(transpose(L),G),L)
    raw=multiply(multiply(multiply(L,inverse(gram)),transpose(L)),G)
    return multiply(multiply(T,raw),Ti)


def full_complement(a,axis,X):
    h=(a,a+2)[axis];other=(a,a+2)[1-axis]
    # Genuine nonzero dense line, independently chosen from a triple projector.
    Qother=fp.frame(other,7,7)
    coords=fp.reversed_coordinates(a);m=len(coords)
    Y=[[X[x[axis]][y[axis]]*Qother[x[1-axis]][y[1-axis]] for y in coords] for x in coords]
    return fp.sub(fp.identity(m),Y)


def run():
    count=0;direct=[];examples=[]
    for a in (5,6):
        for axis in (0,1):
            h=(a,a+2)[axis];full=(1<<h)-1
            labels=[(7,7),(1,full),(1<<h-1,full),(3,full),(1,full^(1<<h-1))]
            for core,cover in labels:
                ts=[t for t in combinations(range(h),3)
                    if all(cover>>i&1 for i in t) and all(i in t for i in range(h) if core>>i&1)]
                actual=gram_projector(h,ts);formula=fp.frame(h,core,cover)
                assert actual==[list(row) for row in formula]
                assert multiply(actual,actual)==actual
                large=full_complement(a,axis,actual)
                direct_piv=fp.pivots_exact(large)
                predicted_profile,predicted_piv=fp.complement_from_local(a,axis,actual)
                assert direct_piv==predicted_piv
                assert fp.profile_pivots(direct_piv)==predicted_profile
                count+=1
                direct.append(dict(dimensions=[a,a+2],axis=axis,core=core,cover=cover,
                                   source_span_rank=len(fp.pivots_exact(actual)),
                                   merged_profile=dict(predicted_profile)))
                if a==5 and axis==0:
                    examples.append(dict(core=core,cover=cover,
                        local_projector=[[str(x) for x in row] for row in actual],
                        global_merged_exit=[[str(x) for x in row] for row in large]))
    # Exact expected local-to-global contraction, plus independent U_c route.
    for axis in (0,1):
        h=(23,25)[axis]
        for c in range(h):
            X=fp.frame(h,1<<c,(1<<h)-1)
            assert fp.complement_profile(23,axis,c)==fp.complement_from_local(23,axis,X)
    # Matrix sensitivity control: omitting one normal rank-one part changes rank.
    I=fp.identity(5);U=fp.frame(5,1,31)
    assert fp.profile(full_complement(5,0,I))!=fp.profile(full_complement(5,0,U))
    return dict(status='PASS exact defining Gram spaces and full rational tensor elimination',
                direct_full_matrix_cases=count,all_48_target_complements_checked=True,
                omitted_normal_direction_detected=True,cases=direct,explicit_examples=examples)


if __name__=='__main__':
    result=run();Path(__file__).with_name('frame-check-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['direct_full_matrix_cases'])
