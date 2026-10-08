"""Exact scalar checks for coupled binary point centres; no frame claim.
Every scalar is a packed row on all public inputs, including shared dirt.
"""
if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')
from itertools import combinations
from pathlib import Path
import json


def rank(rows):
    piv={}
    for row in rows:
        while row:
            k=row.bit_length()-1
            if k in piv:row^=piv[k]
            else:piv[k]=row;break
    return len(piv)


def bits(n):
    while n:
        bit=n&-n;yield bit.bit_length()-1;n^=bit


def motif(h):
    triples=list(combinations(range(h),3))
    B=[sum(1<<i for i,t in enumerate(triples) if c in t) for c in range(h)]
    C=[0]*len(triples)
    for i,t in enumerate(triples):
        for c in t:C[i]^=B[c]
    assert rank(B)==h and rank(C)==h
    return triples,B,C


def product_row(a,b,width):
    out=0
    for i in bits(a):out^=b<<(i*width)
    return out


def linear(rows,values):
    out=[]
    for row in rows:
        v=0
        for i in bits(row):v^=values[i]
        out.append(v)
    return out


def check(a,b):
    ta,Ba,Ca=motif(a);tb,Bb,Cb=motif(b);va,vb=len(ta),len(tb);N=va*vb
    C1=[product_row(Ca[i],1<<j,vb) for i in range(va) for j in range(vb)]
    C2=[Cb[j]<<(i*vb) for i in range(va) for j in range(vb)]
    C12=[product_row(Ca[i],Cb[j],vb) for i in range(va) for j in range(vb)]
    r1,r2,rm=rank(C1),rank(C2),rank(C12)
    assert (r1,r2,rm)==(a*vb,b*va,a*b)
    E=[C1[i]^C2[i]^C12[i]^(C2[i]<<N) for i in range(N)]+C1
    assert rank(E)==r1+r2
    record=dict(a=a,b=b,N=N,first_centre_rank=r1,second_centre_rank=r2,
                mixed_centre_rank=rm,missing_two_corrections_rank=rank(E),
                mixed_totals_do_not_span_the_complete_correction=True)
    if a%4!=2 or b%4!=2:return record
    assert all((u&v).bit_count()%2==0 for B in (Ba,Bb) for u in B for v in B)
    V=[product_row(u,v,vb) for u in Ba for v in Bb]
    H=[sum(1<<(c*b+d) for c in t for d in u) for t in ta for u in tb]
    assert len(V)==rm and len(H)==N
    assert linear(H,V)==C12
    assert linear(V,H)==[0]*rm
    initial=[1<<i for i in range(2*N+rm)]
    X,Y,Z=initial[:N],initial[N:2*N],initial[2*N:]
    # L(C1), U(C2), L(C1), U(C2), in chronological order.
    for bank,other,C in ((Y,X,C1),(X,Y,C2),(Y,X,C1),(X,Y,C2)):
        add=linear(C,other)
        for i,value in enumerate(add):bank[i]^=value
    expectedX=[initial[i]^v for i,v in enumerate(linear(C12,initial[:N]))]
    expectedY=[initial[N+i]^v for i,v in enumerate(linear(C12,initial[N:2*N]))]
    assert X==expectedX and Y==expectedY and Z==initial[2*N:]
    # Realize the same mixed transvection with ONE shared arbitrary Z bank.
    # The four literal block additions restore Z because V H=0.
    X,Y,Z=initial[:N],initial[N:2*N],initial[2*N:]
    for bank in (X,Y):
        for dest,src,M in ((bank,Z,H),(Z,bank,V),(bank,Z,H),(Z,bank,V)):
            add=linear(M,src)
            for i,value in enumerate(add):dest[i]^=value
        assert Z==initial[2*N:]
    assert X==expectedX and Y==expectedY
    # A missing final gather leaves a nonzero arbitrary-input-dependent residue.
    X,Z=initial[:N],initial[2*N:]
    for dest,src,M in ((X,Z,H),(Z,X,V),(X,Z,H)):
        add=linear(M,src)
        for i,value in enumerate(add):dest[i]^=value
    assert Z!=initial[2*N:]
    record.update(centre_commutator_equals_diagonal_mixed_transvection=True,
                  shared_arbitrary_mixed_roles=rm,all_dirty_basis_dimension=2*N+rm,
                  shared_dirty_bank_restored_after_each_data_bank=True,
                  omitted_cleanup_detected=True,
                  scope='Complete scalar identities only; common-frame implementation is not supplied.')
    return record

if __name__=='__main__':
    result=[check(a,b) for a,b in ((5,5),(5,6),(6,6),(6,7),(6,10))]
    Path(__file__).with_name('scalar-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
