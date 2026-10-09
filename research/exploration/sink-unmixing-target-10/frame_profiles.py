"""Exact rational frame profiles for the sink-unmixing proposal.

Independent implementation of Gram projection and north-east elimination.
The reversed rectangular coordinate permutation is credited to the pinned
PR34/PR84 geometry; no upstream executable code is imported.
"""
if not __debug__:
    raise RuntimeError('Assertions are required; do not run with -O')
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from math import gcd, lcm


def identity(h): return [[Q(i == j) for j in range(h)] for i in range(h)]


def sub(A, B): return [[x-y for x,y in zip(a,b)] for a,b in zip(A,B)]


@lru_cache(None)
def frame(h, core, cover):
    """Gram-orthogonal projector, conjugated by the fixed basis I+J.

    An envelope has coordinates outside cover zero, equal core coordinates,
    and total coordinate sum three times a core coordinate. The triple-line
    case is handled separately. Formula checked against its defining spaces
    by check_frames.py; it is not an abstract rank-only label.
    """
    c=core.bit_count(); out=cover & ~core; n=out.bit_count()
    if not core and not cover: return tuple(tuple(Q(0) for _ in range(h)) for _ in range(h))
    assert 1 <= c <= 3 and core & ~cover == 0
    if c == 3:
        assert core == cover
        return tuple(tuple(Q((3+((core>>i)&1))*(3*(h+1)*((core>>j)&1)-10),6*(h+1))
                           for j in range(h)) for i in range(h))
    s=3-c; d=s*s+(c-1)*n; den=3*(h+1)*d
    A=[]
    for i in range(h):
        row=[]
        for j in range(h):
            oi=(out>>i)&1; oj=(out>>j)&1
            wi=3+((core>>i)&1); zj=3*(h+1)*((core>>j)&1)-10
            num=s*oi*zj+3*(h+1)*s*wi*oj+n*wi*zj-3*(h+1)*(c-1)*oi*oj
            row.append(Q(i==j and oi)+Q(num,den))
        A.append(tuple(row))
    return tuple(A)


def pivots_exact(A):
    """Integer row operations; every zero decision is over Q, not a prime."""
    M=[]
    for row in A:
        den=lcm(*(x.denominator for x in row))
        M.append([int(x*den) for x in row])
    piv=[]
    for i,row in enumerate(M):
        j=next((j for j in range(len(row)-1,-1,-1) if row[j]),None)
        if j is None: continue
        piv.append((i,j)); p=row[j]
        for k in range(i+1,len(M)):
            q=M[k][j]
            if not q: continue
            g=gcd(p,q); a=p//g; b=q//g
            new=[a*x-b*y for x,y in zip(M[k],row)]
            factor=gcd(*new)
            if factor>1: new=[x//factor for x in new]
            M[k]=new
    return piv


def profile_pivots(piv):
    ans=[]
    for k,(i,j) in enumerate(piv):
        if k and (i,j)==(piv[k-1][0]+1,piv[k-1][1]+1): ans[-1]+=1
        else: ans.append(1)
    return Counter(ans)


def profile(A): return profile_pivots(pivots_exact(A))


def reversed_coordinates(a):
    """PR34 reversed a by (a+2) basis, independently assembled constraints."""
    b=a+2; m=a*b; d=2*a+1
    left=list(range(a))+[a-1]+list(range(a))
    right=list(range(a))+[0]+list(range(a))
    partial=[{} for _ in range(b)]
    for labels,offset in ((left,0),(right,m-d)):
        for k,label in enumerate(labels):
            alpha,beta=divmod(offset+k,b)
            assert alpha not in partial[beta] or partial[beta][alpha]==label
            assert label not in partial[beta].values() or partial[beta].get(alpha)==label
            partial[beta][alpha]=label
    permutations=[]
    for p in partial:
        unused=iter(sorted(set(range(a))-set(p.values())))
        permutations.append([p[i] if i in p else next(unused) for i in range(a)])
    assert all(sorted(p)==list(range(a)) for p in permutations)
    coords=[(permutations[k%b][k//b],k%b) for k in range(m)]
    assert [i for i,j in coords[:a]]==list(range(a))
    assert [i for i,j in coords[-a:]]==list(range(a))
    assert [j for i,j in coords[:b]]==list(range(b))
    assert [j for i,j in coords[-b:]]==list(range(b))
    return coords


def complement_profile(a, axis, c):
    """Exact NE profile of I - U_c tensor a dense line projector.

    This uses all corner ranks, with the standard idempotent-complement
    formula. Each non-overlap rank reduces to a submatrix of the local U_c.
    A separate small-dimensional direct-Q replay checks this reduction.
    """
    b=a+2; dims=(a,b); h=dims[axis]; m=a*b
    coords=reversed_coordinates(a); labels=[z[axis] for z in coords]
    X=frame(h,1<<c,(1<<h)-1); Y=sub(identity(h),X)
    assert Y[0][0] and all(Y[i][0] for i in range(h)) and all(Y[0][j] for j in range(h))
    u=[Y[i][0] for i in range(h)]; v=[Y[0][j]/Y[0][0] for j in range(h)]
    assert all(Y[i][j]==u[i]*v[j] for i in range(h) for j in range(h))
    weights=[u[i]*v[i] for i in range(h)]; assert sum(weights)==1
    prefix=[0]; suffix=[0]*(m+1)
    for i in labels: prefix.append(prefix[-1]|(1<<i))
    for k in range(m-1,-1,-1): suffix[k]=suffix[k+1]|(1<<labels[k])
    @lru_cache(None)
    def rank(S,T):
        if not S or not T:return 0
        common=(S&T).bit_count()
        if S&~T and T&~S:return common+1
        if S!=T:return common
        return common-int(sum(weights[i] for i in range(h) if S>>i&1)==1)
    full=(1<<h)-1; r=h-1
    # C[i][j] = rank of top i rows and columns j,...,m-1.
    C=[[0]*(m+1) for _ in range(m+1)]
    for i in range(1,m+1):
        for j in range(m):
            if j>=i: C[i][j]=rank(prefix[i],suffix[j])
            else: C[i][j]=i-j-r+rank(prefix[j],full)+rank(full,suffix[i])
            assert 0<=C[i][j]<=min(i,m-j)
    piv=[]
    for i in range(m):
        for j in range(m):
            q=C[i+1][j]-C[i][j]-C[i+1][j+1]+C[i][j+1]
            assert q in (0,1)
            if q:piv.append((i,j))
    assert len(piv)==m-r
    assert len({i for i,j in piv})==len(piv)==len({j for i,j in piv})
    return profile_pivots(piv),piv


def complement_from_local(a,axis,X):
    """General exact complement from the local projector's NE pivot set.

    The first/last h physical coordinates contain every local label once,
    in order. Thus every prefix/suffix set is an initial/final local interval
    or the full set. No large-matrix elimination or prime sampling is used.
    """
    h=(a,a+2)[axis];m=a*(a+2)
    coords=reversed_coordinates(a);labels=[z[axis] for z in coords]
    assert labels[:h]==list(range(h)) and labels[-h:]==list(range(h))
    local=pivots_exact(X);r=len(local)
    # ranks of top i local rows, rightmost h-j local columns
    small=[[sum(p<i and q>=j for p,q in local) for j in range(h+1)] for i in range(h+1)]
    C=[[0]*(m+1) for _ in range(m+1)]
    for i in range(1,m+1):
        for j in range(m):
            if j>=i:C[i][j]=small[min(i,h)][max(0,j-(m-h))]
            else:C[i][j]=i-j-r+small[min(j,h)][0]+small[h][max(0,i-(m-h))]
    piv=[]
    for i in range(m):
        for j in range(m):
            q=C[i+1][j]-C[i][j]-C[i+1][j+1]+C[i][j+1]
            assert q in (0,1)
            if q:piv.append((i,j))
    assert len(piv)==m-r
    return profile_pivots(piv),piv
