"""Exact Delsarte LP upper bound for 5-subsets with no intersection three.
Four variables; enumerate every rational vertex, no external solver required.
"""

if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')

from fractions import Fraction as Q
from itertools import combinations
from math import comb
import json

def choose(n,k):
    return comb(n,k) if n>=k>=0 else 0

def eigen(n,k,i,j):
    return sum((-1)**(i-t)*choose(k-t,i-t)*choose(k-j,t)*choose(n-k+t-j,t)
               for t in range(i+1))

def solve(matrix, rhs):
    n=len(rhs); a=[[Q(x) for x in row]+[Q(y)] for row,y in zip(matrix,rhs)]
    for c in range(n):
        p=next((p for p in range(c,n) if a[p][c]),None)
        if p is None:return None
        a[c],a[p]=a[p],a[c];s=a[c][c];a[c]=[x/s for x in a[c]]
        for r in range(n):
            if r!=c and a[r][c]:
                s=a[r][c];a[r]=[x-s*y for x,y in zip(a[r],a[c])]
    return [row[-1] for row in a]

def check_eigenmatrix(n,k=5):
    val=[choose(k,i)*choose(n-k,i) for i in range(k+1)]
    multiplicity=[choose(n,j)-choose(n,j-1) for j in range(k+1)]
    for j in range(k+1):
        for ell in range(k+1):
            inner=sum(Q(eigen(n,k,i,j)*eigen(n,k,i,ell),val[i]) for i in range(k+1))
            assert inner==(Q(choose(n,k),multiplicity[j]) if j==ell else 0)
        for i in range(k+1):
            left=eigen(n,k,1,j)*eigen(n,k,i,j)
            right=(choose(k-i+1,1)*choose(n-k-i+1,1)*eigen(n,k,i-1,j) if i else 0)
            right+=i*(n-2*i)*eigen(n,k,i,j)
            if i<k:right+=(i+1)**2*eigen(n,k,i+1,j)
            assert left==right

def bound(n,k=5):
    check_eigenmatrix(n,k)
    distances=[1,3,4,5]
    val=[choose(k,i)*choose(n-k,i) for i in range(k+1)]
    # inequalities c*x >= rhs; A_0=1, A_2=0.
    cs=[[Q(eigen(n,k,i,j),val[i]) for i in distances] for j in range(1,k+1)]
    bs=[Q(-1)]*k
    for i in range(4):
        cs.append([Q(int(i==j)) for j in range(4)]);bs.append(Q(0))
    best=None;active=None;xx=None
    for inds in combinations(range(len(cs)),4):
        x=solve([cs[i] for i in inds],[bs[i] for i in inds])
        if x is None or any(sum(a*b for a,b in zip(c,x))<b for c,b in zip(cs,bs)):continue
        objective=1+sum(x)
        if best is None or objective>best:best,active,xx=objective,inds,x
    assert best is not None
    # Exact dual y>=0: sum y*c=-1. Then sum x <= -sum y*rhs.
    y=solve([[cs[i][j] for i in active] for j in range(4)],[-1]*4)
    assert y is not None and all(t>=0 for t in y)
    dual=1-sum(yj*bs[i] for yj,i in zip(y,active))
    assert dual==best
    return dict(h=n,k=k,lp_bound=str(best),integer_bound=best.numerator//best.denominator,
                retained_centre_threshold=2*n*(n-1),can_pass_retained_threshold=best.numerator//best.denominator>2*n*(n-1),
                distances=distances,primal=list(map(str,xx)),active=list(active),dual=list(map(str,y)))

if __name__=='__main__':
    rows=[bound(n) for n in range(10,41)]
    for row in rows:print(row['h'],row['lp_bound'],row['integer_bound'],row['retained_centre_threshold'],row['can_pass_retained_threshold'])
    from pathlib import Path
    Path(__file__).with_name('k5-lp.json').write_text(json.dumps(rows,indent=2)+'\n')
