"""Exact checks for the scoped, rejected five-subset motif. Standard library."""
from fractions import Fraction as Q
from itertools import combinations
from math import comb, isqrt
from pathlib import Path
import json

if not __debug__:
    raise RuntimeError('Run without -O: the certificate uses assertions.')


def rank_mod(rows, prime=1000003):
    pivots = {}
    for row in rows:
        row = [int(x) % prime for x in row]
        for j in range(len(row)):
            if not row[j]:
                continue
            if j in pivots:
                c = row[j]
                row = [(a-c*b) % prime for a,b in zip(row,pivots[j])]
            else:
                inv = pow(row[j],-1,prime)
                pivots[j] = [x*inv % prime for x in row]
                break
    return len(pivots)


def coefficients(s):
    polynomial = Q((s-1)*(s-3),8)
    feature = Q(comb(s,2),4)-Q(3*16*s,128)+Q(3*100,800)
    assert polynomial == feature
    return polynomial


def verify():
    for s in range(6):
        q = coefficients(s)
        if s == 5:
            assert q == 1 and s % 2 == 1
        else:
            assert q == 0 or s % 2 == 0
    assert 78**2-4*3*475 == 384
    assert isqrt(384)**2 != 384 and 38 % 3 != 0

    small_rows = []
    for h in (5,6):
        vertices = tuple(map(frozenset,combinations(range(h),5)))
        pairs = tuple(combinations(range(h),2))
        star = [[int(set(p)<=t) for p in pairs] for t in vertices if 0 in t]
        star_rank = rank_mod(star)
        assert star_rank == (1 if h==5 else 5)
        deficit = 1-Q(2*h*star_rank,len(vertices))
        assert deficit == -9
        small_rows.append(dict(h=h,vertices=len(vertices),star_span=star_rank,deficit=str(deficit)))

    rows = []
    for h in (7,8,9,10):
        pairs = tuple(combinations(range(h),2))
        vertices = tuple(map(frozenset,combinations(range(h),5)))
        feature_rows = [[int(set(p)<=t) for p in pairs] for t in vertices]
        star_rows = [row for t,row in zip(vertices,feature_rows) if 0 in t]
        rank = rank_mod(feature_rows)
        star_rank = rank_mod(star_rows)
        assert rank == comb(h,2)
        assert star_rank == comb(h-1,2)

        form = [[Q(int(p==q),4)-Q(3*len(set(p)&set(q)),128)+Q(3,800)
                 for q in pairs] for p in pairs]
        # Three invariant subspaces of the pair association algebra.
        ones = [Q(1)]*len(pairs)
        difference = [Q(int(0 in p)-int(1 in p)) for p in pairs]
        cycle_coeffs = {(0,1):1,(2,3):1,(0,2):-1,(1,3):-1}
        cycle = [Q(cycle_coeffs.get(p,0)) for p in pairs]
        eigenvalues = (Q(3*h*h-78*h+475,1600),Q(38-3*h,128),Q(1,4))
        for vector,eigen in zip((ones,difference,cycle),eigenvalues):
            assert any(vector) and eigen != 0
            assert [sum((x*y for x,y in zip(row,vector)),Q()) for row in form] == [eigen*x for x in vector]
        # All pair evaluations on these small instances, using the exact
        # reduced incidence counts rather than floating-point matrix products.
        for s in vertices:
            for t in vertices:
                assert coefficients(len(s&t)) == Q(comb(len(s&t),2),4)-Q(3*len(s&t),8)+Q(3,8)
        rows.append(dict(h=h,vertices=len(vertices),rational_rank=rank,star_span=star_rank))

    for h in range(7,40):
        deficit = 1-Q(2*h*comb(h-1,2),comb(h,5))
        assert deficit == 1-Q(120,(h-3)*(h-4))
        assert (deficit>0) == (h>=15)
    # log(105)>4: the positive atanh series is a rigorous lower bound.
    u=Q(105-1,105+1)
    log105_lower=2*sum((u**(2*j+1)/Q(2*j+1) for j in range(80)),Q())
    assert log105_lower>4
    assert comb(15,2)==105 and 4*(4*105-2)==1672
    assert Q(1,1672)<Q(1,1023)
    return dict(status='PASS; rejected candidate, no new exponent',small_cases=small_rows,finite_rank_checks=rows,
                universal_symmetric_bound='a_bit < 1/1672',
                assumptions=['full five-subset polynomial fit','unchanged two-stage data endpoints',
                             'binary point-star centres copied to zero','symmetric factor dimensions'],
                target_attained=False)


if __name__=='__main__':
    result=json.dumps(verify(),indent=2)+'\n'
    print(result,end='')
    Path(__file__).with_name('results.json').write_text(result)
