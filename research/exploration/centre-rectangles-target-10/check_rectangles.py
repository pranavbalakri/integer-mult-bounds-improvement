"""Exact local centre-rectangle identities and binary frame checks.

These checks do not construct a complete recursive word or certify a new saving.
No third-party dependencies. Written with OpenAI Codex assistance.
"""
if not __debug__:
    raise RuntimeError("This exact checker requires assertions; do not run Python with -O.")

from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import json
import sys

RESEARCH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RESEARCH / 'independent' / 'complex-twostage'))
from frames import echelon, perp_in, nondeg, ok_res, dot


def vec(t):
    return sum(1 << i for i in t)


def check(h):
    triples = list(combinations(range(h), 3))
    masks = list(map(vec, triples))
    ambient = [1 << i for i in range(h)]
    ones = (1 << h) - 1
    v = len(triples)
    stars = [list(echelon([t for t in masks if t >> i & 1]).values()) for i in range(h)]
    avoids = [list(echelon([t for t in masks if not t >> i & 1]).values()) for i in range(h)]
    assert all(len(U) == h-1 and nondeg(U) and ok_res(U) for U in stars + avoids)
    assert all(ok_res(perp_in(ambient, U)) for U in stars + avoids)

    # A = (r-1)/2, expressed through h point stars with thirds/sixths.
    # Also A = (1/3) sum P_i P_i - (1/6) sum E_i P_i,
    # where E_i is the exclusion indicator.  This latter decomposition
    # exposes a tempting but invalid alternating residual.
    for a in masks:
        for b in masks:
            r = (a & b).bit_count()
            point = sum((Q(1,3) if a >> i & 1 else Q(-1,6)) for i in range(h) if b >> i & 1)
            split = sum(Q(1,3) * bool(a >> i & 1) * bool(b >> i & 1)
                        - Q(1,6) * (not bool(a >> i & 1)) * bool(b >> i & 1)
                        for i in range(h))
            assert point == split == Q(r-1, 2)
            pair = sum(Q(1,3) for ij in combinations(range(h),2)
                       if all(a >> i & 1 and b >> i & 1 for i in ij))
            assert pair == Q(r*(r-1),6)
            assert (a != b or point == pair == 1)
            assert (r != 1 or point == pair == 0)

    valid_cross_reads = 0
    rejected_alternating_reads = 0
    for i in range(h):
        for j in range(h):
            if i != j:
                # E_i -> targets E_j: read frame span(e_j).
                U, read = avoids[i], [1 << j]
                assert all(dot(z,u)==0 for z in read for u in avoids[j])
                residual = perp_in(U,read)
                assert len(residual)==h-2 and ok_res(read) and ok_res(residual)
                valid_cross_reads += 1
                # P_i -> targets P_j: read frame span(ones+e_j).
                U, read = stars[i], [ones ^ (1 << j)]
                assert all(dot(z,u)==0 for z in read for u in stars[j])
                residual = perp_in(U,read)
                assert len(residual)==h-2 and ok_res(read) and ok_res(residual)
                valid_cross_reads += 1
            else:
                # P_i -> targets E_i has cross-Gram rank h-2, but the
                # remaining phase subspace is nonzero alternating.
                residual = perp_in(stars[i],[1 << i])
                assert len(residual)==h-2 and nondeg(residual)
                assert not any(z.bit_count() & 1 for z in residual)
                assert not ok_res(residual)
                rejected_alternating_reads += 1

    # Pair stars: every member has norm 1, distinct members pair to 0;
    # complements and every proper source/side residual are also allowed.
    for pair in combinations(range(h),2):
        U = [vec(pair + (i,)) for i in range(h) if i not in pair]
        assert all(dot(a,b)==(i==j) for i,a in enumerate(U) for j,b in enumerate(U))
        assert ok_res(U) and ok_res(perp_in(ambient,U))
        for t in U:
            residual=[u for u in U if u != t]
            assert all(dot(t,u)==0 for u in residual) and ok_res(residual)

    return dict(h=h, v=v, ordered_scalar_pairs=v*v,
                point_star_centre_rank_cost=h*(h-1),
                inherited_centre_rank_cost=h*h-h+1,
                pair_star_centre_rank_cost=3*v,
                copied_two_stage_deficit_point_stars=v*(v-2*h*(h-1)),
                valid_cross_read_examples=valid_cross_reads,
                alternating_cross_reads_rejected=rejected_alternating_reads,
                scope='Local identities and frame facts only; no complete word or new saving')


def main():
    results=[check(h) for h in (6,8,10,12,14)]
    print(json.dumps(results,indent=2))


# Separate finite certificate for the complete point-star/exclusion dictionary.
def inverse(M):
    n=len(M)
    rows=[list(map(Q,row))+[Q(i==j) for j in range(n)] for i,row in enumerate(M)]
    for j in range(n):
        p=next((i for i in range(j,n) if rows[i][j]),None)
        if p is None: return None
        rows[j],rows[p]=rows[p],rows[j]
        x=rows[j][j]
        rows[j]=[a/x for a in rows[j]]
        for i in range(n):
            if i!=j and rows[i][j]:
                x=rows[i][j]
                rows[i]=[a-x*b for a,b in zip(rows[i],rows[j])]
    return [row[n:] for row in rows]


def dictionary_certificate(h=14):
    assert h==14
    P=[[Q(i==j) for j in range(h)] for i in range(h)]
    E=[[Q(1,3)-Q(i==j) for j in range(h)] for i in range(h)]
    H=[[Q(i==j,2)-Q(1,18) for j in range(h)] for i in range(h)]
    dictionary=P+E
    def allowed(v):
        for u in dictionary:
            j=next(i for i,x in enumerate(u) if x)
            c=v[j]/u[j]
            if c and all(a==c*b for a,b in zip(v,u)): return True
        return False
    def reject(F):
        inv=inverse(F)
        if inv is None: return 'singular'
        G=[[sum(inv[k][i]*H[k][j] for k in range(h)) for j in range(h)] for i in range(h)]
        assert any(not allowed(row) for row in G)
        return 'dual outside dictionary'
    assert all(allowed(v) for v in dictionary)
    # An identity matrix has the all-P factorization: the detector must
    # accept these target rows instead of rejecting every possible input.
    assert all(allowed(row) for row in P)
    one_each=[]
    for k in range(h+1):
        one_each.append(dict(exclusions=k,result=reject([E[i] if i<k else P[i] for i in range(h)])))
    doubled=[]
    # Double index 0, omit index 1. The remaining choices are equivalent
    # under permutations, classified by their number of exclusions.
    for k in range(h-1):
        F=[P[0],E[0]]+[E[i] if i<2+k else P[i] for i in range(2,h)]
        doubled.append(dict(other_exclusions=k,result=reject(F)))
    assert inverse(H) is not None
    # At >=16 terms even the cheapest possible ranks exceed the target.
    assert 16*(h-2)>h*(h-1)
    # At 15 terms with cost <182, at most one term can have cost >=13.
    assert 15*(h-2)+2==h*(h-1)
    # Fourteen cheap terms cover <=14 pair invariants, and one expensive
    # term covers <=h-1; there are C(h,2) nonzero target invariants.
    assert 14+(h-1)<h*(h-1)//2
    return dict(h=h, lower_bound=182, rank=14,
                minimum_term_basis_classes=one_each,
                double_index_basis_classes=doubled,
                fifteen_term_pair_coverage_max=14+h-1,
                required_pair_coverage=h*(h-1)//2,
                identity_factorization_acceptance_control=True,
                scope='P_i/E_i rank-one rectangles only; rational or complex coefficients allowed')


if __name__=='__main__':
    if '--dictionary' in sys.argv:
        print(json.dumps(dictionary_certificate(),indent=2))
    else:
        main()
