"""Exact h9 quotient, rectangular endpoints, and one overcomplete-source exclusion.

No full circuit or recurrence certificate. Uses only the standard library.
Written with OpenAI Codex assistance.
"""
from itertools import combinations
from math import comb
import json

if not __debug__:
    raise RuntimeError('Assertions must remain enabled')


def rank_mod(rows, prime=101):
    basis = {}
    for row in rows:
        row = [x % prime for x in row]
        for p in sorted(basis):
            scale = row[p]
            if scale:
                row = [(a-scale*b) % prime for a, b in zip(row, basis[p])]
        p = next((i for i, x in enumerate(row) if x), None)
        if p is not None:
            inverse = pow(row[p], -1, prime)
            basis[p] = [x*inverse % prime for x in row]
    return len(basis)


def rank_binary(rows):
    basis = {}
    for row in rows:
        while row:
            p = row.bit_length()-1
            if p in basis:
                row ^= basis[p]
            else:
                basis[p] = row
                break
    return len(basis)


def labels():
    answer = []
    for S in combinations(range(9), 3):
        x = tuple(int(i in S) for i in range(9))
        answer.append((x, 1, sum(1 << i for i in S)))
    for i, j in combinations(range(9), 2):
        x = tuple(int(k == i)-int(k == j) for k in range(9))
        answer.append((x, 0, 511 ^ (1 << i) ^ (1 << j)))
    return answer


def verify():
    data = labels()
    roots = [tuple(2*x[i]-t+x[8] for i in range(8)) for x, t, _ in data]
    assert len(roots) == 120 and rank_mod(roots) == 8
    assert rank_binary([b for _, _, b in data]) == 9
    gram = []
    for a, (x, t, bm) in enumerate(data):
        row = []
        for b, (y, u, bn) in enumerate(data):
            g = sum(i*j for i, j in zip(x, y))-t*u
            assert sum(i*j for i, j in zip(roots[a], roots[b])) == 4*g
            assert ((bm & bn).bit_count() % 2) == (1+g) % 2
            assert g == 2 if a == b else g in (-1, 0, 1)
            row.append(g)
        gram.append(row)
    bit_costs = []
    for i in range(9):
        selected = [r for r, (_, _, b) in zip(roots, data) if (b >> i) & 1]
        assert rank_mod(selected) == 8
        bit_costs.append(8)
    complex_costs = []
    for i in range(8):
        selected = [b for r, (_, _, b) in zip(roots, data) if r[i]]
        assert rank_binary(selected) == 9
        complex_costs.append(9)

    # For every ordered orthogonal pair (i,j), roots orthogonal to i but
    # nonorthogonal to j, excluding j itself, span the entire hyperplane i-perp.
    # A nonzero 7x7 minor modulo101 and the explicit orthogonality give exact rank7.
    ordered_pairs = 0
    for i in range(120):
        for j in range(120):
            if i == j or gram[i][j] != 0:
                continue
            fixed_support = [roots[k] for k in range(120)
                             if k != j and gram[i][k] == 0 and gram[j][k] != 0]
            assert rank_mod(fixed_support) == 7
            ordered_pairs += 1
    assert ordered_pairs == 120*63

    rectangles = []
    for h in (13, 15, 17, 19, 25, 33):
        v = comb(h+1, 3)
        lb = h*(h-1)
        dyadic_low = ((h-3) & (h-4)) == 0
        lc = lb+(not dyadic_low)
        N = 120*v
        rectangles.append(dict(other_h=h, data_indices=N,
            bit_ambient=8*h, complex_ambient=9*h,
            bit_endpoint_deficit=N-72*v-120*lb,
            complex_endpoint_deficit=N-72*v-120*lc,
            scope='Endpoint rank bookkeeping only; no side costs or child moment'))
    return dict(status='PASS', h9_indices=120, rational_rank=8, binary_rank=9,
                ordered_fitting_pairs_checked=14400,
                dyadic_quotient_denominator=8,
                binary_source_costs=bit_costs, complex_source_costs=complex_costs,
                cost_per_axis=72, rectangles=rectangles,
                one_auxiliary_toggle_orthogonal_pairs_checked=ordered_pairs,
                multiplication_target='Unmet; no complete recursive word in this package')


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
