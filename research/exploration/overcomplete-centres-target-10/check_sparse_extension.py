"""Exact certificate for one sparse auxiliary direction in the E8 source code.

All source rows may be outside the original code. This certifies the scoped
cost bound proved in README.md, not unrestricted overcomplete optimality.
Standard library only; written with OpenAI Codex assistance.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement, product
from math import gcd, lcm
import json

if not __debug__:
    raise RuntimeError('Assertions must remain enabled')


def roots():
    result = []
    for i, j in combinations(range(8), 2):
        for sign in (-1, 1):
            row = [0]*8
            row[i], row[j] = 2, 2*sign
            result.append(tuple(row))
    for tail in product((-1, 1), repeat=7):
        row = (1,)+tail
        if sum(x < 0 for x in row) % 2 == 0:
            result.append(row)
    return result


def canonical(row):
    row = tuple(row)
    first = next(x for x in row if x)
    return row if first > 0 else tuple(-x for x in row)


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def binary_basis(rows):
    pivots, selected = {}, []
    for original in rows:
        row = original
        while row:
            p = row.bit_length()-1
            if p in pivots:
                row ^= pivots[p]
            else:
                pivots[p] = row
                selected.append(original)
                break
    return selected


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
            inv = pow(row[p], -1, prime)
            basis[p] = [x*inv % prime for x in row]
    return len(basis)


def inverse(matrix):
    n = len(matrix)
    a = [[Q(x) for x in row]+[Q(i == j) for j in range(n)]
         for i, row in enumerate(matrix)]
    for j in range(n):
        p = next(i for i in range(j, n) if a[i][j])
        a[j], a[p] = a[p], a[j]
        scale = a[j][j]
        a[j] = [x/scale for x in a[j]]
        for i in range(n):
            if i != j:
                scale = a[i][j]
                a[i] = [x-scale*y for x, y in zip(a[i], a[j])]
    return [row[n:] for row in a]


def primitive(row):
    denominator = lcm(*(x.denominator for x in row))
    integers = [int(x*denominator) for x in row]
    divisor = gcd(*integers)
    return canonical([x//divisor for x in integers])


def verify():
    R = roots()
    assert len(R) == 120 and len(set(R)) == 120
    gram = [[dot(a, b)//4 for b in R] for a in R]
    assert all(dot(a, b) % 4 == 0 for a in R for b in R)
    assert all(gram[i][i] == 2 for i in range(120))
    assert all(gram[i][j] in (-1, 0, 1) for i in range(120) for j in range(i))
    natural = [sum(((1+x) % 2) << j for j, x in enumerate(row)) for row in gram]
    basis = binary_basis(natural)
    assert len(basis) == 9
    code = []
    for choice in range(512):
        row = 0
        for i, value in enumerate(basis):
            if (choice >> i) & 1:
                row ^= value
        code.append(row)
    assert len(set(code)) == 512
    for row in code[1:]:
        assert rank_mod([r for i, r in enumerate(R) if (row >> i) & 1]) == 8
    weights = Counter(row.bit_count() for row in code)
    assert min(row.bit_count() for row in code[1:]) == 56

    # An explicit simple root basis, in doubled coordinates. Its reflections
    # permute the root lines. Every root has coefficients of one sign in it.
    S = [(1,-1,-1,-1,-1,-1,-1,1), (2,2,0,0,0,0,0,0),
         (-2,2,0,0,0,0,0,0), (0,-2,2,0,0,0,0,0),
         (0,0,-2,2,0,0,0,0), (0,0,0,-2,2,0,0,0),
         (0,0,0,0,-2,2,0,0), (0,0,0,0,0,-2,2,0)]
    assert all(dot(s, s) == 8 for s in S)
    assert all(dot(a, b) <= 0 for a, b in combinations(S, 2))
    root_set = set(R)
    assert all(canonical(s) in root_set for s in S)
    for a in S:
        for r in R:
            value = dot(a, r)
            assert value % 4 == 0
            reflected = tuple(x-(value//4)*y for x, y in zip(r, a))
            assert canonical(reflected) in root_set
    Sinv = inverse(S)
    for r in R:
        coefficients = [sum(Sinv[j][i]*r[j] for j in range(8)) for i in range(8)]
        assert all(x.denominator == 1 for x in coefficients)
        assert all(x >= 0 for x in coefficients) or all(x <= 0 for x in coefficients)

    types = []
    for i in range(8):
        normal = primitive([Sinv[j][i] for j in range(8)])
        inside = [j for j, r in enumerate(R) if dot(normal, r) == 0]
        assert rank_mod([R[j] for j in inside]) == 7
        outside = sum(1 << j for j in range(120) if j not in inside)
        punctured = Counter((word & outside).bit_count() for word in code)
        assert punctured[0] == 1
        minimum = min((word & outside).bit_count() for word in code[1:])
        if len(inside) == 63:
            root_index = next(j for j, r in enumerate(R) if primitive(tuple(map(Q, r))) == normal)
            assert (natural[root_index] & outside) == (1 << root_index)
            assert min((word & outside).bit_count() for word in code
                       if word not in (0, natural[root_index])) == 24
        else:
            assert minimum >= 14
        types.append(dict(deleted_simple_root=i+1, primitive_normal=normal,
                          root_lines_in_hyperplane=len(inside),
                          smallest_nonzero_punctured_weight=minimum,
                          punctured_weight_histogram=dict(sorted(punctured.items()))))

    # Exact second/fourth moments on every E7 root hyperplane. These imply
    # that removing fewer than 27 root lines cannot lower its rational rank.
    moment_entries = 0
    for root_index, root in enumerate(R):
        plane = [r for j, r in enumerate(R) if gram[root_index][j] == 0]
        assert len(plane) == 63
        form = [[8*int(a == b)-root[a]*root[b] for b in range(8)] for a in range(8)]
        for a, b in combinations_with_replacement(range(8), 2):
            assert sum(r[a]*r[b] for r in plane) == 9*form[a][b]
            moment_entries += 1
        for a, b, c, d in combinations_with_replacement(range(8), 4):
            actual = sum(r[a]*r[b]*r[c]*r[d] for r in plane)
            expected = (form[a][b]*form[c][d]+form[a][c]*form[b][d]
                        +form[a][d]*form[b][c])
            assert actual == expected
            moment_entries += 1

    # Sharp overcomplete example: replace one natural row by a singleton and
    # its sum with that row. Ten independent arbitrary rows span the code,
    # and their cost is exactly72, tying the nine-row minimum.
    old_basis = binary_basis([natural[0]]+natural)
    assert len(old_basis) == 9 and old_basis[0] == natural[0]
    dictionary = [1, natural[0]^1]+old_basis[1:]
    assert len(binary_basis(dictionary)) == 10
    costs = [rank_mod([r for i, r in enumerate(R) if (row >> i) & 1]) for row in dictionary]
    assert costs == [1, 7]+[8]*8 and sum(costs) == 72
    return dict(status='PASS', root_lines=120, code_dimension=9,
                all_nonzero_codeword_support_ranks=8,
                code_weight_histogram=dict(sorted(weights.items())),
                simple_reflections_checked=960, simple_basis_root_coefficients_checked=120,
                hyperplane_orbit_representatives=types,
                E7_moment_tensor_entries_checked=moment_entries,
                E7_minimum_rank_dropping_deletion_bound=27,
                sparse_extra_direction_weight_limit=13,
                lower_bound_on_total_source_rank=72,
                sharp_overcomplete_dictionary=dict(rows_hex=[hex(x) for x in dictionary],
                    source_ranks=costs, independent_rows=10, total_cost=72),
                scope='Rows in C+<f> for one f of weight at most13; arbitrary source supports allowed',
                target='Unmet; denser or multiple extra directions and other motifs remain open')


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
