"""Exact obstruction to one delayed cleanup using only mixed point totals.

This does not exclude a clean second bank, extra observations, or a different
word. Uses only the standard library. Written with OpenAI Codex assistance.
"""
from itertools import combinations
from math import comb
import json

if not __debug__:
    raise RuntimeError('Assertions must remain enabled')


def rank(rows):
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


def check(a, b):
    assert a >= 5 and b >= 4
    ta = list(combinations(range(a), 3))
    tb = list(combinations(range(b), 3))
    ma = [sum(1 << i for i in t) for t in ta]
    mb = [sum(1 << i for i in t) for t in tb]
    assert rank(ma) == a and rank(mb) == b
    kernel_triples = [(0, 1, 2), (0, 1, 3), (0, 2, 4), (0, 3, 4)]
    kernel_indices = [ta.index(t) for t in kernel_triples]
    total = 0
    for i in kernel_indices:
        total ^= ma[i]
    assert total == 0
    chosen = 0
    cb_column = sum(((mask & mb[chosen]).bit_count() % 2) << j
                    for j, mask in enumerate(mb))
    assert cb_column & (1 << chosen)  # its diagonal entry is one
    va, vb = len(ta), len(tb)
    residue = sum(cb_column << (i*vb) for i in kernel_indices)
    assert residue != 0
    # The point-product totals vanish on this same input r tensor e_chosen.
    mixed = 0
    for index in kernel_indices:
        for left_point in ta[index]:
            for right_point in tb[chosen]:
                mixed ^= 1 << (left_point*b+right_point)
    assert mixed == 0
    return dict(a=a, b=b, data_coordinates=va*vb,
                witness_first_triples=[list(t) for t in kernel_triples],
                witness_second_triple=list(tb[chosen]),
                witness_mixed_totals_zero=True,
                witness_cleanup_residue_weight=residue.bit_count(),
                late_residue_map_rank=va*b,
                mixed_total_map_rank=a*b,
                minimum_additional_linear_observations=b*(va-a),
                scope='Arbitrary second-bank data; correction sees only mixed totals')


if __name__ == '__main__':
    print(json.dumps(dict(status='PASS', cases=[check(a, b) for a, b in
        ((6, 6), (10, 14), (12, 12), (18, 18), (23, 25))]), indent=2))
