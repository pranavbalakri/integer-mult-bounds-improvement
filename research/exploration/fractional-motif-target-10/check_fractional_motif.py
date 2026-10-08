"""Exact checks for a symmetry-breaking h=8 fit and a scoped block-rank result.

No numerical arithmetic or third-party dependencies. This is a motif check,
not a recursive multiplication certificate. Written with OpenAI Codex assistance.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from math import comb
import json

if not __debug__:
    raise RuntimeError('Assertions must remain enabled')


def binary_rank(rows):
    basis = {}
    for row in rows:
        while row:
            j = row.bit_length()-1
            if j in basis:
                row ^= basis[j]
            else:
                basis[j] = row
                break
    return len(basis)


def dot2(a, b):
    return (a & b).bit_count() % 2


def affine_fit():
    triples = list(combinations(range(8), 3))
    masks = [sum(1 << x for x in t) for t in triples]
    colors = []
    for a, b, c in triples:
        # In F2^3 a triple spans a unique affine plane. Its direction is
        # the kernel of this unique nonzero linear functional.
        choices = [ell for ell in range(1, 8)
                   if dot2(ell, a ^ b) == dot2(ell, a ^ c) == 0]
        assert len(choices) == 1
        colors.append(choices[0])
    assert Counter(colors) == Counter({i: 8 for i in range(1, 8)})
    binary_gram = []
    for i, left in enumerate(masks):
        row = 0
        for j, right in enumerate(masks):
            intersection = (left & right).bit_count()
            rational_entry = int(colors[i] == colors[j])
            binary_entry = intersection % 2
            if i == j:
                assert rational_entry == binary_entry == 1
            else:
                assert rational_entry * binary_entry == 0
                if intersection == 1:
                    assert rational_entry == 0
            row |= binary_entry << j
        binary_gram.append(row)
    assert binary_rank(binary_gram) == 8

    fano = [(0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5),
            (1, 4, 6), (2, 3, 6), (2, 4, 5)]
    assert all(len(set(a) & set(b)) == 1 for a, b in combinations(fano, 2))
    assert len({colors[triples.index(t)] for t in fano}) == 7

    # The displayed support-frame costs. Each color class spans F2^8;
    # consequently every nonzero binary incidence functional meets every
    # color, proving the same costs for every minimal scalar factorization.
    point_support_color_ranks = []
    for point in range(8):
        point_support_color_ranks.append(len({colors[i] for i, t in enumerate(triples)
                                             if point in t}))
    color_support_binary_ranks = [binary_rank(
        [masks[i] for i, c in enumerate(colors) if c == color])
        for color in range(1, 8)]
    assert point_support_color_ranks == [7]*8
    assert color_support_binary_ranks == [8]*7
    for functional in range(1, 1 << 8):
        assert len({colors[i] for i, mask in enumerate(masks)
                    if (functional & mask).bit_count() % 2}) == 7
    for color_support in range(1, 1 << 7):
        assert binary_rank([masks[i] for i, color in enumerate(colors)
                            if color_support & (1 << (color-1))]) == 8
    L = sum(point_support_color_ranks)
    assert L == sum(color_support_binary_ranks) == 56
    v = len(triples)
    return dict(status='PASS: exact rank7 fit, optimal by the Fano clique', h=8,
                indices=v, ordered_pairs_checked=v*v, rational_rank=7, binary_rank=8,
                colors=[dict(triple=list(t), color=c) for t, c in zip(triples, colors)],
                color_class_sizes=dict(Counter(colors)),
                point_support_color_ranks=point_support_color_ranks,
                color_support_binary_ranks=color_support_binary_ranks,
                nonzero_binary_functionals_checked=255,
                nonempty_rational_color_supports_checked=127,
                displayed_centre_cost_per_axis=L,
                data_indices=v*v, two_stage_rank_deficit=v*v-2*v*L,
                scope='Every minimal scalar factorization of these fixed matrices has centre cost at least 56 per axis; no positive two-stage deficit')


def scheme_parameters(h):
    multiplicity = [1, h-1, comb(h, 2)-h, comb(h, 3)-comb(h, 2)]
    disjoint = [(-1)**j * comb(h-3-j, 3-j) for j in range(4)]
    overlap_two = [(3-j)*(h-3-j)-j for j in range(4)]
    det = disjoint[2]*overlap_two[3]-disjoint[3]*overlap_two[2]
    assert det == -2*(h-4)
    alpha = []
    for j in (0, 1):
        beta = Q(disjoint[j]*overlap_two[3]-disjoint[3]*overlap_two[j], det)
        gamma = Q(disjoint[2]*overlap_two[j]-disjoint[j]*overlap_two[2], det)
        alpha.append(1-beta-gamma)
    assert alpha == [Q((h-1)*(h-2)*(9-h), 12), Q((h-2)*(h-3), 4)]
    lower = h-1 if h == 9 else h
    assert multiplicity[2] >= lower and multiplicity[3] >= lower
    assert alpha[1] != 0 and ((alpha[0] == 0) == (h == 9))
    assert sum(multiplicity) == comb(h, 3)
    return dict(h=h, multiplicities=multiplicity,
                disjoint_eigenvalues=disjoint, overlap_two_eigenvalues=overlap_two,
                residual_identity_coefficients=list(map(str, alpha)),
                minimum_invariant_normalized_rank=lower)


def incidence_check(h):
    triples = list(map(frozenset, combinations(range(h), 3)))
    checked = 0
    # Verify the quotient-action identities behind the four eigenspaces by
    # direct counting, independently of the closed eigenvalue formulas.
    for i in range(4):
        for small_tuple in combinations(range(h), i):
            small = frozenset(small_tuple)
            for source in triples:
                for wanted in (0, 2):
                    actual = sum(small <= target and len(source & target) == wanted
                                 for target in triples)
                    if wanted == 0:
                        predicted = comb(h-3-i, 3-i)*int(not (source & small))
                    else:
                        common = len(source & small)
                        predicted = ((3-i)*(h-3) if common == i else
                                     (4-i) if common == i-1 else 0)
                    assert actual == predicted
                    checked += 1
    return dict(h=h, incidence_identities_checked=checked)


def verify():
    fit = affine_fit()
    examples = [scheme_parameters(h) for h in range(7, 101)]
    identities = [incidence_check(h) for h in (7, 8, 9, 10)]
    return dict(status='PASS', affine_rank7_fit=fit,
                invariant_block_check=dict(
                    scope='Blocks depend only on intersection size; arbitrary block size and coefficients',
                    exact_statement='Minimum rank/block-size is h, except h=9 where it is 8, for h>=7',
                    parameter_values_checked=len(examples),
                    examples=[x for x in examples if x['h'] in (7, 8, 9, 10, 23)],
                    incidence_checks=identities),
                multiplication_target='Unmet: no new full child histogram or global exponent')


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
