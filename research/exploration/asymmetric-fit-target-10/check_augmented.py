"""Exact augmented-triple motif and centre-frame checks; no circuit claim.

Uses only the Python standard library. Modular minors certify rational
independence; explicit annihilators certify matching upper bounds.
Written with OpenAI Codex assistance.
"""
from fractions import Fraction as Q
from itertools import combinations
from math import comb
import json

if not __debug__:
    raise RuntimeError('Assertions must remain enabled')


def independent_indices(rows, prime=101):
    basis, indices = {}, []
    for index, row in enumerate(rows):
        row = [x % prime for x in row]
        for pivot in sorted(basis):
            scale = row[pivot]
            if scale:
                row = [(a-scale*b) % prime for a, b in zip(row, basis[pivot])]
        pivot = next((i for i, value in enumerate(row) if value), None)
        if pivot is not None:
            inverse = pow(row[pivot], -1, prime)
            basis[pivot] = [(x*inverse) % prime for x in row]
            indices.append(index)
    return indices


def binary_basis(rows):
    basis, selected = {}, []
    for original in rows:
        row = original
        while row:
            pivot = row.bit_length()-1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                selected.append(original)
                break
    return selected


def rational_gram(left, right):
    # x^T (I-J/9)y; root sums are zero or three.
    return sum(a*b for a, b in zip(left, right))-sum(left)*sum(right)//9


def roots(h):
    answer = []
    for points in combinations(range(h), 3):
        x = tuple(int(i in points) for i in range(h))
        answer.append((x, sum(1 << i for i in points), 1))
    full = (1 << h)-1
    for i, j in combinations(range(h), 2):
        x = tuple(int(k == i)-int(k == j) for k in range(h))
        answer.append((x, full ^ (1 << i) ^ (1 << j), 0))
    return answer


def central_identity(h, all_low_sources):
    # Verify the bilinear identity on a coordinate basis, with exact fractions.
    for i in range(h):
        for j in range(h):
            x = [Q(int(k == i)) for k in range(h)]
            y = [Q(int(k == j)) for k in range(h)]
            t, u = sum(x)/3, sum(y)/3
            wanted = (sum(a*b for a, b in zip(x, y))-t*u)/2
            if all_low_sources:
                actual = sum((-x[k]/2+t/(h-3))*(u-y[k]) for k in range(h))
            else:
                actual = sum((x[-1]-x[k])*(u-y[k])/2 for k in range(h-1))
                actual += (t+Q(3-h, 2)*x[-1])*u
            assert actual == wanted


def subset_frames(h):
    coefficient = Q(9*(h-5), 9-h)
    normals = [tuple(1-3*int(i == j) for i in range(h)) for j in range(h)]
    for i, a in enumerate(normals):
        for j, b in enumerate(normals):
            actual = sum(x*y for x, y in zip(a, b))+Q(sum(a)*sum(b), 9-h)
            assert actual == 9*int(i == j)+coefficient
    # Every principal k-by-k block has eigenvalues 9 (k-1 times) and this one.
    eigenvalues = [9+size*coefficient for size in range(1, h+1)]
    assert all(value != 0 for value in eigenvalues)
    return list(map(str, eigenvalues))


def check(h):
    assert h >= 11 and h % 2
    labels = roots(h)
    vectors = [x for x, _, _ in labels]
    masks = [b for _, b, _ in labels]
    v = len(labels)
    assert v == comb(h+1, 3)
    assert len(independent_indices(vectors)) == h
    assert len(binary_basis(masks)) == h
    assert all(mask.bit_count() % 2 for mask in masks)

    # Check every ordered pair, including both fitting-matrix diagonals.
    for i, (left, bm, _) in enumerate(labels):
        for j, (right, bn, _) in enumerate(labels):
            gram = rational_gram(left, right)
            binary = (bm & bn).bit_count() % 2
            assert binary == (1+gram) % 2
            if i == j:
                assert gram == 2 and binary == 1
            else:
                assert gram in (-1, 0, 1)
                assert not gram or not binary

    bit_source_ranks, complex_source_ranks = [], []
    for i in range(h):
        # f_i = binary coordinate i. Its support spans ker(1-3e_i).
        selected = [x for x, b, _ in labels if (b >> i) & 1]
        assert all(sum(x)-3*x[i] == 0 for x in selected)
        basis = [selected[j] for j in independent_indices(selected)]
        assert len(basis) == h-1
        gram = [[rational_gram(a, b) for b in basis] for a in basis]
        assert len(independent_indices(gram)) == h-1
        bit_source_ranks.append(len(basis))

        # z_i = t-x_i. Its support spans the nondegenerate coordinate hyperplane.
        selected = [b for x, b, t in labels if t-x[i]]
        assert all(not ((b >> i) & 1) for b in selected)
        basis = binary_basis(selected)
        assert len(basis) == h-1
        gram = [sum(((a & b).bit_count() % 2) << j for j, b in enumerate(basis))
                for a in basis]
        assert len(binary_basis(gram)) == h-1
        complex_source_ranks.append(len(basis))

    # All binary linear functionals split into h permutation orbits, by weight.
    # A deficient rank has an explicit annihilator, so the modular lower bound
    # proves exact rational rank. Every other case has full ambient rank.
    orbit_ranks = []
    for weight in range(1, h+1):
        functional = (1 << weight)-1
        selected = [x for x, b, _ in labels if (b & functional).bit_count() % 2]
        rank = len(independent_indices(selected))
        assert rank == (h-1 if weight == 1 else h)
        orbit_ranks.append(dict(weight=weight, rational_support_rank=rank,
                                functionals_in_orbit=comb(h, weight)))
    assert sum(row['functionals_in_orbit'] for row in orbit_ranks) == (1 << h)-1

    total_support = [b for _, b, t in labels if t]
    assert len(binary_basis(total_support)) == h
    central_identity(h, False)
    central_identity(h, True)
    all_low_dyadic = ((h-3) & (h-4)) == 0
    selected_complex_cost = h*(h-1)+(not all_low_dyadic)
    bit_cost = h*(h-1)
    return dict(h=h, indices=v, ordered_pairs_checked=v*v,
                rational_rank=h, binary_rank=h, tensor_ambient=h*h,
                binary_source_frame_ranks=bit_source_ranks,
                complex_low_source_frame_ranks=complex_source_ranks,
                all_source_frames_nondegenerate=True,
                subset_normal_exceptional_eigenvalues=subset_frames(h),
                minimal_scalar_centre_cost=bit_cost,
                dyadic_flag_complex_centre_cost=h*(h-1)+1,
                all_low_source_factorization_dyadic=all_low_dyadic,
                selected_dyadic_complex_centre_cost=selected_complex_cost,
                bit_endpoint_deficit=v*v-2*v*bit_cost,
                complex_endpoint_deficit=v*v-2*v*selected_complex_cost,
                binary_functional_orbits=orbit_ranks,
                scope='Finite fitting matrices and source frames only; no side network or child histogram')


def verify():
    return dict(status='PASS', profiles=[check(h) for h in (11, 13, 15, 17, 19)],
                multiplication_target='Unmet: full compiler and recurrence costs have not been certified',
                precision_scope='Central coefficients remain dyadic; no new global precision bound is claimed')


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
