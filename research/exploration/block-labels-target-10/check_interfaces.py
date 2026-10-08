"""Exact endpoint-encoding and block-label checks; no new exponent claim."""

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse
import json

if not __debug__:
    raise RuntimeError("This exact checker requires assertions; do not run Python with -O.")


@dataclass(frozen=True)
class QI:
    r: F = F(0)
    i: F = F(0)

    def __add__(self, other):
        if not isinstance(other, QI):
            other = QI(F(other))
        return QI(self.r + other.r, self.i + other.i)

    __radd__ = __add__

    def __neg__(self):
        return QI(-self.r, -self.i)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if not isinstance(other, QI):
            other = QI(F(other))
        return QI(self.r * other.r - self.i * other.i,
                  self.r * other.i + self.i * other.r)

    __rmul__ = __mul__

    def serial(self):
        return [str(self.r), str(self.i)]


ZERO, ONE = QI(), QI(F(1))
A, B = QI(F(1, 2), F(1, 2)), QI(F(1, 2), F(-1, 2))
N = 4  # Two independent address directions; all arithmetic is in Q(i).


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def neg(x):
    return tuple(-a for a in x)


def mul(x, y):
    return tuple(a * b for a, b in zip(x, y))


def phase(x, axis, inverse=False):
    a, b = (B, A) if inverse else (A, B)
    return tuple(a * x[j] + b * x[j ^ (1 << axis)] for j in range(N))


def shift(x, axis):
    return tuple(x[j ^ (1 << axis)] for j in range(N))


def endpoint(z, axis):
    """T(P)(a,b)=(-b,P^-2 a+P^-1 b), P=C on the selected axis."""
    a, b = z
    return neg(b), add(shift(a, axis), phase(b, axis, inverse=True))


def repeat(op, z, count):
    for _ in range(count):
        z = op(z)
    return z


def role_r(z):
    a, b = z
    return neg(b), add(a, b)


def role_s(z, inverse=False):
    a, b = z
    coeff = B if inverse else -B
    return a, add(b, tuple(coeff * x for x in a))


def basis(banks):
    for j in range(banks * N):
        yield tuple(tuple(ONE if j == bank * N + i else ZERO for i in range(N))
                    for bank in range(banks))


def encoded_product(z, w, axis):
    u, v = z
    s, t = w
    us = mul(u, s)
    return neg(us), add(mul(add(v, phase(u, axis, inverse=True)),
                           add(t, phase(s, axis, inverse=True))),
                        phase(us, axis, inverse=True))


def first_difference(a, b):
    for bank, (x, y) in enumerate(zip(a, b)):
        for address, (u, v) in enumerate(zip(x, y)):
            if u != v:
                return {"bank": bank, "address": address,
                        "left": u.serial(), "right": v.serial()}
    return None


def encoding_checks():
    pairs = list(basis(2))
    for axis, z in product(range(2), pairs):
        # T = P^-1 D(P) R D(P)^-1; the outside phase acts on both banks.
        q = (phase(z[0], axis, inverse=True), z[1])
        q = role_r(q)
        q = (phase(q[0], axis), q[1])
        q = tuple(phase(x, axis, inverse=True) for x in q)
        assert endpoint(z, axis) == q
        assert repeat(role_r, z, 3) == tuple(neg(x) for x in z)
        assert repeat(lambda v: endpoint(v, axis), z, 3) == tuple(
            neg(phase(x, axis)) for x in z)
        assert repeat(lambda v: endpoint(v, axis), z, 6) == tuple(
            shift(x, axis) for x in z)
        assert repeat(lambda v: endpoint(v, axis), z, 12) == z
        # For a SINGLE phase direction, relative encoding is a fixed
        # scalar conjugate of a controlled address translation.
        q = role_s(role_r(z), inverse=True)
        q = (q[0], shift(q[1], axis))
        assert role_s(q) == endpoint(z, axis)
    naive_failure = None
    for axis, z, w in product(range(2), pairs, pairs):
        expected = endpoint(tuple(mul(x, y) for x, y in zip(z, w)), axis)
        ez, ew = endpoint(z, axis), endpoint(w, axis)
        assert encoded_product(ez, ew, axis) == expected
        naive = tuple(mul(x, y) for x, y in zip(ez, ew))
        if naive_failure is None:
            naive_failure = first_difference(naive, expected)
    assert naive_failure is not None

    # The preceding reduction does not extend to two active phase factors.
    two_block_failure = None
    for column, z in enumerate(pairs):
        a, b = z
        x_a = shift(shift(a, 0), 1)
        q_b = phase(phase(b, 0, inverse=True), 1, inverse=True)
        true_value = (neg(b), add(x_a, q_b))
        q = role_s(role_r(z), inverse=True)
        q = (q[0], shift(shift(q[1], 0), 1))
        proposed_value = role_s(q)
        diff = first_difference(proposed_value, true_value)
        if diff is not None:
            two_block_failure = {"input_column": column, **diff}
            break
    assert two_block_failure is not None

    # Parent scalar mixer adds the second role pair into the first.
    def mixer(z):
        a, b, c, d = z
        return add(a, c), add(b, d), c, d

    def enc(z):
        return endpoint(z[:2], 0) + endpoint(z[2:], 1)

    witness = None
    for column, z in enumerate(basis(4)):
        diff = first_difference(enc(mixer(z)), mixer(enc(z)))
        if diff is not None:
            witness = {"input_column": column, **diff}
            break
    assert witness is not None
    # Positive control: uniform encodings do commute with this mixer.
    def uniform_enc(z):
        return endpoint(z[:2], 0) + endpoint(z[2:], 0)
    for z in basis(4):
        assert uniform_enc(mixer(z)) == mixer(uniform_enc(z))
    return {"address_dimension": 2, "endpoint_basis_checks": 16,
            "encoded_product_bilinear_basis_checks": 128,
            "naive_product_negative_control": naive_failure,
            "single_direction_translation_factorization_checks": 16,
            "two_block_translation_factorization_negative_control": two_block_failure,
            "different_direction_parent_mixer_commutator": witness,
            "uniform_encoding_parent_mixer_positive_control": True}


def clique_checks():
    result = []
    for h in range(4, 65):
        labels = []
        for start in range(0, h - 3, 4):
            block = sum(1 << i for i in range(start, start + 4))
            labels.extend(block ^ (1 << i) for i in range(start, start + 4))
        assert len(labels) == 4 * (h // 4)
        assert all(x.bit_count() == 3 for x in labels)
        assert all(((x & y).bit_count() % 2) == (i == j)
                   for i, x in enumerate(labels) for j, y in enumerate(labels))
        result.append({"h": h, "orthogonal_clique": len(labels)})
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {"status": "PASS", "new_multiplication_bound_claimed": False,
              "exact_field": "Q(i)", "encoding": encoding_checks(),
              "triple_graph_block_label_cliques": clique_checks()}
    value = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(value)
    print(value, end="")


if __name__ == "__main__":
    main()
