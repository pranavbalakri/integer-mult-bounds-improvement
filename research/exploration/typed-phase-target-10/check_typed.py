"""Exact checks for typed complex phases. No new multiplication bound."""
from dataclasses import dataclass
from fractions import Fraction as F
from math import comb
from pathlib import Path
import argparse
import json

if not __debug__:
    raise RuntimeError("Run without -O: the exact checker uses assertions.")


@dataclass(frozen=True)
class G:
    r: F = F(0)
    i: F = F(0)

    def __add__(self, y):
        y = y if isinstance(y, G) else G(F(y))
        return G(self.r + y.r, self.i + y.i)

    __radd__ = __add__

    def __neg__(self):
        return G(-self.r, -self.i)

    def __sub__(self, y):
        return self + (-y)

    def __mul__(self, y):
        y = y if isinstance(y, G) else G(F(y))
        return G(self.r*y.r-self.i*y.i, self.r*y.i+self.i*y.r)

    __rmul__ = __mul__

    def inverse(self):
        norm = self.r*self.r+self.i*self.i
        assert norm
        return G(self.r/norm, -self.i/norm)

    def __truediv__(self, y):
        y = y if isinstance(y, G) else G(F(y))
        return self * y.inverse()

    def __pow__(self, n):
        if n < 0:
            return self.inverse() ** (-n)
        value = ONE
        for _ in range(n):
            value = value * self
        return value


ZERO, ONE, I = G(), G(F(1)), G(F(0), F(1))
A, B = (ONE+I)/2, (ONE-I)/2


def plus(x, y):
    return tuple(a+b for a, b in zip(x, y))


def neg(x):
    return tuple(-a for a in x)


def phase(x, axes, power=1):
    for axis, sign in axes:
        k = (power*sign) % 4
        shifted = tuple(x[j ^ (1 << axis)] for j in range(len(x)))
        if k == 1:
            x = tuple(A*a+B*b for a, b in zip(x, shifted))
        elif k == 2:
            x = shifted
        elif k == 3:
            x = tuple(B*a+A*b for a, b in zip(x, shifted))
    return x


def endpoint(z, axes, inverse=False):
    a, b = z
    if inverse:
        return plus(phase(a, axes), phase(b, axes, 2)), neg(a)
    return neg(b), plus(phase(a, axes, 2), phase(b, axes, -1))


def role_fourier(z, inverse=False):
    sign = 1 if inverse else -1
    return tuple(tuple(sum((I ** ((sign*k*j) % 4) * z[j][address]
                            for j in range(4)), ZERO)/2
                       for address in range(len(z[0]))) for k in range(4))


def regular(z, axes):
    z = role_fourier(z)
    z = tuple(phase(x, axes, k) for k, x in enumerate(z))
    return role_fourier(z, inverse=True)


def basis(banks, n):
    for j in range(banks*n):
        yield tuple(tuple(ONE if k*n+x == j else ZERO for x in range(n))
                    for k in range(banks))


def operator_checks():
    p0, p1, relative = ((0, 1),), ((1, 1),), ((0, 1), (1, -1))
    for z in basis(2, 4):
        assert endpoint(endpoint(z, p1, inverse=True), p1) == z
        u, v = z
        h_u = plus(phase(phase(u, p1), p0, 2), neg(phase(u, p0, -1)))
        expected = (u, plus(h_u, phase(phase(v, p1, 2), p0, 2)))
        assert endpoint(endpoint(z, p1, inverse=True), p0) == expected
    for z in basis(4, 4):
        assert role_fourier(role_fourier(z), inverse=True) == z
        assert regular(regular(z, ((1, -1),)), p0) == regular(z, relative)
        assert regular(regular(z, p0), ((0, -1),)) == z

    # F T(P) = D(P) [(F P^-1) R] D(P)^-1, including a batched P.
    full = ((0, 1), (1, 1))
    for selected in (p0, full):
        for z in basis(2, 4):
            actual = tuple(phase(x, full) for x in endpoint(z, selected))
            a, b = phase(z[0], selected, -1), z[1]
            q = (neg(b), plus(a, b))  # R
            q = tuple(phase(phase(x, selected, -1), full) for x in q)
            expected = (phase(q[0], selected), q[1])
            assert actual == expected
    # The regular representation includes controlled increment mod four
    # (use the inverse representation for the positive-shift convention).
    def controlled_increment(bits):
        z, r0, r1 = bits & 1, (bits >> 1) & 1, (bits >> 2) & 1
        return z | ((r0 ^ z) << 1) | ((r1 ^ (z & r0)) << 2)
    values = [controlled_increment(x) for x in range(8)]
    assert len(set(values)) == 8
    affine_failure = values[0] ^ values[1] ^ values[2] ^ values[3]
    assert affine_failure == 4
    return {"relative_two_bank_basis_checks": 8,
            "four_channel_regular_closure_basis_checks": 16,
            "free_chart_similarity_basis_checks": 16,
            "controlled_increment_truth_table": values,
            "non_affine_second_difference": affine_failure}


def mod5(x):
    """Q(i) -> F5, with i mapped to 2; denominators here are powers of 2."""
    def rational(v):
        return v.numerator * pow(v.denominator, -1, 5) % 5
    return (rational(x.r)+2*rational(x.i)) % 5


def rank_mod(rows):
    pivots = {}
    for original in rows:
        row = [x % 5 for x in original]
        for j in range(len(row)):
            value = row[j]
            if not value:
                continue
            if j in pivots:
                row = [(a-value*b) % 5 for a, b in zip(row, pivots[j])]
            else:
                inv = pow(value, -1, 5)
                pivots[j] = [a*inv % 5 for a in row]
                break
    return len(pivots)


def phase_kernel(f, x, y, power):
    k = power % 4
    if k == 0:
        return ONE if x == y else ZERO
    if k == 2:
        return ONE if x ^ y == (1 << f)-1 else ZERO
    w = (x ^ y).bit_count()
    a, b = (A, B) if k == 1 else (B, A)
    return a ** (f-w) * b ** w


def channel_rank_checks():
    receipts = []
    for f in range(1, 6):
        n = 1 << f
        signs = [(-1) ** x.bit_count() for x in range(n)]
        c, d = (-I) ** f, A ** (-f)
        # Flattened normalized kernel of F T(P), indexed by bank and address.
        normalized = [[ZERO for _ in range(2*n)] for _ in range(2*n)]
        for x in range(n):
            for y in range(n):
                base = phase_kernel(f, x, y, 1)
                assert base != ZERO
                normalized[x][n+y] = -ONE
                normalized[n+x][y] = c*signs[x]*signs[y]
                normalized[n+x][n+y] = d if x == y else ZERO
                assert normalized[n+x][y] == phase_kernel(f, x, y, 3)/base
                assert normalized[n+x][n+y] == phase_kernel(f, x, y, 0)/base
                # An explicit n-channel factorization A(x) B(y).
                left = ((-d.inverse(),)*n,
                        tuple(ONE if k == x else ZERO for k in range(n)))
                right = (tuple(c*signs[y]*signs[k] for k in range(n)),
                         tuple(d if k == y else ZERO for k in range(n)))
                for a in range(2):
                    for b in range(2):
                        assert sum((s*t for s, t in zip(left[a], right[b])), ZERO) == normalized[a*n+x][b*n+y]
        rank = rank_mod([[mod5(x) for x in row] for row in normalized])
        assert rank == n
        mode_ranks = []
        for mode in range(4):
            rows = [[mod5(phase_kernel(f, x, y, mode+1)/phase_kernel(f, x, y, 1))
                     for y in range(n)] for x in range(n)]
            mode_ranks.append(rank_mod(rows))
        assert mode_ranks == [1, n, 1, n]
        assert sum(mode_ranks) == 2*n+2
        receipts.append({"active_directions": f,
                         "minimum_common_F_channels_for_two_bank_endpoint": n,
                         "literal_four_channel_lift_common_F_kernel_rank": 2*n+2,
                         "mode_kernel_ranks": mode_ranks})
    assert receipts[2]["minimum_common_F_channels_for_two_bank_endpoint"] == 8
    return receipts


def budget_checks():
    h = 18
    v, m = comb(h, 3), h*h
    n = v*v
    centre = 2*v*(h*h-h+1)
    old = n-centre
    uncorrected = 2*n-centre
    # Costs are normalized by f, the active directions per endpoint.
    assert old == 164832 and uncorrected == 830688
    assert uncorrected-2*n == -centre
    return {"h": h, "m": m, "data_pairs_N": n, "copied_centre_rank": centre,
            "corrected_full_phase_rank_deficit": old,
            "uncorrected_if_incorrectly_counted_as_full_phase": uncorrected,
            "proper_free_chart_full_phase_children_deficit": -centre,
            "literal_relative_T_rank_cost_per_unequal_mixer": 2,
            "literal_regular_lift_rank_cost_per_unequal_mixer": 4,
            "two_bank_mixer_count_needed_to_keep_positive_rank_deficit_strictly_below": uncorrected//2,
            "qualification": "Conditional rank budgets only; no closed global encoded word is claimed."}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {"status": "PASS", "new_multiplication_bound_claimed": False,
              "operator_identities": operator_checks(),
              "pointwise_channel_lift_ranks": channel_rank_checks(),
              "conditional_rank_budget": budget_checks()}
    value = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(value)
    print(value, end="")


if __name__ == "__main__":
    main()
