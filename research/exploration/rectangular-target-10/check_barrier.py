"""Exact obstruction for a specified family, not all multiplication algorithms.

The proof and assumptions are in README.md. No floating-point arithmetic enters
the certificate. Run from any directory with Python 3.10 or later.
"""
from fractions import Fraction as Q


def log_lower(n):
    """A rational lower bound from positive atanh terms after scaling by two."""
    shift = n.bit_length() - 1
    reduced = Q(n, 1 << shift)

    def atanh_series(z):
        return 2 * sum((z ** (2 * j + 1) / (2 * j + 1) for j in range(24)), Q())

    lower = shift * atanh_series(Q(1, 3)) + atanh_series((reduced - 1) / (reduced + 1))
    # Rounding down retains the rigorous lower bound and keeps checks small.
    return Q((lower * 10**9).__floor__(), 10**9)


def verify():
    ceiling_denominator = 1300
    logs = {h: log_lower(h) for h in range(9, 145)}
    smallest_slack = None
    checked = 0
    for a in range(9, 145):
        for b in range(a, 145):
            deficit = 1 - Q(6, a - 2) - Q(6, b - 2)
            if deficit <= 0:
                continue
            entropy_lower = 3 * (a * (logs[b] + 1 - Q(1, b))
                                 + b * (logs[a] + 1 - Q(1, a)))
            slack = entropy_lower - ceiling_denominator * deficit
            assert slack > 0, (a, b, slack)
            checked += 1
            if smallest_slack is None or slack < smallest_slack[0]:
                smallest_slack = (slack, a, b)

    # Infinite tail: log(9)>3*(2/3)+2/17>19/9. Thus for b>=145,
    # 3*b*(log(a)+1-1/a)>9*b>=1305>1300*deficit.
    assert 3 * Q(2, 3) + Q(2, 17) > Q(19, 9)
    assert 9 * 145 > ceiling_denominator
    assert Q(1, ceiling_denominator) < Q(1, 1024) < Q(1, 1023)
    print({"status": "PASS", "finite_pairs": checked,
           "minimum_exact_slack": str(smallest_slack[0]),
           "minimum_slack_pair": list(smallest_slack[1:]),
           "infinite_tail": "max(a,b)>=145",
           "necessary_saving_bound": "a_bit < 1/1300"})


if __name__ == "__main__":
    verify()
