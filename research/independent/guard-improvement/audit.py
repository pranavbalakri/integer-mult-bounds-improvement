"""Independent exact audits of the complex producer and its scalar prefix bound.

This checks the actual producer at all ordered triple pairs, the complete
small-instance two-stage scalar circuit with arbitrary dirty scratch and copied
retained totals, and a positive-envelope bound on every scalar prefix.  It does
not claim to prove the global multiplication theorem or the inherited layout
interfaces.

Run from any directory: python3 path/to/audit.py
"""
from fractions import Fraction as Q
from pathlib import Path
import json
import random
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "independent" / "complex-twostage"))
from producer import NStar3, compile_roles


def scalar_program(c, kk):
    """Signed elementary additions in the ten groups of a forward invocation."""
    h, triples = c.h, c.triples
    tid = {t: i for i, t in enumerate(triples)}
    groups = []

    def mix(inverse=False):
        ops = []
        gates = reversed(kk["gates"]) if inverse else kk["gates"]
        for _, ins, outs in gates:
            pivot = ("s", ins[0])
            if inverse:
                for q in reversed(outs[1:]):
                    ops.append((("s", q), pivot, Q(-1)))
                for q in reversed(ins[1:]):
                    ops.append((pivot, ("s", q), Q(-1)))
            else:
                for q in ins[1:]:
                    ops.append((pivot, ("s", q), Q(1)))
                for q in outs[1:]:
                    ops.append((("s", q), pivot, Q(1)))
        return ops

    def retained(sign):
        ops = []
        for i, t in enumerate(triples):
            last = h - 1 in t
            coeff = Q(5 - h, 2) if last else Q(1)
            ops.append((("y", i), ("s", kk["rout"][("*",)]), sign * coeff))
            points = [j for j in range(h - 1) if j not in t] if last else t
            for j in points:
                ops.append((("y", i), ("s", kk["rout"][("E", j)]), sign * Q(1 if last else -1, 2)))
        return ops

    def inject(sign):
        return [(("y", tid[c.pieces[i][0]]), ("s", sl), sign * Q(c.pieces[i][2], 2))
                for i, sl in kk["pout"].items()]

    def copy(sign):
        return [(("s", sl), ("x", tid[t]), Q(sign)) for t, sl in kk["src"].items()]

    groups = [mix(), retained(-1), inject(-1), mix(True), copy(1),
              mix(), retained(1), inject(1), mix(True), copy(-1)]
    return groups


def matrix_identity(h):
    """Check 2*(retained scatter + side injection) equals 2I, entry by entry."""
    c = NStar3(h)
    v = len(c.triples)
    retained = dict(c.retained)
    pieces = {t: [] for t in c.triples}
    for t, node, coeff in c.pieces:
        pieces[t].append((node, coeff))

    def add_support(row, node, coeff):
        bits = c.sup[node]
        while bits:
            bit = bits & -bits
            row[bit.bit_length() - 1] += coeff
            bits ^= bit

    negative_control = None
    for i, t in enumerate(c.triples):
        row = [0] * v
        last = h - 1 in t
        add_support(row, retained[("*",)], 5 - h if last else 2)
        points = [j for j in range(h - 1) if j not in t] if last else t
        for j in points:
            add_support(row, retained[("E", j)], 1 if last else -1)
        for node, coeff in pieces[t]:
            add_support(row, node, coeff)
        assert row == [2 * (i == j) for j in range(v)], (h, t)
        if negative_control is None and pieces[t]:
            bad = row[:]
            node, coeff = pieces[t][0]
            add_support(bad, node, -2 * coeff)
            negative_control = bad != row
    assert negative_control
    return {"h": h, "ordered_pairs": v * v, "identity": True,
            "flipped_piece_negative_control": negative_control}


def scalar_envelope(h):
    """Absolute elementary-update envelopes; no numerical approximations."""
    c = NStar3(h)
    kk = compile_roles(c)
    groups = scalar_program(c, kk)
    v = len(c.triples)
    bounds = []
    for inverse in (False, True):
        vals = {"s": [Q(1)] * kk["size"], "x": [Q(1)] * v, "y": [Q(1)] * v}
        peak = Q(1)
        ops = [op for group in groups for op in group]
        if inverse:
            ops = list(reversed(ops))
        for dst, src, coeff in ops:
            vals[dst[0]][dst[1]] += abs(coeff) * vals[src[0]][src[1]]
            peak = max(peak, vals[dst[0]][dst[1]])
        bits = 0
        while peak > 2 ** bits:
            bits += 1
        bounds.append({"inverse": inverse, "peak": str(peak), "ceil_log2": bits})
    return {"h": h, "roles": kk["size"], "elementary_updates": sum(map(len, groups)), "bounds": bounds}


def run_scalar_invocation(groups, retained_slots, x, y, scratch, modulus, inverse=False, early_copy=False):
    """Execute literal copied-total scatter on a prime-field scalar instance."""
    vals = {"s": list(scratch), "x": list(x), "y": list(y)}
    temp = None
    order = reversed(range(len(groups))) if inverse else range(len(groups))
    for gi in order:
        if early_copy and not inverse and gi == 4:
            temp = {sl: vals["s"][sl] for sl in retained_slots}
        if gi == 6 and temp is None:
            temp = {sl: vals["s"][sl] for sl in retained_slots}
        ops = reversed(groups[gi]) if inverse else groups[gi]
        for dst, src, coeff in ops:
            if inverse:
                coeff = -coeff
            value = temp[src[1]] if gi == 6 and src[0] == "s" and src[1] in retained_slots else vals[src[0]][src[1]]
            a = coeff.numerator * pow(coeff.denominator, -1, modulus) % modulus
            vals[dst[0]][dst[1]] = (vals[dst[0]][dst[1]] + a * value) % modulus
        if gi == 6:
            temp = None  # copied stream is erased; all original dirty slots remain.
    return vals["x"], vals["y"], vals["s"]


def dirty_scalar_network(h=8, seed=109):
    """Run both complete scalar stages on all v^2 data pairs and arbitrary scratch."""
    c = NStar3(h)
    kk = compile_roles(c)
    groups = scalar_program(c, kk)
    v, p = len(c.triples), 1000000007
    rng = random.Random(seed)
    x0 = [rng.randrange(p) for _ in range(v * v)]
    y0 = [rng.randrange(p) for _ in range(v * v)]
    x, y = x0[:], y0[:]
    retained_slots = set(kk["rout"].values())
    restored = 0
    for stage in (1, 2):
        for fixed in range(v):
            positions = [i * v + fixed for i in range(v)] if stage == 1 else [fixed * v + i for i in range(v)]
            sx = [x[i] if stage == 1 else y[i] for i in positions]
            sy = [y[i] if stage == 1 else x[i] for i in positions]
            dirty = [rng.randrange(p) for _ in range(kk["size"])]
            xx, yy, out = run_scalar_invocation(groups, retained_slots, sx, sy, dirty, p, inverse=stage == 2)
            assert xx == sx and out == dirty
            assert yy == [(b + (1 if stage == 1 else -1) * a) % p for a, b in zip(sx, sy)]
            restored += len(dirty)
            for pos, a, b in zip(positions, xx, yy):
                if stage == 1:
                    x[pos], y[pos] = a, b
                else:
                    y[pos], x[pos] = a, b
    assert x == [(-z) % p for z in y0]
    assert y == [(a + b) % p for a, b in zip(x0, y0)]
    # Scalar part of the charged phase correction: B += copied A; negate A;
    # bank reassembly sends (B, -A) to the original role order.
    corrected_x = [(b + a) % p for a, b in zip(x, y)]
    corrected_y = [(-a) % p for a in x]
    assert corrected_x == x0 and corrected_y == y0

    sx = [rng.randrange(p) for _ in range(v)]
    sy = [rng.randrange(p) for _ in range(v)]
    dirty = [rng.randrange(p) for _ in range(kk["size"])]
    _, bad, _ = run_scalar_invocation(groups, retained_slots, sx, sy, dirty, p, early_copy=True)
    assert bad != [(a + b) % p for a, b in zip(sx, sy)]
    return {"h": h, "data_pairs": v * v, "invocations": 2 * v,
            "dirty_scratch_values_restored": restored, "two_stage_map": "(-y,x+y)",
            "corrected_scalar_map": "(x,y)", "early_copy_negative_control": True}


def main():
    out = {"h16_all_pair_identity": matrix_identity(16),
           "h16_scalar_prefix_bound": scalar_envelope(16),
           "h8_full_dirty_scalar_network": dirty_scalar_network()}
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
