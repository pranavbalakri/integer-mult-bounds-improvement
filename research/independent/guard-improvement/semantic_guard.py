"""Exact guard constants for the inherited two-stage complex phase network.

This certifies finite scalar-prefix envelopes and the numeric guard induction. It
is not an independent proof of the phase network or the multiplication assembly.
Run: python3 independent/guard-improvement/semantic_guard.py [h ...]
"""
from fractions import Fraction as Q
from pathlib import Path
import json
import sys

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / 'complex-twostage'))
from producer import NStar3, compile_roles


def prefix_envelope(c, kk, inverse=False):
    """Every row's absolute coefficient sum, bounded without cancellations.

    Each scratch/data input has row norm 1. The actual values may be arbitrary.
    Mix gates use the exact reversible addition order; the inverse reverses it.
    Signs are irrelevant to an upper envelope. Every scalar intermediate,
    including integer multiplication before halving, is included.
    """
    h = c.h
    val = {('s', i): Q(1) for i in range(kk['size'])}
    val.update({(a, t): Q(1) for a in ('x', 'y') for t in c.triples})
    den = {k: 0 for k in val}
    largest = Q(1)
    largest_den = 0

    def add(dst, src, coeff=Q(1)):
        nonlocal largest, largest_den
        coeff = abs(coeff)
        assert coeff.denominator in (1, 2)
        shift = coeff.denominator.bit_length() - 1
        # Form integer numerator * input first, then divide by 1 or 2.
        # Repeated additions / binary doubling never exceed this numerator norm.
        largest = max(largest, coeff.numerator * val[src])
        largest_den = max(largest_den, den[src] + shift)
        val[dst] += coeff * val[src]
        den[dst] = max(den[dst], den[src] + shift)
        largest = max(largest, val[dst])
        largest_den = max(largest_den, den[dst])

    def mix(reverse):
        for _, ins, outs in (reversed(kk['gates']) if reverse else kk['gates']):
            pivot = ('s', ins[0])
            if not reverse:
                for sl in ins[1:]: add(pivot, ('s', sl))
            for sl in outs[1:]: add(('s', sl), pivot)
            if reverse:
                for sl in ins[1:]: add(pivot, ('s', sl))

    def inject():
        for i, (t, _, coeff) in enumerate(c.pieces):
            add(('y', t), ('s', kk['pout'][i]), Q(coeff, 2))

    def scatter():
        for t in c.triples:
            if h - 1 not in t:
                add(('y', t), ('s', kk['rout'][('*',)]))
                points = t
            else:
                add(('y', t), ('s', kk['rout'][('*',)]), Q(5 - h, 2))
                points = [i for i in range(h - 1) if i not in t]
            for i in points:
                add(('y', t), ('s', kk['rout'][('E', i)]), Q(1, 2))

    def copy():
        for t, sl in kk['src'].items(): add(('s', sl), ('x', t))

    # The logical forward word of hist.forward_ops. The copied second scatter
    # reads the same logical retained value, so it has this same envelope.
    word = [('mix', False), ('scatter', None), ('inject', None), ('mix', True),
            ('copy', None), ('mix', False), ('scatter', None), ('inject', None),
            ('mix', True), ('copy', None)]
    for name, rev in (reversed(word) if inverse else word):
        if name == 'mix': mix(rev ^ inverse)
        elif name == 'scatter': scatter()
        elif name == 'inject': inject()
        else: copy()
    return largest, largest_den


def certify(h=16):
    c = NStar3(h)
    kk = compile_roles(c)
    gf, bf = prefix_envelope(c, kk)
    gi, bi = prefix_envelope(c, kk, inverse=True)
    # Stage-one invocations have disjoint banks, as do stage-two invocations.
    # Uniform incoming envelopes compose multiplicatively. Endpoint correction
    # copies A, applies the inverse rank-one phase child, adds its copy to B,
    # and erases the copy: its only scalar growth is this extra addition.
    G = 2 * gf * gi
    B = bf + bi
    assert G < 2**48 and B <= 8
    # A parent cut has e phase bits + B denominator bits and e/2+log2(G)
    # magnitude bits. We deliberately use the looser common charge 2e+64.
    J = 64
    m = h * h
    rho = Q(h - 1, h)
    C = max(8, (2 * h + (J + h - 1) // h))
    assert Q(C) * (1 - rho) >= 2 + Q(J, m)
    assert C * m >= 8 * m
    # h16: D(e)<=36e per recursive invocation. Input to a disjoint piece
    # has at most d completed C-kernels; preprocessing and tails at most8d,
    # and normalized butterfly conversion at most d+9 extra bits.
    assert C + 10 <= 128
    return dict(h=h, m=m, roles=kk['size'], forward_G=str(gf), inverse_G=str(gi),
                G=str(G), B=B, prefix_charge='2e+64', contraction=str(rho),
                recursive_guard_coefficient=C, layer_guard='128(d+1)', C1=1,
                digit_coefficient=128, precision_coefficient=768, precision_exponent=3,
                guard_numeric_cutoff_b=1, gaussian_numeric_lower_bound='15/8',
                assumptions=['Exact inherited phase-frame network identity',
                             'Valid copied-retained-centre schedule',
                             'Whole-residual children of width w*floor(e/m)',
                             'Bounded remainder axes evaluated individually',
                             'Disjoint top-level active-axis pieces'])


if __name__ == '__main__':
    for h in map(int, sys.argv[1:] or ['16']):
        print(json.dumps(certify(h), indent=2))
