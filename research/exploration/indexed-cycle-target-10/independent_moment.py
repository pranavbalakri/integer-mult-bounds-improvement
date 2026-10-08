"""Independent exact moment audit of the credited, pinned PR84 profile.

Original verifier written with OpenAI Codex assistance. Apache-2.0.
The finite numerical inputs are copied/extracted from Chafik Boukhalfa's PR84;
see README.md, inputs/source.json and the retained UPSTREAM-NOTICE.
This file imports no upstream arithmetic or compiler implementation.
"""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from math import comb, factorial
from pathlib import Path
import json

if not __debug__:
    raise RuntimeError('Run with Python assertions enabled; do not use -O')

HERE = Path(__file__).resolve().parent
GRID = 10**70


def check_source(source):
    pinned = json.loads((HERE/'inputs/source.json').read_text())
    for name, digest in pinned['files'].items():
        assert sha256((source/name).read_bytes()).hexdigest() == digest, name
    closure = json.loads((source/'research/indexed-cycle/SOURCE.json').read_text())
    assert len(closure['files']) == 580
    for name, digest in closure['files'].items():
        assert sha256((source/name).read_bytes()).hexdigest() == digest, name
    return pinned['commit']


def profile():
    m, n = 23*25, comb(23, 3)*comb(25, 3)
    rows = Counter({1:19*n, 21:2*n, 17:2*n, 481:2*n})
    width, loss = 2*n, 0
    pinned = json.loads((HERE/'inputs/source.json').read_text())
    for h in (23, 25):
        name = f'indexed-cycle-profiles-{h}.json'
        raw = (HERE/'inputs'/name).read_bytes()
        assert sha256(raw).hexdigest() == pinned['files']['certificates/'+name]
        p = json.loads(raw)
        repeat, bank = n//p['v'], n//p['v']*p['R']
        width += bank
        loss += repeat*p['loss']
        assert p['h'] == h and p['v'] == comb(h, 3) and p['loss'] == h*(h-1)
        assert sum(i*c for i, c in enumerate(p['blocks'])) == h*p['R']+p['loss']
        for i, count in enumerate(p['blocks']):
            if i and count: rows[i] += repeat*count
        rows[h] += bank
        rows[m-2*h] += bank
        rows[1] += 2*n
        rows[h-2] += 2*n
    mass = sum(i*c for i, c in rows.items())
    assert (width, mass, m*width-mass) == (132466108, 76166165200, 1846900)
    assert mass == m*width-n+loss and max(rows) == 529
    pub = json.loads((HERE/'inputs/published-bit-profile.json').read_text())
    assert {str(k):v for k, v in sorted(rows.items())} == pub['child_multiplicities']
    assert (pub['W'], pub['m'], pub['total_rank']) == (width, m, mass)
    return dict(m=m, W=width, mass=mass, deficit=m*width-mass, rows=rows, published=pub)


def down(x): return Q(x.numerator*GRID//x.denominator, GRID)
def up(x): return -down(-x)


def atanh_bounds(u):
    lo, term = Q(0), u
    for j in range(80):
        lo += term/Q(2*j+1)
        term *= u*u
    return 2*lo, 2*(lo+term/(Q(161)*(1-u*u)))


LN2 = atanh_bounds(Q(1, 3))


def log_bounds(x):
    assert x >= 1
    power = 0
    while x >= 2:
        x /= 2
        power += 1
    lo, hi = atanh_bounds((x-1)/(x+1))
    return down(power*LN2[0]+lo), up(power*LN2[1]+hi)


def exp_bounds(x):
    assert 0 <= x < 1
    lo = sum(x**j/factorial(j) for j in range(13))
    hi = lo + x**13/factorial(13)/(1-x/14)
    return lo, hi


def moment(p, saving):
    lo, hi = Q(0), Q(0)
    for rank, count in p['rows'].items():
        loglo, loghi = log_bounds(Q(p['m'], rank))
        weight = Q(rank*count, p['m']*p['W'])
        lo += weight*exp_bounds(saving*loglo)[0]
        hi += weight*exp_bounds(saving*loghi)[1]
    return down(lo), up(hi)


def audit(source=None):
    p = profile()
    if source is not None:
        check_source(source)
        pub = json.loads((source/'certificates/indexed-cycle-kappa.json').read_text())
        assert pub['bit']['child_multiplicities'] == p['published']['child_multiplicities']
    accepted = Q(p['published']['bit_saving'])
    assert accepted == Q(408390793141, 7812500000000000)
    lo, hi = moment(p, accepted)
    next_lo, next_hi = moment(p, accepted+Q(1, 10**18))
    assert hi < 1 < next_lo
    backoff = Q(522, 10**7)
    back_lo, back_hi = moment(p, backoff)
    assert back_hi < 1
    return dict(source_commit='88ca39571907343a49e97f328971ec7bcd26fbfd',
                arithmetic_only=True, m=p['m'], W=p['W'], total_rank=p['mass'],
                deficit=p['deficit'], child_ranks=len(p['rows']),
                published_saving=str(accepted), published_accept_gap=str(1-hi),
                next_grid_rejection_gap=str(next_lo-1),
                backed_off_saving=str(backoff), backed_off_moment_gap=str(1-back_hi),
                published_gap_float=float(1-hi), backed_off_gap_float=float(1-back_hi))


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, help='Optional pinned PR84 checkout for hash checks')
    args = parser.parse_args()
    print(json.dumps(audit(args.source), indent=2))
