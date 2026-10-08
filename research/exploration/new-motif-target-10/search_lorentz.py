"""Finite Lorentz-root search and independent validation of its saved candidates.

This is a heuristic search, not an exhaustive maximum-clique certificate.
All candidate checks and ranks are exact. Written with OpenAI Codex assistance.
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from random import Random

from check_e8 import binary_basis, rank_mod

if not __debug__:
    raise RuntimeError('Assertions must remain enabled')

HERE = Path(__file__).resolve().parent


def catalog(h, name):
    """Keep the original search's enumeration order for reproducible trials."""
    out = []
    if name == 'positive-level3':
        for i, j in combinations(range(h), 2):
            row = [0]*h
            row[i], row[j] = 1, -1
            out.append((tuple(row), 0))
        for w, k in ((3, 1), (6, 2)):
            for support in combinations(range(h), w):
                row = tuple(int(i in support) for i in range(h))
                out.append((row, k))
        for p in range(h):
            for support in combinations([i for i in range(h) if i != p], 7):
                row = tuple(2 if i == p else int(i in support) for i in range(h))
                out.append((row, 3))
    elif name == 'complete-level4':
        # These are all coefficient patterns through height4 for 3<=h<=12,
        # modulo simultaneous negation. We claim no larger-dimensional coverage.
        assert 3 <= h <= 12
        patterns = {
            0: [(1, -1)],
            1: [(1,)*3],
            2: [(1,)*6],
            3: [(2,)+(1,)*7, (-1,)+(1,)*10],
            4: [(3,)+(1,)*9, (2,)*3+(1,)*6, (2,)*2+(-1,)+(1,)*9],
        }
        for k, group in patterns.items():
            for pattern in group:
                if len(pattern) > h:
                    continue
                coefficients = sorted(set(pattern))
                def visit(j, available, row):
                    if j == len(coefficients):
                        if k == 0 and row.index(1) > row.index(-1):
                            return
                        out.append((tuple(row), k))
                        return
                    c = coefficients[j]
                    for slots in combinations(available, pattern.count(c)):
                        updated = row[:]
                        for a in slots:
                            updated[a] = c
                        visit(j+1, [a for a in available if a not in slots], updated)
                visit(0, list(range(h)), [0]*h)
    else:
        raise ValueError(name)
    assert len(out) == len(set(out))
    assert all(sum(row) == 3*k and sum(x*x for x in row)-k*k == 2
               for row, k in out)
    return out


def inner(left, right):
    x, k = left
    y, ell = right
    return sum(a*b for a, b in zip(x, y))-k*ell


def conflict_graph(roots):
    adjacent = [0]*len(roots)
    for i, x in enumerate(roots):
        for j in range(i):
            t = inner(x, roots[j])
            if t and t % 2 == 0:
                adjacent[i] |= 1 << j
                adjacent[j] |= 1 << i
    return adjacent


def fill(adjacent, rng, initial=0):
    possible = ((1 << len(adjacent))-1) ^ initial
    answer = initial
    current = initial
    while current:
        bit = current & -current
        current -= bit
        possible &= ~adjacent[bit.bit_length()-1]
    while possible:
        current = possible
        choices = []
        best = len(adjacent)+1
        while current:
            bit = current & -current
            current -= bit
            i = bit.bit_length()-1
            degree = (adjacent[i] & possible).bit_count()
            if degree < best:
                best, choices = degree, [i]
            elif degree == best:
                choices.append(i)
        i = rng.choice(choices)
        bit = 1 << i
        answer |= bit
        possible &= ~(adjacent[i] | bit)
    return answer


def validate(h, roots):
    assert len(roots) == len(set(roots))
    assert all(len(row) == h and sum(row) == 3*k and
               sum(x*x for x in row)-k*k == 2 for row, k in roots)
    gram = [[inner(x, y) for y in roots] for x in roots]
    assert all(gram[i][i] == 2 for i in range(len(roots)))
    off = Counter(gram[i][j] for i, j in combinations(range(len(roots)), 2))
    assert all(t == 0 or t % 2 for t in off)
    binary = [sum(((value+1) & 1) << j for j, value in enumerate(row))
              for row in gram]
    assert all((row >> i) & 1 for i, row in enumerate(binary))
    for i, j in combinations(range(len(roots)), 2):
        assert (gram[i][j] == 0) == bool((binary[i] >> j) & 1)
    # G=X(I-J/9)X^T. For h!=9, full rank(X)=h implies rank(G)=h.
    # The nonzero minor modulo101 proves rank_Q(X)>=h.
    coordinate_rank = rank_mod([row for row, k in roots])
    assert h != 9 and coordinate_rank == h
    binary_rank = len(binary_basis(binary)[0])
    canonical = json.dumps(sorted((list(row), k) for row, k in roots),
                           separators=(',', ':')).encode()
    return dict(count=len(roots), rational_rank=h, binary_rank=binary_rank,
                centre_density=f'{len(roots)}/{h*binary_rank}',
                unordered_pairs_checked=len(roots)*(len(roots)-1)//2,
                off_diagonal_histogram={str(k): v for k, v in sorted(off.items())},
                candidate_sha256=sha256(canonical).hexdigest())


def search(h, name, trials, seed):
    roots = catalog(h, name)
    adjacent = conflict_graph(roots)
    rng = Random(seed)
    best = sum(1 << i for i, (row, k) in enumerate(roots) if k in (0, 1))
    for trial in range(trials):
        if trial % 10 == 0:
            initial = 0
        else:
            drop = rng.choice((.03, .08, .15, .3))
            initial = 0
            current = best
            while current:
                bit = current & -current
                current -= bit
                if rng.random() > drop:
                    initial |= bit
            if trial % 3 == 0:
                i = rng.randrange(len(roots))
                initial &= ~adjacent[i]
                initial |= 1 << i
        candidate = fill(adjacent, rng, initial)
        if candidate.bit_count() > best.bit_count():
            best = candidate
    selected = [root for i, root in enumerate(roots) if (best >> i) & 1]
    return dict(dimension=h, catalog=name, catalog_size=len(roots),
                seed=seed, trials=trials, selected_indices=[
                    i for i in range(len(roots)) if (best >> i) & 1],
                validation=validate(h, selected),
                scope='Heuristic lower bound only; no maximum-size or multiplication claim')


def verify_receipts(path, rerun=False):
    saved = json.loads(path.read_text())
    checks = []
    for receipt in saved['runs']:
        h, name = receipt['dimension'], receipt['catalog']
        roots = catalog(h, name)
        assert len(roots) == receipt['catalog_size']
        chosen = [roots[i] for i in receipt['selected_indices']]
        checked = validate(h, chosen)
        assert checked == receipt['validation']
        if rerun:
            fresh = search(h, name, receipt['trials'], receipt['seed'])
            assert fresh == receipt
        checks.append(dict(dimension=h, catalog=name, status='PASS',
                           count=len(chosen), search_repeated=rerun))
    return dict(status='Exact candidate validation PASS', checks=checks,
                scope='Candidate validation does not certify search optimality')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-receipts', action='store_true')
    parser.add_argument('--rerun-search', action='store_true')
    parser.add_argument('--dimension', type=int, default=10)
    parser.add_argument('--catalog', choices=('positive-level3', 'complete-level4'),
                        default='complete-level4')
    parser.add_argument('--trials', type=int, default=1000)
    parser.add_argument('--seed', type=int)
    args = parser.parse_args()
    if args.verify_receipts:
        answer = verify_receipts(HERE/'lorentz-search-receipts.json', args.rerun_search)
    else:
        assert args.trials >= 0
        seed = args.seed if args.seed is not None else (
            1031 if args.catalog == 'positive-level3' else 102938)+args.dimension
        answer = search(args.dimension, args.catalog, args.trials, seed)
    print(json.dumps(answer, indent=2))


if __name__ == '__main__':
    main()
