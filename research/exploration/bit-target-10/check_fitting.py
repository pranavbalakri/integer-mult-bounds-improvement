"""Stdlib-only exact validation of an unsuccessful scalar fitting candidate.
This does not implement or certify an address-frame transfer network.
"""
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import json


def rank_q(matrix):
    rows = [[Q(x) for x in row] for row in matrix]
    rank = 0
    for column in range(len(rows[0]) if rows else 0):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][column]
        rows[rank] = [x/scale for x in rows[rank]]
        for i in range(rank+1, len(rows)):
            scale = rows[i][column]
            if scale:
                rows[i] = [x-scale*y for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def rank_f2(rows):
    basis = {}
    for row in rows:
        while row:
            pivot = row.bit_length()-1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return len(basis)


def check_candidate(path):
    data = json.loads(path.read_text())
    h = data['h']
    triples = list(combinations(range(h), 3))
    assert h == 7 and len(triples) == 35 and data['rank'] == 6
    index = {triple: i for i, triple in enumerate(triples)}
    matrix = [0] * len(triples)
    costs = []
    for rectangle in data['rectangles']:
        sources = [tuple(t) for t in rectangle['sources']]
        targets = [tuple(t) for t in rectangle['targets']]
        assert len(sources) == len(set(sources)) and len(targets) == len(set(targets))
        assert all(t in index for t in sources+targets)
        source_mask = sum(1 << index[t] for t in sources)
        for target in targets:
            matrix[index[target]] ^= source_mask
        gram = [[len(set(s) & set(t))-1 for s in sources] for t in targets]
        rank = rank_q(gram)
        assert rank == rectangle['cross_gram_rank']
        costs.append(rank)
    assert len(costs) == 6 and rank_f2(matrix) == 6
    for i, target in enumerate(triples):
        assert matrix[i] >> i & 1
        for j, source in enumerate(triples):
            if i != j and len(set(target) & set(source)) != 1:
                assert not (matrix[i] >> j & 1)
    cost = sum(costs)
    assert cost == 35
    v = len(triples)
    deficit = v*v - 2*v*cost
    assert deficit == -1225
    # Negative control: flipping an off-diagonal forbidden entry must violate fitting.
    i, j = next((i, j) for i in range(v) for j in range(v)
                if i != j and len(set(triples[i]) & set(triples[j])) != 1)
    corrupted = matrix[:]
    corrupted[i] ^= 1 << j
    assert corrupted[i] >> j & 1
    return dict(h=h, dimension=v, scalar_rank=6, cross_gram_ranks=costs,
                total_cross_gram_rank=cost, optimistic_rank_deficit=deficit,
                conclusion='Scalar fitting succeeds; retained two-stage rank budget fails.')


def check_adjacency_shortcut():
    triples = list(combinations(range(6), 3))
    adjacency = [sum(1 << j for j, t in enumerate(triples) if len(set(s)&set(t)) == 1)
                 for s in triples]
    for i, row in enumerate(adjacency):
        squared = 0
        for j in range(len(triples)):
            if row >> j & 1:
                squared ^= adjacency[j]
        assert squared == 1 << i
    # But distinct adjacent labels have different rank-one projectors. A mixes
    # those labels, so it does not commute with their role-dependent frames.
    i, j = next((i, j) for i in range(len(triples)) for j in range(len(triples))
                if adjacency[i] >> j & 1)
    def projector(triple):
        vector = [Q(int(k in triple)) for k in range(6)]
        covector = [(x-Q(1,3))/2 for x in vector]
        return [[x*y for y in covector] for x in vector]
    assert projector(triples[i]) != projector(triples[j])
    return dict(h=6, adjacency_square='identity over F2',
                endpoint_frame_commutation=False,
                conclusion='Inverting scalar adjacency is not a free endpoint-frame correction.')


if __name__ == '__main__':
    directory = Path(__file__).resolve().parent
    print(json.dumps(check_candidate(directory/'fitting-h7-rank6.json'), indent=2))
    print(json.dumps(check_adjacency_shortcut(), indent=2))
