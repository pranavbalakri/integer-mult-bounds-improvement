"""Exact support, common-point, retained-total, rank, and frame-chain checks.

For common point c, the Gram matrix of t_{cij} under G=I-J/9 is
B^T B, where B is the unsigned incidence matrix of the pairs {i,j}.
It is positive semidefinite with kernel exactly ker B. Thus every source span
is nondegenerate, its dimension is rank(B), and nested disjoint-sum supports
produce the inherited orthogonal frame chains. The checks below recompute all
supports independently, check every target complement and retained point total,
and verify all local span ranks directly modulo an odd prime.
"""
import sys
from itertools import combinations
from paired_loo import PairedTreeExcl, histogram
import side_chains


def bit_indices(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length()-1
        mask ^= bit


def rank_mod(vectors, dimension, prime=1000000007):
    basis = {}
    for vector in vectors:
        row = [int(i in vector) for i in range(dimension)]
        for k in range(dimension):
            if not row[k]:
                continue
            if k in basis:
                a = row[k]
                row = [(u-a*v) % prime for u, v in zip(row, basis[k])]
            else:
                a = pow(row[k], prime-2, prime)
                basis[k] = [(a*x) % prime for x in row]
                break
    return len(basis)


def check_leave_one_out():
    cases = 0
    for n in range(11):
        for zero_mask in range(1 << n):
            circuit = PairedTreeExcl.__new__(PairedTreeExcl)
            circuit.support = [0] + [1 << i for i in range(n)]
            circuit.args = [None] * (n+1)
            circuit.lookup = {s: i for i, s in enumerate(circuit.support)}
            values = [0 if zero_mask >> i & 1 else i+1 for i in range(n)]
            total, omitted = circuit.loo(values)
            target = sum(circuit.support[x] for x in values)
            assert circuit.support[total] == target
            assert len(omitted) == n
            for i, node in enumerate(omitted):
                assert circuit.support[node] == target ^ circuit.support[values[i]]
            for node, operands in enumerate(circuit.args):
                if operands:
                    a, b = operands
                    assert not circuit.support[a] & circuit.support[b]
                    assert circuit.support[node] == circuit.support[a] | circuit.support[b]
            cases += 1
    return cases


def check_local_ranks(h):
    local = PairedTreeExcl(h-1)
    count = 0
    for node in sorted(local.active):
        pairs = [local.inputs[i] for i in bit_indices(local.support[node])]
        expected = side_chains.graph_rank(pairs)
        assert rank_mod(pairs, h-1) == expected, (h, node)
        count += 1
    return count


def check(h, builder=histogram):
    certificate, graph, compiled = builder(h)
    triples = graph['trip']
    v = len(triples)
    masks = [sum(1 << p for p in triple) for triple in triples]
    point = [sum(1 << i for i, triple in enumerate(triples) if p in triple) for p in range(h)]
    supports = {}
    cores = {}
    for node in sorted(graph['active']):
        operands = graph['args'][node]
        if operands is None:
            assert 1 <= node <= v
            supports[node] = 1 << (node-1)
            cores[node] = masks[node-1]
        else:
            a, b = operands
            assert a < node and b < node
            assert not supports[a] & supports[b]
            supports[node] = supports[a] | supports[b]
            cores[node] = cores[a] & cores[b]
        assert supports[node] == graph['sup'][node]
        assert cores[node], ('missing common point', node)
        common = (cores[node] & -cores[node]).bit_length()-1
        pairs = [tuple(p for p in triples[i] if p != common) for i in bit_indices(supports[node])]
        assert graph['dn'][node] == side_chains.graph_rank(pairs)
    assert len(graph['outputs']) == 3*v
    for (common, target), node in graph['outputs'].items():
        a, b = (p for p in target if p != common)
        assert supports[node] == point[common] & ~point[a] & ~point[b]
    assert len(graph['retained']) == h
    for common, node in graph['retained'].items():
        assert supports[node] == point[common]
        assert graph['dn'][node] == h-1
    for slot, held in enumerate(compiled['hold']):
        for left, right in zip(held, held[1:]):
            assert supports[left] & ~supports[right] == 0
        if slot in compiled['out']:
            common, target = compiled['out'][slot]
            assert all(len(set(triples[i]) & set(target)) == 1 for i in bit_indices(supports[held[-1]]))
    assert compiled['roles'] == compiled['adds'] + 3*v + h
    assert sum(w*n for w, n in certificate['hist'].items()) == certificate['s']
    assert certificate['W']*certificate['m']-certificate['s'] == v*(v-2*h*(h-1))
    return dict(h=h, roles=compiled['roles'], additions=compiled['adds'], outputs=3*v,
                retained=h, local_ranks=check_local_ranks(h), source_nodes=len(supports),
                slots=len(compiled['hold']))


if __name__ == '__main__':
    print('leave-one-out cases', check_leave_one_out(), flush=True)
    arguments = sys.argv[1:]
    builder = histogram
    if '--reassociate' in arguments:
        from reassociate import histogram as builder
        arguments.remove('--reassociate')
    for h in map(int, arguments or ['23']):
        print(check(h, builder), flush=True)
