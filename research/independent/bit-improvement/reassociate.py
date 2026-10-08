"""Reassociate exact existing disjoint sums to improve side-edge child widths.

For a node of source-span dimension d, replacing its operands of dimensions a,b
changes the complete rank-chain multiset only through
  inner(a), inner(d-a), inner(b), inner(d-b),
provided the active node set is held fixed. These terms include both the changed
input transitions and the changed fan-out slots at the operand nodes.

We select a pair with the same disjoint support union that maximizes the sum of
squared child widths, an integer heuristic for concentrating the child widths.
This heuristic is not a proof of moment optimality: the final exact histogram
and moment certificate establish the actual saving. Unused nodes are pruned.
The emitted DAG is topologically renumbered. No new support or frame is needed.
"""
from collections import Counter
import paired_loo
import rtgm
import cert_rt
from certificate_round5 import inner_children


def improve(graph):
    supports = graph['sup']
    dimensions = graph['dn']
    arguments = list(graph['args'])
    active = graph['active']
    h = graph['h']
    cardinality = [bin(s).count('1') for s in supports]
    by_support = {supports[n]: n for n in active}
    point_stars = [sum(1 << i for i, triple in enumerate(graph['trip']) if c in triple) for c in range(h)]
    buckets = [[] for _ in range(h)]
    common_point = {}
    for node in sorted(active):
        for c, star in enumerate(point_stars):
            if supports[node] & ~star == 0:
                buckets[c].append(node)
                common_point.setdefault(node, c)
    for bucket in buckets:
        bucket.sort(key=lambda n: (-cardinality[n], n))
    square = [sum(w*w for w in inner_children(r, h)) for r in range(h+1)]
    changes = 0
    for node in sorted(active):
        if arguments[node] is None:
            continue
        d = dimensions[node]
        score = lambda pair: sum(square[dimensions[a]] + square[d-dimensions[a]] for a in pair)
        chosen = arguments[node]
        best = score(chosen)
        for a in buckets[common_point[node]]:
            if cardinality[a] >= cardinality[node]:
                continue
            if 2*cardinality[a] < cardinality[node]:
                break
            if supports[a] & ~supports[node]:
                continue
            b = by_support.get(supports[node] ^ supports[a])
            if b is None:
                continue
            candidate = (a, b)
            value = score(candidate)
            if value > best:
                chosen, best = candidate, value
        if chosen != arguments[node]:
            arguments[node] = chosen
            changes += 1
    roots = list(graph['outputs'].values()) + list(graph['retained'].values())
    used = set()
    stack = list(roots)
    while stack:
        node = stack.pop()
        if node in used:
            continue
        used.add(node)
        if arguments[node]:
            stack.extend(arguments[node])
    v = len(graph['trip'])
    sequence = list(range(1, v+1)) + sorted((n for n in used if n > v), key=lambda n: (cardinality[n], n))
    remap = {n: i+1 for i, n in enumerate(sequence)}
    result = dict(graph)
    result['args'] = [None] + [tuple(remap[a] for a in arguments[n]) if arguments[n] else None for n in sequence]
    result['sup'] = [0] + [supports[n] for n in sequence]
    result['dn'] = [0] + [dimensions[n] for n in sequence]
    result['outputs'] = {key: remap[n] for key, n in graph['outputs'].items()}
    result['retained'] = {key: remap[n] for key, n in graph['retained'].items()}
    result['active'] = {remap[n] for n in used}
    result['reassociated_nodes'] = changes
    return result


def build(h):
    return improve(paired_loo.build(h))


def histogram(h):
    graph = build(h)
    compiled = rtgm.compile_(graph)
    ranks = Counter()
    for slot in range(compiled['roles']):
        chain = rtgm.chain_dims(graph, compiled, slot)
        ranks.update(b-a for a, b in zip(chain, chain[1:]))
    certificate = cert_rt.rt_histogram(h, compiled['roles'], ranks)
    return certificate, graph, compiled


if __name__ == '__main__':
    import sys
    import moment
    for h in map(int, sys.argv[1:] or ['23']):
        certificate, graph, compiled = histogram(h)
        saving, root = moment.certify(certificate)
        print(h, compiled['roles'], graph['reassociated_nodes'], saving, root, flush=True)
