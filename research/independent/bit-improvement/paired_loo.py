"""A paired-tree leave-one-out producer for the retained-total bit network.

Only the addition DAG changes: inherited frame and moment contracts are unchanged.
Every operation is addition of disjoint positive supports. Empty inputs are removed,
nontrivial lists are rotated last-first, adjacent pairs are reduced, the pair totals
are recursively excluded, and the surviving mate is restored in each output.

The original circuit is not modified. build() installs the subclass only while the
upstream builder runs, then restores the original class, even if construction fails.
"""
from pathlib import Path
import sys
from collections import Counter

UPSTREAM = Path(__file__).resolve().parents[1] / 'two-stage-bit'
sys.path.insert(0, str(UPSTREAM))
import sidegen
import rtgm
import cert_rt

_OriginalExcl = sidegen.Excl


class PairedTreeExcl(_OriginalExcl):
    def loo(self, values):
        nonzero = [(i, x) for i, x in enumerate(values) if x]
        total, omitted = self._nonzero_loo([x for _, x in nonzero])
        output = [total] * len(values)
        for (i, _), x in zip(nonzero, omitted):
            output[i] = x
        return total, output

    def _nonzero_loo(self, values):
        n = len(values)
        if n <= 2:
            return self.total(values), [self.total(values[:i] + values[i+1:]) for i in range(n)]
        rotated = values[-1:] + values[:-1]
        pairs = [rotated[i:i+2] for i in range(0, n, 2)]
        total, outside = self.loo([self.total(pair) for pair in pairs])
        omitted = []
        for pair, remainder in zip(pairs, outside):
            omitted.extend(self.add(remainder, self.total(pair[:j] + pair[j+1:])) for j in range(len(pair)))
        return total, omitted[1:] + omitted[:1]


def build(h):
    previous = sidegen.Excl
    try:
        sidegen.Excl = PairedTreeExcl
        return rtgm.build(h, retain='search')
    finally:
        sidegen.Excl = previous


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
    import moment
    for arg in sys.argv[1:] or ['23']:
        h = int(arg)
        certificate, graph, compiled = histogram(h)
        saving, numerical_root = moment.certify(certificate)
        print(h, compiled['roles'], saving, numerical_root, flush=True)
