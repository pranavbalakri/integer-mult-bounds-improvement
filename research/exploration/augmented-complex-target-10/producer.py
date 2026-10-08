"""A signed monotone side producer for the augmented triple/difference motif.

Original triple-exclusion DAG: research/independent/complex-twostage/producer.py.
New pieces use incident-edge intervals and deleted-vertex pair sums.
"""
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import importlib.util
import sys

if not __debug__:
    raise RuntimeError("Run without -O: exact assertions are required.")

ROOT = Path(__file__).resolve().parents[2]


def imported(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


old = imported("augmented_old_producer", ROOT / "independent/complex-twostage/producer.py")
fr = imported("augmented_frames", ROOT / "independent/complex-twostage/frames.py")
motif = imported("augmented_motif", ROOT / "exploration/asymmetric-fit-target-10/check_augmented.py")


class Producer:
    def __init__(self, h):
        assert h >= 7 and h % 2
        self.h = h
        self.roots = motif.roots(h)
        self.v = len(self.roots)
        self.triples = list(combinations(range(h), 3))
        self.edges = list(combinations(range(h), 2))
        self.tid = {t: i for i, t in enumerate(self.triples)}
        self.eid = {e: len(self.triples)+i for i, e in enumerate(self.edges)}
        self.pos = [0]+[1 << i for i in range(self.v)]
        self.neg = [0]*(self.v+1)
        self.args = [None]*(self.v+1)
        self.kind = [None]*(self.v+1)
        self.lookup = {}
        self.pieces = []
        self.retained = []
        self.labels = {}
        self.cache = {}
        self.demand = {}
        self.full = tuple(1 << i for i in range(h))
        self.target = [self.canonical(fr.perp_in(self.full, [b])) for _, b, _ in self.roots]

        base = old.NStar3(h+1)
        remap = [0]*len(base.sup)
        for n in range(1, len(base.sup)):
            if base.args[n]:
                a, b = base.args[n]
                remap[n] = self.add(remap[a], remap[b], "triple")
            else:
                t = base.triples[n-1]
                remap[n] = 0 if h in t else self.tid[t]+1
        for target, n, coefficient in base.pieces:
            if h not in target and remap[n]:
                self.pieces.append((self.tid[target], remap[n], Q(coefficient, 2)))
        exclusions = {i: remap[base.results[(i,)]] for i in range(h)}

        # An incident edge at i has the signed value x_(i,j) if i<j,
        # and -x_(j,i) otherwise. Its binary labels are orthonormal.
        intervals = {}
        incident = {}
        for i in range(h):
            points = [j for j in range(h) if j != i]
            leaves = []
            for j in points:
                source = self.eid[tuple(sorted((i,j)))]+1
                leaves.append(source if i < j else self.minus(source))
            table = {}
            for lo in range(len(points)):
                running = 0
                for hi in range(lo+1, len(points)+1):
                    running = self.add(running, leaves[hi-1], ("incident", i))
                    table[lo, hi] = running
            intervals[i] = (points, table)
            incident[i] = table[0, len(points)]
            centre = self.add(exclusions[i], self.minus(incident[i]), ("centre", i))
            self.retained.append((i, centre))

        # Triple targets: retain the old triple part and add signed cuts.
        for target, triple in enumerate(self.triples):
            for i in triple:
                points, table = intervals[i]
                p, q = sorted(points.index(j) for j in triple if j != i)
                for lo, hi in ((0,p), (p+1,q), (q+1,len(points))):
                    if lo < hi:
                        self.pieces.append((target, table[lo,hi], Q(-1,2)))

        # Deleted-vertex pair totals in the fixed-i triple star.
        # prefix[j] contains all pairs below j; suffix[j+1] all pairs above.
        # For the crossing term, row[a][j] contains (a,b) with b>j.
        deleted = {}
        for i in range(h):
            points = [j for j in range(h) if j != i]
            count = len(points)
            edge = {(a,b): self.tid[tuple(sorted((i,points[a],points[b])))]+1
                    for a,b in combinations(range(count),2)}
            kind = ("deleted", i)
            prefix = [0]
            for j in range(count):
                prefix.append(self.add(prefix[-1], self.total([edge[a,j] for a in range(j)], kind), kind))
            suffix = [0]*(count+1)
            for j in range(count-1,-1,-1):
                suffix[j] = self.add(suffix[j+1], self.total([edge[j,b] for b in range(j+1,count)], kind), kind)
            row = {}
            for a in range(count):
                running = 0
                for j in range(count-1,a,-1):
                    row[a,j] = running
                    running = self.add(running, edge[a,j], kind)
            for j, omitted in enumerate(points):
                cross = self.total([row[a,j] for a in range(j)], kind)
                node = self.total([prefix[j],suffix[j+1],cross],kind)
                deleted[i,omitted] = node
                self.demand[node] = self.demand.get(node,0) | (1 << omitted)

        # Difference targets: signed deleted triple stars and edge stars.
        for i,j in self.edges:
            target = self.eid[i,j]
            self.pieces.extend(((target,deleted[i,j],Q(-1,2)),
                                (target,deleted[j,i],Q(1,2))))
            for centre, omitted, coefficient in ((i,j,Q(-1,2)), (j,i,Q(1,2))):
                points, table = intervals[centre]
                k = points.index(omitted)
                parts = [table[lo,hi] for lo,hi in ((0,k),(k+1,len(points))) if lo < hi]
                node = self.total(parts, ("incident",centre))
                self.pieces.append((target,node,coefficient))

        self.prune()
        for n in sorted(self.active, reverse=True):
            if self.kind[n] and self.kind[n][0] == "deleted":
                demand = self.demand.get(n,0)
                assert demand
                for a, _ in self.args[n] or ():
                    if self.kind[a] == self.kind[n]:
                        self.demand[a] = self.demand.get(a,0) | demand
        raw = self.pieces
        self.pieces = []
        for target, n, coefficient in raw:
            self.emit(target,n,coefficient)
        self.prune()
        # Filtering an even-arity triple DAG to odd arity can expose an
        # alternating residual from an odd-sized pair star into its exact
        # coordinate cover. At that consumer, read the star's proper
        # summands separately. Their labels have strictly smaller covers.
        # This preserves the scalar node and its label, but can add uses.
        self.inlined_arguments = 0
        for n in sorted(self.active):
            if not self.args[n]: continue
            pending = list(self.args[n])
            repaired = []
            while pending:
                child, coefficient = pending.pop(0)
                if self.valid_edge(self.label(child),self.label(n)):
                    repaired.append((child,coefficient))
                else:
                    assert self.args[child], (child,n,"unsplittable invalid input edge")
                    self.inlined_arguments += 1
                    pending[:0] = [(a,coefficient*c) for a,c in self.args[child]]
            self.args[n] = tuple(repaired)
        self.prune()

    @staticmethod
    def canonical(vectors):
        return tuple(sorted(fr.echelon(vectors).values()))

    def store(self, pos, neg, args, kind):
        assert not pos & neg
        if not neg and pos.bit_count() == 1:
            return pos.bit_length()
        key = (pos,neg,kind)
        if key not in self.lookup:
            n = len(self.pos)
            self.lookup[key] = n
            self.pos.append(pos); self.neg.append(neg)
            self.args.append(tuple(args)); self.kind.append(kind)
        return self.lookup[key]

    def add(self,a,b,kind):
        if not a: return b
        if not b: return a
        assert not (self.pos[a]|self.neg[a]) & (self.pos[b]|self.neg[b])
        return self.store(self.pos[a]|self.pos[b],self.neg[a]|self.neg[b],((a,1),(b,1)),kind)

    def minus(self,a):
        if not a: return 0
        return self.store(self.neg[a],self.pos[a],((a,-1),),("negative",))

    def total(self, values, kind):
        values = [x for x in values if x]
        if not values: return 0
        if len(values) == 1: return values[0]
        k = len(values)//2
        return self.add(self.total(values[:k],kind),self.total(values[k:],kind),kind)

    def support(self,n):
        value = self.pos[n]|self.neg[n]
        while value:
            low = value & -value
            yield low.bit_length()-1
            value ^= low

    def label(self,n):
        if n in self.labels: return self.labels[n]
        indices = list(self.support(n))
        kind = self.kind[n]
        if len(indices) == 1:
            vectors = [self.roots[indices[0]][1]]
        elif kind == "triple":
            common = set(self.triples[indices[0]])
            covered = set()
            for index in indices:
                common &= set(self.triples[index]); covered.update(self.triples[index])
            if len(common) >= 2 and len(indices) < self.h-2:
                vectors = [self.roots[index][1] for index in indices]
            else:
                vectors = [1 << i for i in covered]
        elif kind[0] == "incident":
            vectors = [self.roots[index][1] for index in indices]
        elif kind[0] == "centre":
            vectors = [1 << j for j in range(self.h) if j != kind[1]]
        elif kind[0] == "negative":
            vectors = self.label(self.args[n][0][0])
        elif kind[0] == "deleted":
            i = kind[1]
            normal = [self.roots[self.eid[tuple(sorted((i,j))) ]][1]
                      for j in range(self.h) if (self.demand[n] >> j) & 1]
            vectors = fr.perp_in(self.full,normal)
        else:
            raise AssertionError(kind)
        value = self.canonical(vectors)
        self.labels[n] = value
        return value

    def valid_edge(self,u,v):
        key = (u,v)
        if key not in self.cache:
            bv = fr.echelon(v)
            self.cache[key] = all(fr.red(bv,x)==0 for x in u) and fr.ok_res(fr.perp_in(v,u))
        return self.cache[key]

    def emit(self,target,n,coefficient):
        if self.valid_edge(self.label(n),self.target[target]):
            self.pieces.append((target,n,coefficient))
            return
        assert self.args[n], (target,n,"invalid source-line incidence")
        for child, scalar in self.args[n]:
            self.emit(target,child,coefficient*scalar)

    def prune(self):
        active = set()
        stack = [n for _,n,_ in self.pieces]+[n for _,n in self.retained]
        while stack:
            n = stack.pop()
            if n in active: continue
            active.add(n)
            stack.extend(a for a,_ in self.args[n] or ())
        self.active = active

    def compile(self):
        users = {n:[] for n in self.active}
        for n in sorted(self.active):
            for position,(child,_) in enumerate(self.args[n] or ()):
                users[child].append(("gate",n,position))
        for i,(_,n,_) in enumerate(self.pieces): users[n].append(("piece",i))
        for i,n in self.retained: users[n].append(("retained",i))
        edge,source,pout,rout,gates = {},{},{},{},[]
        size = 0
        for n in sorted(self.active):
            if self.args[n]:
                inputs = tuple(edge[n,j] for j in range(len(self.args[n])))
                pivot = inputs[0]
            else:
                pivot = size; size += 1
                inputs = (pivot,)
                source[n-1] = pivot
            outputs = (pivot,)+tuple(range(size,size+len(users[n])-1))
            size += len(users[n])-1
            assert set(inputs) & set(outputs) == {pivot}
            assert len(set(inputs)) == len(inputs)
            gates.append((n,inputs,outputs))
            for use,slot in zip(users[n],outputs):
                if use[0] == "gate": edge[use[1],use[2]] = slot
                elif use[0] == "piece": pout[use[1]] = slot
                else: rout[use[1]] = slot
        expected = sum(len(self.args[n])-1 for n in self.active if self.args[n])+len(self.pieces)+len(self.retained)
        assert size == expected and len(source) == self.v
        return {"size":size,"source":source,"pout":pout,"rout":rout,"gates":gates}
