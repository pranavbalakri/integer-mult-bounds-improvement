"""Sparse references for the sink-unmixing scalar word.

This emits finite cleanup plans, not an exponent certificate.  All arithmetic
in the scalar replay is exact over F_2.  A zero reference plan is supported.
"""
if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')

from collections import Counter
from pathlib import Path
import argparse
import hashlib
import json
import random
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'sink-unmixing-target-10'))
import scalar_word as old


def model(h):
    G = old.rtgm.build(h, retain='search')
    C, mid, sources = old.producer(G)
    triples = G['trip']; v = len(triples); R = C['roles']
    signal = [0] * R
    for i, s in sources.items():
        signal[s] = 1 << i
    for a, b, n in mid:
        signal[a] ^= signal[b]
    candidates = []
    for s in range(R):
        n = C['hold'][s][-1]
        label = G['sup'][n]
        assert signal[s] and signal[s] & ~label == 0
        if s in C['out']:
            continue
        common = set(range(h))
        for i in old.bits(label):
            common.intersection_update(triples[i])
        for c in sorted(common):
            candidates.append((s, c, tuple(old.bits(signal[s]))))
    return dict(h=h, G=G, C=C, mid=mid, sources=sources, triples=triples,
                v=v, R=R, signal=signal, candidates=candidates)


def optimize_primary(M, start=None, weights=None):
    """Deterministic improving single-source changes; no optimality claim."""
    v = M['v']; candidates = M['candidates']
    weights = weights or [1] * M['R']
    primary = list(start or [min(t) for t in M['triples']])
    by_source = [{} for _ in range(v)]
    missing = []
    for j, (s, c, signal) in enumerate(candidates):
        missing.append(sum(primary[i] != c for i in signal))
        for i in signal:
            by_source[i].setdefault(c, []).append(j)
    flips = 0
    while True:
        changed = False
        for i in range(v):
            oldc = primary[i]
            lost = sum(weights[candidates[j][0]]
                       for j in by_source[i].get(oldc, []) if missing[j] == 0)
            best = (0, oldc)
            for c in M['triples'][i]:
                if c == oldc:
                    continue
                gain = sum(weights[candidates[j][0]]
                           for j in by_source[i].get(c, []) if missing[j] == 1)
                if gain - lost > best[0]:
                    best = gain - lost, c
            if best[0] > 0:
                newc = best[1]
                for j in by_source[i].get(oldc, []):
                    missing[j] += 1
                for j in by_source[i].get(newc, []):
                    missing[j] -= 1
                primary[i] = newc; flips += 1; changed = True
        if not changed:
            break
    return primary, flips


def choose(M, primary, refs=()):
    refs = set(refs)
    anchors = {}
    for s, c, signal in M['candidates']:
        if all(primary[i] == c or (c, i) in refs for i in signal):
            anchors.setdefault(s, c)
    # No reference is retained unless the selected cleanup word reads it.
    used = {(c, i) for s, c in anchors.items() for i in old.bits(M['signal'][s])
            if primary[i] != c}
    assert used <= refs
    used_primary = {i for s, c in anchors.items() for i in old.bits(M['signal'][s])
                    if primary[i] == c}
    return dict(primary=primary, anchors=anchors, refs=sorted(used),
                used_primary=sorted(used_primary))


def sparse_greedy(M, primary, budget, weights=None):
    """Add a candidate row's missing-reference bundle by marginal density.

    This optimizes a scalar incidence surrogate, not a certified child moment.
    Full lookup of all newly covered rows makes reference sharing explicit.
    """
    weights = weights or [1] * M['R']
    required = [(s, c, frozenset((c, i) for i in signal if primary[i] != c))
                for s, c, signal in M['candidates']]
    refs = set()
    while len(refs) < budget:
        eligible = {s for s, c, req in required if req <= refs}
        pending = [(s, c, req - refs) for s, c, req in required
                   if s not in eligible and len(req - refs) <= budget - len(refs)]
        if not pending:
            break
        # Each entry is a genuine completion bundle.  Evaluate every distinct
        # singleton; cap larger bundles to the shortest candidates per anchor.
        bundles = set()
        for c in range(M['h']):
            opts = sorted((req for s, anchor, req in pending if anchor == c),
                          key=lambda z: (len(z), tuple(sorted(z))))
            bundles.update(opts[:16])
        best = None
        for bundle in sorted(bundles, key=lambda z: (len(z), tuple(sorted(z)))):
            assert bundle
            new = {s for s, c, req in pending if req <= bundle}
            score = sum(weights[s] for s in new)
            key = (score / len(bundle), score, -len(bundle))
            if best is None or key > best[0]:
                best = (key, bundle)
        refs.update(best[1])
    return choose(M, primary, refs)


def audit_plan(M, P, full=False):
    G = M['G']; C = M['C']; triples = M['triples']
    v = M['v']; R = M['R']; primary = P['primary']; anchors = P['anchors']
    refs = P['refs']; B = len(refs)
    assert len(primary) == v and all(c in triples[i] for i, c in enumerate(primary))
    assert refs == sorted(set(refs))
    assert all(c in triples[i] and c != primary[i] for c, i in refs)
    bref = {key: 2 * v + R + j for j, key in enumerate(refs)}
    size = 2 * v + R + B
    middle = [(2 * v + a, 2 * v + b) for a, b, n in M['mid']]
    inject = [(2 * v + s, i) for i, s in sorted(M['sources'].items())]
    injectB = [(bref[c, i], i) for c, i in refs]
    lookup = {t: i for i, t in enumerate(triples)}
    Jcentre = [(v + i, 2 * v + s) for s, c in sorted(C['ret'].items())
               for i, t in enumerate(triples) if c in t]
    Jside = [(v + lookup[t], 2 * v + s) for s, (_, t) in sorted(C['out'].items())]
    J = Jcentre + Jside
    cleanout = [0] * v
    for a, b in J:
        cleanout[a - v] ^= M['signal'][b - 2 * v]
    assert cleanout == [1 << i for i in range(v)]
    early_clean = []; late_clean = []; KB = []
    local_nodes = P.get('local_nodes', {})
    role_order = sorted(range(R), key=lambda s: (
        G['dn'][local_nodes[s]] if s in local_nodes else M['h']-1, s))
    for s in role_order:
        if s in anchors:
            c = anchors[s]
            assert s not in C['out']
            last = G['sup'][C['hold'][s][-1]]
            assert all(c in triples[i] for i in old.bits(last))
            for i in old.bits(M['signal'][s]):
                control = i if primary[i] == c else bref[c, i]
                early_clean.append((2 * v + s, control))
                if control >= 2 * v + R:
                    KB.append((2 * v + s, control))
        else:
            # At D1 every unchanged X is available.  Late cleanup needs no B.
            late_clean.extend((2 * v + s, i) for i in old.bits(M['signal'][s]))
    used = sorted({b for a, b in early_clean if b >= 2 * v + R})
    assert used == list(range(2 * v + R, size))
    actual_primary = sorted({b for a, b in early_clean if b < v})
    assert actual_primary == P['used_primary']
    summary = dict(h=M['h'], v=v, producer_roles=R, reference_roles=B,
                   common_point_cleanup_roles=len(anchors),
                   full_stage_cleanup_roles=R-len(anchors),
                   primary_X_detours=len(actual_primary),
                   reference_dirt_unmix_cnot_count=len(KB),
                   cleanup_rank_upper_bound=M['h']-1,
                   cleanup_by_last_rank=dict(sorted(Counter(
                       G['dn'][C['hold'][s][-1]] for s in anchors).items())),
                   clean_scatter_identity_JAP_equals_I=True,
                   every_early_cleanup_has_legal_common_point_support=True,
                   every_reference_is_used_at_one_common_point=True,
                   late_cleanups_use_unchanged_X_at_full_local_frame=True,
                   status='Exact scalar/incidence plan; no fixed-basis histogram or exponent certified.')
    canonical = json.dumps(dict(primary=primary, anchors=sorted(anchors.items()), refs=refs),
                           sort_keys=True, separators=(',', ':')).encode()
    summary['plan_sha256'] = hashlib.sha256(canonical).hexdigest()
    if not full:
        summary['full_dirty_basis_replayed'] = False
        return summary
    initial = [1 << i for i in range(size)]
    dirt = initial[:]
    for a, b in middle:
        dirt[a] ^= dirt[b]
    D = [0] * v
    for a, b in J:
        D[a-v] ^= dirt[b]
    rmask = ((1 << R) - 1) << (2 * v)
    assert all(row & ~rmask == 0 for row in D)
    early = [(v+i, s) for i, row in enumerate(D) for s in old.bits(row)]
    prefix = early + injectB + inject + middle + J + early_clean + late_clean + injectB
    sink = KB + list(reversed(middle))
    word = prefix + sink
    for transpose in (False, True):
        state = initial[:]
        for a, b in word:
            if transpose:
                a, b = b, a
            state[a] ^= state[b]
        expected = initial[:]
        for i in range(v):
            if transpose:
                expected[i] ^= initial[v+i]
            else:
                expected[v+i] ^= initial[i]
        assert state == expected
    state = initial[:]
    for a, b in prefix:
        state[a] ^= state[b]
    datamask = (1 << (2 * v)) - 1
    assert not any(row & datamask for row in state[2*v:])
    assert state[2*v+R:] == initial[2*v+R:]
    expectedR = dirt[2*v:2*v+R]
    for a, b in KB:
        expectedR[a-2*v] ^= initial[b]
    assert state[2*v:2*v+R] == expectedR
    # A nonzero required cleanup cannot be dropped even with arbitrary dirt.
    bad = initial[:]
    cleanup_start = len(early + injectB + inject + middle + J)
    for index, (a, b) in enumerate(word):
        if index == cleanup_start:
            continue
        bad[a] ^= bad[b]
    forward_expected = initial[:]
    for i in range(v):
        forward_expected[v+i] ^= initial[i]
    assert bad != forward_expected
    if KB:
        bad = initial[:]
        for a, b in prefix + list(reversed(middle)):
            bad[a] ^= bad[b]
        assert bad[2*v:] != initial[2*v:]
    summary.update(full_dirty_basis_replayed=True,
                   complete_public_basis_dimension=size,
                   forward_all_dirty_identity=True,
                   same_chronology_transposed_identity=True,
                   intermediate_auxiliary_map_independent_of_data=True,
                   negative_missing_cleanup_rejected=True,
                   negative_missing_reference_unmix_rejected=bool(KB),
                   word_cnot_count=len(word))
    return summary


def run(h, full=False, sparse=False):
    M = model(h)
    starts = [[min(t) for t in M['triples']], [max(t) for t in M['triples']]]
    starts += [[random.Random(1000+i).choice(t) for i, t in enumerate(M['triples'])]]
    best = None
    for start in starts:
        primary, flips = optimize_primary(M, start)
        P = choose(M, primary)
        key = (len(P['anchors']), -len(P['used_primary']))
        if best is None or key > best[0]:
            best = key, P, flips
    rows = []
    plans = [('least-point-zero', choose(M, starts[0])), ('optimized-zero', best[1])]
    if sparse:
        for budget in sorted(set([max(1,M['v']//16), max(1,M['v']//4)])):
            plans.append((f'sparse-budget-{budget}', sparse_greedy(M, best[1]['primary'], budget)))
    for name, P in plans:
        receipt = audit_plan(M, P, full)
        receipt['plan_name'] = name
        rows.append(receipt)
    return rows


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('h', type=int, nargs='*', default=[5,6,7,8,23])
    parser.add_argument('--full-max', type=int, default=8)
    parser.add_argument('--sparse', action='store_true')
    parser.add_argument('--output', type=Path, default=HERE/'results.json')
    args = parser.parse_args()
    results = [row for h in args.h for row in run(h, h <= args.full_max, args.sparse)]
    args.output.write_text(json.dumps(results, indent=2)+'\n')
    print(json.dumps(results, indent=2))
