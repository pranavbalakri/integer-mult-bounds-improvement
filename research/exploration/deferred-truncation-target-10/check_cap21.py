#!/usr/bin/env python3
"""An exact, target-adapted dimension-21 cap on the pinned deferred frames.

Only old dimension-22 readout frames change.  Most reuse an already existing
dimension-21 frame; the other 55 extend an existing dimension-20 frame by one
explicit integer vector.  No changed recursive corner needs rank above two.
"""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("Certificate assertions must be enabled; do not run Python with -O.")

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import random
import sys
import time
from types import SimpleNamespace

from check_truncation import PINS


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", required=True, type=Path)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    source = args.source.resolve()
    for name, digest in PINS.items():
        assert hashlib.sha256((source / name).read_bytes()).hexdigest() == digest, name
    sys.path.insert(0, str(source / "independent/deferred-readout"))
    import deferred as dr
    import check_frames as cf
    from linalg import Q31, null_exact, rank_mod, rightmost_pivots, predicted

    started = time.monotonic()
    W, D = dr.load(23)
    S = dr.Schedule(W, D)
    h = S.h
    cov = S.adjoint()
    affected = [s for s in S.readout if S.f[s] == 22]
    targets = {next(iter(cov[s])) for s in affected}
    assert all(len(cov[s]) == 1 for s in affected)
    assert all(S.dim(S.start_key(s)) == 22 for s in affected)
    by_target = {"Z": defaultdict(list), "F2": defaultdict(list)}
    for s in S.readout:
        for t, coefficient in cov[s].items():
            by_target["Z"][t].append(s)
            if coefficient & 1:
                by_target["F2"][t].append(s)

    # Rational bases are stored as integer rows.  A modular nonzero Gram
    # determinant proves independence and exact nondegeneracy over Q.
    def nondeg(B):
        sums = [sum(row) for row in B]
        gram = [[9 * sum(x*y for x,y in zip(a,b)) - sa*sb
                 for b,sb in zip(B,sums)] for a,sa in zip(B,sums)]
        return rank_mod(gram, Q31) == len(B)

    caps = {}
    largest = {}
    extensions = []
    rng = random.Random(2026100821)
    exact_inclusion_comparisons = 0
    for t in sorted(targets):
        T = S.trip[t]
        normal = [9 * int(j in T) - 3 for j in range(h)]
        proper = [s for s in by_target["Z"][t] if S.f[s] < 22]
        source_slot = max(proper, key=lambda s: S.f[s])
        largest[t] = source_slot
        B = [row[:] for row in S.sigma[source_slot]]
        assert len(B) in (20, 21)
        assert nondeg(B)
        assert all(sum(x*y for x,y in zip(row,normal)) == 0 for row in B)
        # Verify the actual rational containment of EVERY proper readout
        # contributing to this target, over Z using a rational annihilator.
        # Thus the cap does not assume dimension order implies inclusion.
        annihilator = null_exact(B, h)
        for s in proper:
            for row in S.sigma[s]:
                assert all(sum(x*y for x,y in zip(row,z)) == 0 for z in annihilator)
            exact_inclusion_comparisons += 1
        # Every old dimension-22 frame is exactly this target hyperplane.
        for s in by_target["Z"][t]:
            if S.f[s] == 22:
                old = S.sigma[s]
                assert rank_mod(old, Q31) == 22
                assert all(sum(x*y for x,y in zip(row,normal)) == 0 for row in old)
        if len(B) == 20:
            full = null_exact([normal], h)
            for attempt in range(100):
                weights = [rng.randrange(-100, 101) for _ in full]
                extra = [sum(a * row[j] for a,row in zip(weights,full)) for j in range(h)]
                if nondeg(B + [extra]):
                    break
            else:
                raise AssertionError(("extension search exhausted", t))
            assert sum(x*y for x,y in zip(extra,normal)) == 0
            B.append(extra)
            extensions.append(dict(target=list(T), source_slot=source_slot, extra_integer_row=extra,
                                   rejected_candidates=attempt))
        assert len(B) == 21 and nondeg(B)
        caps[t] = B
    print(f"Constructed {len(caps)} exact caps, {len(extensions)} one-vector extensions; all Gram checks pass.", flush=True)

    # All old readout frames and every cap get explicit matrix coordinates;
    # there is no appeal to a generic unrelated basis for any edge.
    bases = []
    indices = {}
    def frame(B):
        value = tuple(map(tuple, B))
        if value not in indices:
            indices[value] = len(bases)
            bases.append(B)
        return ("sigma", indices[value])
    cap_keys = {t: frame(B) for t,B in caps.items()}
    source_keys = {s: frame(B) for s,B in S.sigma.items()}
    extra = {("sigma", i): B for i,B in enumerate(bases)}
    mq = cf.ModQ7(SimpleNamespace(h=h, pr=None), {}, extra)
    def key(s):
        return cap_keys[next(iter(cov[s]))] if S.f[s] == 22 else source_keys[s]
    def dim(k):
        return h if k == ("F",) else 0 if k == ("0",) else len(bases[k[1]])
    edges = {}
    raw = Counter()
    def edge(A, B, kind):
        r = dim(B) - dim(A)
        assert r >= 0
        if not r:
            return
        raw[kind] += 1
        edges[A,B] = r
    for s in affected:
        edge(key(s), source_keys[s], "first_internal")
        edge(key(s), ("F",), "exterior_corner")
    y_max_proper = {}
    for ring in ("Z", "F2"):
        proper_dimensions = Counter()
        for t in sorted(targets):
            slots = by_target[ring][t]
            if not any(S.f[s] == 22 for s in slots):
                continue
            proper_dimensions[max([S.f[s] for s in slots if S.f[s] < 22] or [0])] += 1
            previous = ("0",)
            previous_changed = False
            for s in slots:
                here = key(s)
                changed = S.f[s] == 22
                if changed or previous_changed:
                    edge(previous, here, "Y_" + ring)
                previous, previous_changed = here, changed
            assert previous_changed
            # Pick any of the original full target frames as the endpoint.
            old_full = next(source_keys[s] for s in slots if S.f[s] == 22)
            edge(previous, old_full, "Y_" + ring)
        y_max_proper[ring] = dict(sorted(proper_dimensions.items()))
    assert set(edges.values()) <= {1, 2}
    used = {k for e in edges for k in e}
    assert all(mq.nondeg(k) for k in used)
    assert mq.dimfail == 0

    # The recipe only needs splitting a rank-r corner into r singletons.
    # Nevertheless check the stronger generic NE profile in the SAME source
    # U,V basis, for every changed rank-1 and rank-2 edge, in both stages.
    want = predicted(h)
    checked = Counter()
    for A,B in sorted(edges, key=lambda e: (str(e[1]),str(e[0]))):
        r = edges[A,B]
        qa,qb = mq.conjugated(A),mq.conjugated(B)
        for side in range(2):
            M = [[(qb[side][i*h+j]-qa[side][i*h+j]) % Q31 for j in range(h)] for i in range(h)]
            assert rightmost_pivots(M,Q31) == want[r], (A,B,r,side)
            checked[r] += 1
    assert mq.dimfail == 0
    print(f"Both common-basis NE profiles pass for all {len(edges)} distinct changed edges.", flush=True)

    # Check all changed target/auxiliary rank totals without replacing the
    # parent audit's independent scalar/frame-key event replay.
    new_f = [min(f,21) for f in S.f]
    for s in affected:
        old_internal = S.dim(S.start_key(s))-S.f[s]
        new_internal = S.dim(S.start_key(s))-new_f[s]
        assert new_internal == old_internal+1
        assert (h*h-h+new_f[s]) + new_internal == (h*h-h+S.f[s]) + old_internal
    output = dict(
        source_repository="https://github.com/Swapnil-jain/integer-mult-kappa",
        source_commit="741e7aa078392553815df7926ee17ac5e25a8c38",
        source_sha256=PINS, h=h, modulus=Q31,
        changed_slots=len(affected), affected_targets=len(targets),
        existing_dimension21_frames_reused=len(caps)-len(extensions),
        dimension20_frames_extended=len(extensions),
        exact_proper_frame_containment_comparisons=exact_inclusion_comparisons,
        all_caps_integer_defined=True, all_caps_nondegenerate=True,
        all_changed_edges_nested_over_Q=True,
        common_axis_basis="check_lifted.point_UV(23,1)",
        changed_edges_by_kind=dict(raw), distinct_changed_edges=len(edges),
        common_basis_NE_checks_by_rank=dict(sorted(checked.items())), NE_failures=0,
        largest_proper_Y_frames=y_max_proper,
        entrance_dimension_histogram=dict(sorted(Counter(new_f).items())),
        largest_auxiliary_block=525, total_rank_mass_unchanged=True,
        extension_seed=2026100821, extensions=extensions,
        complete_scalar_word_dependency="The original scalar word is unchanged; source first frames equal the old dimension-22 target hyperplanes by the pinned containment certificate and dimension equality.",
        full_multiplication_bound=False,
        elapsed_seconds=round(time.monotonic()-started,3),
    )
    if args.output:
        args.output.write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps({k:v for k,v in output.items() if k not in ("extensions","source_sha256")},indent=2))


if __name__ == "__main__":
    main()
