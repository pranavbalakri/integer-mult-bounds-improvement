#!/usr/bin/env python3
"""Realized common-cut frames for pinned Jain round seven, using stdlib only.

This is a certificate for the changed frame geometry and inner profiles, not
a new full multiplication theorem.  The unchanged word and frame inclusions
are the pinned source's separately auditable input.  See README.md for the
exact lifting argument from modular computations to rational subspaces.
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


PINS = {
    "certificates/round7/witness_23.json.gz": "b2486aa2bb222bacea6e52162a3920780e45dd8edc5d7cb03eabea336a1c6588",
    "certificates/round7/deferred_23.json.gz": "8c38e947ff9e021e308000dd82bb5e2194eb265d8b75194e22ce783d3e75317c",
    "independent/deferred-readout/deferred.py": "e92ac1db9fce8608cf5302336f05750fa327eda99f7170ed32e7323220ffc164",
    "independent/deferred-readout/linalg.py": "56cc779911a3c5c901bc5194ccaf97f720157a2c4ba602b07ee8642af2676b8b",
    "independent/deferred-readout/check_lifted.py": "c3e87ad6b31f377df549421f53d22b99131d2228dcd08a1956913211383cbd6c",
    "independent/deferred-readout/check_frames.py": "1cb2b8cde2a8b3f12008da0d335dcdd38005cc692367423d9589aa1f4d6306e8",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", type=Path, required=True,
                    help="pinned swapnil-round7 directory, with certificates and independent/")
    ap.add_argument("--codim", type=int, nargs="+", default=[1, 2, 3, 4], choices=[1, 2, 3, 4, 22])
    ap.add_argument("--seed", type=int, default=20261008)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    assert all(q in (1, 2, 3, 4, 22) for q in args.codim)
    source = args.source.resolve()
    for name, expected in PINS.items():
        assert hashlib.sha256((source / name).read_bytes()).hexdigest() == expected, name
    sys.path.insert(0, str(source / "independent/deferred-readout"))
    import deferred as dr
    import check_lifted as cl
    import check_frames as cf
    from linalg import Q31, Echelon, null_mod, rank_mod, predicted, rightmost_pivots

    start = time.monotonic()
    def log(message):
        print(f"[{time.monotonic() - start:7.1f}s] {message}", flush=True)

    W, D = dr.load(23)
    S = dr.Schedule(W, D)
    h = S.h
    pr = cl.Prog(W)
    bas, dn = cl.spans(pr)
    users = cl.users_of(pr, dn)
    fr = cl.Frames(pr, dn, bas, users)
    assert fr.check_structure() == (0, 0, 0, True)
    fr.cobases()
    dims, _ = fr.dims()
    ld, _ = cl.late_frames(fr, pr, users)
    assert dims == S.node_dims
    assert {k[1:]: v for k, v in ld.items()} == S.late_dims
    log("Reconstructed unchanged node and late-copy frame dimensions exactly.")

    # Four integer rows define nested rational ambient cuts.  Each candidate
    # uses the SAME first q rows for every deferred frame and every target.
    rng = random.Random(args.seed)
    cut = [[rng.randrange(-10**6, 10**6 + 1) for _ in range(h)] for _ in range(4)]
    assert rank_mod(cut, Q31) == 4
    old_frames = []
    frame_index = {}
    slot_index = {}
    for s in S.readout:
        key = tuple(map(tuple, S.sigma[s]))
        if key not in frame_index:
            frame_index[key] = len(old_frames)
            old_frames.append(S.sigma[s])
        slot_index[s] = frame_index[key]
    assert all(rank_mod(B, Q31) == len(B) for B in old_frames)

    cov = S.adjoint()
    by_target = {"Z": defaultdict(list), "F2": defaultdict(list)}
    for s in S.readout:
        for t, c in cov[s].items():
            by_target["Z"][t].append(s)
            if c & 1:
                by_target["F2"][t].append(s)
    original_extra = {("v", s): B for s, B in S.vstart.items()}
    expected = predicted(h)
    # One matrix U and one matrix V, from the same source integer point,
    # are reused for every new frame, edge, and q; no per-edge basis choices.
    U, V = cl.point_UV(h, 1)
    basis_hash = hashlib.sha256(json.dumps([U, V], separators=(",", ":")).encode()).hexdigest()
    results = []
    for q in args.codim:
        # q=22 is a distinct deterministic endpoint: H=span(e_0).
        # Every target normal has e_0 coefficient 6 or -3, so every sigma
        # contained in a target hyperplane intersects this line trivially.
        applied_cut = ([[int(i == j) for j in range(h)] for i in range(1,h)]
                       if q == 22 else cut[:q])
        def intersection(B, C=applied_cut):
            restricted = [[sum(a * b for a, b in zip(c, row)) for row in B] for c in C]
            # This equals the algebraic maximum min(q, dim sigma), so it is
            # also the exact rational rank.  A pivot minor is invertible in
            # Z_(Q31); null_mod is therefore a reduction of a rational basis.
            assert rank_mod(restricted, Q31) == min(len(C), len(B))
            coefficients = null_mod(restricted, len(B), Q31)
            out = [[sum(x * B[i][j] for i, x in enumerate(z)) % Q31
                    for j in range(h)] for z in coefficients]
            assert len(out) == max(len(B) - len(C), 0)
            assert rank_mod(out, Q31) == len(out)
            assert all(sum(x * y for x, y in zip(row, c)) % Q31 == 0
                       for row in out for c in C)
            return out

        truncated = [intersection(B) for B in old_frames]
        extra = dict(original_extra)
        extra.update({("sigma", i): B for i, B in enumerate(truncated)})
        mq = cf.ModQ7(fr, ld, extra)
        # A determinant modulo Q31 certifies exact nondegeneracy, not just
        # rank of a symbolic formal frame.  Check even low-rank steps' frames.
        for i in range(len(truncated)):
            assert mq.nondeg(("sigma", i)), (q, "degenerate", i)
        assert mq.dimfail == 0
        log(f"q={q}: all {len(truncated)} distinct intersections have the expected dimension and nonzero Gram determinant.")

        def sigma(s):
            i = slot_index[s]
            return ("sigma", i) if truncated[i] else ("0",)

        def dim(key):
            if key[0] == "sigma":
                return len(truncated[key[1]])
            return S.dim(key)

        edges = {}
        raw_counts = Counter()
        same_dimension = 0
        def edge(A, B, kind):
            nonlocal same_dimension
            r = dim(B) - dim(A)
            assert r >= 0, (q, kind, A, B)
            # Exact nesting is preserved by a common intersection.  An equal
            # dimension edge is thus equality and needs no recursive child.
            if r == 0:
                same_dimension += 1
                return
            raw_counts[kind] += 1
            if (A, B) in edges:
                assert edges[A, B][0] == r
                edges[A, B][1].add(kind)
            else:
                edges[A, B] = (r, {kind})

        for s in S.readout:
            edge(sigma(s), S.start_key(s), "first_internal")
            edge(sigma(s), ("F",), "exterior_corner")
        for ring, targets in by_target.items():
            for t, slots in targets.items():
                previous = ("0",)
                for s in slots:
                    here = sigma(s)
                    edge(previous, here, "Y_" + ring)
                    previous = here
                T = S.trip[t]
                edge(previous, ("out", T[0], T), "Y_" + ring)

        # Nondegenerate nested A <= B have P_B-P_A an idempotent of rank
        # dim B-dim A.  Complement/time reversal gives the same difference.
        used = {k for A, B in edges for k in (A, B)}
        assert all(mq.nondeg(k) for k in used)
        assert mq.dimfail == 0
        checks = Counter()
        high = [(A, B, r, kinds) for (A, B), (r, kinds) in edges.items() if 2 * r > h]
        high.sort(key=lambda x: (str(x[1]), str(x[0])))
        for A, B, r, kinds in high:
            QA, QB = mq.conjugated(A), mq.conjugated(B)
            for side in range(2):
                a, b = QA[side], QB[side]
                M = [[(b[i * h + j] - a[i * h + j]) % Q31 for j in range(h)] for i in range(h)]
                actual = rightmost_pivots(M, Q31)
                assert actual == expected[r], (q, A, B, r, side, actual, expected[r])
                checks[r] += 1
        assert mq.dimfail == 0
        log(f"q={q}: all {len(high)} distinct high-rank changed edges pass both common-basis NE profiles.")

        # Negative controls demonstrate that dimensions alone, arbitrary cuts,
        # and arbitrary coordinate bases are insufficient certificates.
        bad_cut = [cut[0], cut[0]]
        rank_failure = next(i for i, B in enumerate(old_frames) if len(B) >= 2 and
            rank_mod([[sum(a * b for a, b in zip(c, row)) for row in B] for c in bad_cut], Q31) < 2)
        independent_cut = [cut[0][:]]
        independent_cut[0][0] += 1
        broken_chain = None
        for slots in by_target["Z"].values():
            for a, b in zip(slots, slots[1:]):
                if slot_index[a] == slot_index[b] or len(S.sigma[a]) <= 1:
                    continue
                A = intersection(S.sigma[a], cut[:1])
                B = intersection(S.sigma[b], independent_cut)
                if rank_mod(A + B, Q31) > len(B):
                    broken_chain = [a, b]
                    break
            if broken_chain:
                break
        assert broken_chain is not None
        identity_basis_failures = 0
        control = cf.ModQ7(fr, ld, extra, control=True)
        for A, B, r, _ in high[:50]:
            if r == h:
                continue
            qa, qb = control.conjugated(A), control.conjugated(B)
            M = [[(qb[0][i*h+j] - qa[0][i*h+j]) % Q31 for j in range(h)] for i in range(h)]
            identity_basis_failures += rightmost_pivots(M, Q31) != expected[r]
        assert identity_basis_failures > 0
        f_hist = Counter(max(f - q, 0) for f in S.f)
        result = dict(
            codimension=q, deferred_slots=len(S.readout), distinct_old_frames=len(old_frames),
            dimension_formula="max(dim(sigma)-q,0)",
            all_rational_intersections_realized=True, all_changed_frames_nondegenerate=True,
            changed_edges_by_kind=dict(raw_counts), distinct_changed_edges=len(edges),
            omitted_equal_frame_edges=same_dimension,
            high_rank_changed_edges=len(high), common_basis_conjugations_per_edge=2,
            NE_checks_by_rank=dict(sorted(checks.items())), NE_failures=0,
            entrance_dimension_histogram=dict(sorted(f_hist.items())),
            largest_auxiliary_block=max(h*h-2*(h-f) for f in f_hist),
            negative_controls=dict(duplicate_cut_rows_detected_at_frame=rank_failure,
                                   per_readout_cuts_break_Y_chain=broken_chain,
                                   identity_basis_profile_failures=identity_basis_failures),
        )
        if q == 22:
            result["zero_entrance_cut"] = "H=span(e_0)"
            result["integer_cut_rows"] = applied_cut
            assert all(not B for B in truncated)
        results.append(result)

    out = dict(
        source_repository="https://github.com/Swapnil-jain/integer-mult-kappa",
        source_commit="741e7aa078392553815df7926ee17ac5e25a8c38",
        source_sha256=PINS, h=h, modulus=Q31, seed=args.seed,
        integer_cut_rows=cut, common_axis_basis="check_lifted.point_UV(23,1)",
        common_axis_basis_sha256=basis_hash, profiles=results,
        unchanged_geometry_dependency="Pinned check_lifted.py and check_frames.py; exact inclusions are inherited, never inferred from modular containment.",
        full_multiplication_bound=False,
        elapsed_seconds=round(time.monotonic() - start, 3),
    )
    if args.output:
        args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
