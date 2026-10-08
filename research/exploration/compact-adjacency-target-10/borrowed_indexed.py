"""Audit a borrowed-source wrapper around the pinned, credited PR84 producer.

This is an exploratory finite network certificate, not a complete multiplication
theorem. The producer, its frames and the old profile are Chafik Boukhalfa's;
the new dirty-cancellation wrapper and its independent replay are described in
borrowed-sources.tex. See ../indexed-cycle-target-10/UPSTREAM-NOTICE.
No external compiler is imported. Run with --source /path/to/pinned/PR84.
"""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import gzip
import json
import sys

if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'indexed-cycle-target-10'))
import independent_literal as literal
import independent_moment as arithmetic


def producer(root, h):
    original = literal.check(root, h)
    raw = gzip.decompress((root / f'certificates/indexed-cycle-word-{h}.json.gz').read_bytes())
    d = json.loads(raw)
    triples = list(combinations(range(h), 3))
    v, old_r = len(triples), d['R']
    sources = {int(k): val for k, val in d['sources'].items()}
    source_slots = set(sources.values())
    assert len(source_slots) == v
    rename = {slot: i for i, slot in sources.items()}
    rename.update({slot: 2*v+i for i, slot in enumerate(s for s in range(old_r) if s not in source_slots)})
    roles = old_r + v
    assert set(rename.values()) == set(range(v)) | set(range(2*v, roles))
    middle = [(rename[a], rename[b]) for a, b, _ in d['ops']]
    centre_slots = {slot for slot, frame, common, target in d['outputs'] if len(target) == 1}
    centre = [(a, rename[b-2*v]) for a, b in d['scatter'] if b-2*v in centre_slots]
    side = [(a, rename[b-2*v]) for a, b in d['scatter'] if b-2*v not in centre_slots]
    assert len(centre_slots) == h and all(v <= a < 2*v for a, _ in centre+side)

    # Every integer is a complete row on ALL original inputs, including dirt.
    initial = [1 << i for i in range(roles)]
    state = initial[:]
    for a, b in middle:
        state[a] ^= state[b]
    output = [0]*v
    for a, b in centre+side:
        output[a-v] ^= state[b]
    xmask = (1 << v)-1
    ymask = xmask << v
    assert [r & xmask for r in output] == [1 << i for i in range(v)]
    assert not any(r & ymask for r in output)
    dirt = [r & ~((1 << (2*v))-1) for r in output]

    # Forward chronological word: early D, M, retained J, side J, M^{-1}.
    state = initial[:]
    for i, row in enumerate(dirt):
        state[v+i] ^= row
    for a, b in middle+centre+side+middle[::-1]:
        state[a] ^= state[b]
    expected = initial[:]
    for i in range(v):
        expected[v+i] ^= initial[i]
    assert state == expected

    # Stage two uses the SAME chronology and inverse-transposes each gate.
    # All CNOTs are involutions. Aggregate only the simultaneous early fanout.
    state = initial[:]
    early_terms = 0
    for i, row in enumerate(dirt):
        while row:
            bit = row & -row
            state[bit.bit_length()-1] ^= initial[v+i]
            row ^= bit
            early_terms += 1
    for a, b in middle+centre+side+middle[::-1]:
        state[b] ^= state[a]
    expected = initial[:]
    for i in range(v):
        expected[i] ^= initial[v+i]
    assert state == expected

    # A specific omitted dirty correction must change the forward scalar map.
    assert any(dirt)
    victim = next(i for i, row in enumerate(dirt) if row)
    deleted = dirt[victim] & -dirt[victim]
    state = initial[:]
    for i, row in enumerate(dirt):
        state[v+i] ^= row ^ (deleted if i == victim else 0)
    for a, b in middle+centre+side+middle[::-1]:
        state[a] ^= state[b]
    assert state[v+victim] == initial[v+victim] ^ initial[victim] ^ deleted

    # Initial frame entrances are already carried by the borrowed X lines.
    # Every subsequent common gate/frame and output frame is unchanged.
    source_events = Counter()
    for slot, prior, frame in d['events']:
        if prior == -1 and slot in source_slots:
            core, cover = d['frames'][frame]
            assert core == cover and core.bit_count() == 3
            source_events[slot] += 1
    assert source_events == Counter({s: 1 for s in source_slots})
    old_profile = json.loads((root / f'certificates/indexed-cycle-profiles-{h}.json').read_text())
    local = Counter({r: n for r, n in enumerate(old_profile['blocks']) if r and n})
    local[1] -= v
    assert local[1] >= 0
    assert sum(r*n for r, n in local.items()) == h*(old_r-v) + (h-1)*v + h*(h-1)
    digest = sha256()
    width = (roles+7)//8
    for row in dirt:
        digest.update(row.to_bytes(width, 'little'))
    return dict(h=h, v=v, old_aux_roles=old_r, new_aux_roles=old_r-v,
                borrowed_source_roles=v, all_input_dimension=roles,
                dirty_cancellation_nonzeros=early_terms,
                dirty_matrix_sha256=digest.hexdigest(),
                producer_word_sha256=sha256(raw).hexdigest(),
                full_dirty_basis_forward=True,
                full_dirty_basis_same_chronology_inverse_transpose=True,
                omitted_dirty_gate_detected=True,
                inherited_frame_events_checked=original['actual_events'],
                local_blocks={str(r): n for r, n in sorted(local.items()) if n})


def profile():
    old = arithmetic.profile()
    rows = Counter(old['rows'])
    n, m = 4073300, old['m']
    removed = Counter()
    for h in (23, 25):
        # Per invocation/source: external selected pair, leaf entrance,
        # and the superseded source data's final (h-1)-dimensional raise.
        local = Counter([m-2*h, h, 1, 1, h-2])
        assert sum(r*c for r, c in local.items()) == m
        for r, count in local.items():
            removed[r] += n*count
            rows[r] -= n*count
    assert all(c >= 0 for c in rows.values())
    rows = Counter({r: c for r, c in rows.items() if c})
    w = old['W']-2*n
    mass = sum(r*c for r, c in rows.items())
    assert mass == old['mass']-2*n*m
    assert w*m-mass == old['deficit']
    p = dict(m=m, W=w, rows=rows)
    # Deliberately leave numerical room; this is not an optimized exponent.
    saving = Q(5488, 10**8)
    low, high = arithmetic.moment(p, saving)
    assert high < 1
    target_low, _ = arithmetic.moment(p, Q(1, 1023))
    assert target_low > 1
    return dict(m=m, W=w, total_rank=mass, rank_deficit=w*m-mass,
                old_W=old['W'], removed_roles=2*n,
                removed_children={str(r): c for r, c in sorted(removed.items())},
                child_multiplicities={str(r): c for r, c in sorted(rows.items())},
                bit_saving=str(saving), exact_moment_gap=str(1-high),
                target_10_bit_moment_excess_lower=str(target_low-1),
                scope='Finite new-wrapper profile. Conditional on documented lifted frames and inherited global interfaces; no new complete multiplication certificate.')


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=HERE/'borrowed-indexed-audit.json')
    args = parser.parse_args()
    commit = arithmetic.check_source(args.source)
    result = dict(source_commit=commit,
                  producers=[producer(args.source, h) for h in (23, 25)],
                  profile=profile())
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result['profile'].items()
                      if k not in ('child_multiplicities', 'removed_children')}, indent=2))
