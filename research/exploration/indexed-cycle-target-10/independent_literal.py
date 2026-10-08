"""Independent, import-free audit of PR84's serialized framed XOR words.

This checks literal scatter linkage as well as the scalar map and frame path.
Arbitrary-dirty correctness follows algebraically from the emitted word:
  M J M^-1 V M J M^-1 V,
where M acts only on auxiliary roles, J adds auxiliary controls to y, and V
adds x to distinct source roles. The two J contributions from initial dirty
auxiliaries cancel over F2, and the remaining contribution is S M P x.
Run: python3 independent_literal.py --source /path/to/pinned/pr84-checkout

Original independent verifier written with OpenAI Codex assistance. Apache-2.0.
All external construction credit belongs to the authors retained in
UPSTREAM-NOTICE. No external compiler or replay implementation is imported.
"""
from argparse import ArgumentParser
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import gzip
import json
import sys
from independent_moment import check_source


def check(root, h):
    raw = gzip.decompress((root / f'certificates/indexed-cycle-word-{h}.json.gz').read_bytes())
    d = json.loads(raw)
    triples = list(combinations(range(h), 3))
    masks = [sum(1 << i for i in t) for t in triples]
    v, roles = len(triples), d['R']
    assert (d['h'], d['v']) == (h, v)
    frames = [tuple(f) for f in d['frames']]
    lookup = {f: i for i, f in enumerate(frames)}
    assert len(lookup) == len(frames)
    rank = []
    for core, cover in frames:
        assert 0 < core <= cover < 1 << h and not core & ~cover
        assert core.bit_count() in (1, 2) or core == cover and core.bit_count() == 3
        rank.append(1 if core == cover else cover.bit_count() - core.bit_count())
    state = [None] * roles
    vector = [0] * roles
    events = []
    transitions = Counter()

    def arrive(slot, frame):
        assert 0 <= slot < roles and 0 <= frame < len(frames)
        prior = state[slot]
        if prior == frame:
            return
        core, cover = frames[frame]
        if prior is None:
            previous_rank = 0
        else:
            oldcore, oldcover = frames[prior]
            assert not core & ~oldcore and not oldcover & ~cover
            previous_rank = rank[prior]
        assert rank[frame] > previous_rank
        transitions[rank[frame] - previous_rank] += 1
        events.append([slot, -1 if prior is None else prior, frame])
        state[slot] = frame

    sources = {int(k): val for k, val in d['sources'].items()}
    assert set(sources) == set(range(v)) and len(set(sources.values())) == v
    for i, slot in sources.items():
        arrive(slot, lookup[masks[i], masks[i]])
        vector[slot] = 1 << i
    for dest, source, frame in d['ops']:
        assert dest != source
        arrive(dest, frame)
        arrive(source, frame)
        vector[dest] ^= vector[source]
    expected_scatter = Counter()
    scalar_rows = [0] * v
    outputs = set()
    retained = 0
    for slot, frame, common, target in d['outputs']:
        assert slot not in outputs
        outputs.add(slot)
        arrive(slot, frame)
        assert common in target and len(target) in (1, 3)
        exclude = set(target) - {common}
        desired = sum(1 << i for i, t in enumerate(triples)
                      if common in t and not exclude.intersection(t))
        assert vector[slot] == desired
        expected_frame = (1 << common, ((1 << h) - 1) ^ sum(1 << j for j in exclude))
        assert frames[frame] == expected_frame
        if len(target) == 1:
            retained += 1
            targets = [i for i, t in enumerate(triples) if common in t]
            transitions[rank[frame]] += 1  # copied retained centre
            transitions[h - rank[frame]] += 1  # original cleanup
        else:
            targets = [triples.index(tuple(target))]
            transitions[1] += h - rank[frame]  # side growth split into singletons
        for i in targets:
            expected_scatter[v + i, 2 * v + slot] += 1
            scalar_rows[i] ^= vector[slot]
    assert retained == h and len(outputs) == h + 3 * v
    assert Counter(map(tuple, d['scatter'])) == expected_scatter
    assert scalar_rows == [1 << i for i in range(v)]
    # The compiler may preallocate an unused slot before its first operand
    # appearance and may record idempotent incidences. Verify its per-slot
    # chronology, then compare all nontrivial paid moves as a multiset.
    event_state = [-1] * roles
    for slot, prior, frame in d['events']:
        assert 0 <= slot < roles and event_state[slot] == prior
        assert 0 <= frame < len(frames)
        event_state[slot] = frame
    assert event_state == state
    assert Counter(map(tuple, events)) == Counter(tuple(e) for e in d['events'] if e[1] != e[2])
    for slot, frame in enumerate(state):
        assert frame is not None
        if slot not in outputs:
            transitions[h - rank[frame]] += 1
    mass = sum(k * n for k, n in transitions.items())
    assert mass == h * roles + h * (h - 1)
    return dict(h=h, source_sha256=sha256(raw).hexdigest(), roles=roles,
                unique_output_roles=len(outputs),
                literal_xors=len(d['ops']), retained_centres=retained,
                scatter_xors=len(d['scatter']), actual_events=len(events),
                rank_mass=mass, scalar_identity=True, literal_scatter_linked=True,
                arbitrary_dirty_identity='M J M^-1 V M J M^-1 V over F2')


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    root = parser.parse_args().source.resolve()
    check_source(root)
    print(json.dumps([check(root, h) for h in (23, 25)], indent=2))
