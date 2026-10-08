# Audited external profile and exploratory lower-overhead combinations

**This directory does not establish κ > 2⁻¹⁰, and does not change this
repository's main result.** It independently checks the finite profile from
[Chafik Boukhalfa's PR84](https://github.com/CrocSwap/integer-mult-bounds/pull/84)
and records candidate combinations with this repository's smaller complex
circuit and precision guard. The **combined global interface remains
unreviewed**. Every candidate exponent below is smaller than PR84's credited
headline κ = 0.0000522712890914.

The external source is pinned to
[`88ca39571907343a49e97f328971ec7bcd26fbfd`](https://github.com/chafreaky/integer-mult-bounds/tree/88ca39571907343a49e97f328971ec7bcd26fbfd).
Its [proof and scope](https://github.com/chafreaky/integer-mult-bounds/blob/88ca39571907343a49e97f328971ec7bcd26fbfd/research/indexed-cycle/PROOF.md)
describe reclaimed dirty roles, pending controls, carry exchanges and
two-cycles. Those changes are external contributions, not discoveries of this
repository. The changed role formula requires checking a ceiling's hypotheses
before transferring it. In particular, these serialized words still have
3v+h distinct output slots, loss h(h−1), and the relevant two-part exterior
rank profile. They therefore satisfy the hypotheses of the separately checked
rectangular ceiling in the sibling `rectangular-target-10` directory; this
compiler does not evade that scoped barrier to the requested target.

## What was checked

The external focused verification completed successfully: it regenerated the
h23 and h25 physical words byte-for-byte, replayed every ordinary and arbitrary
dirty basis column in both orientations, reconstructed every charged frame
transition, regenerated the fixed-basis CRT profiles, ran five compiler audits,
and passed 23 mutation/unit tests. Its source closure contains 580 pinned
files. The saved log is [results/upstream-focused.log](results/upstream-focused.log).
This is execution of the external authors' verification package, not a claim
that its implementation was independently rederived.

The additional local checks import no external compiler or arithmetic code:

- `independent_literal.py` reconstructs the scalar map from the serialized XOR
  operands, checks all frame incidences, verifies every scalar output, links
  the literal scatter list to those outputs, and checks all nontrivial charged
  frame moves. The identity `M J M⁻¹ V M J M⁻¹ V` restores arbitrary dirty
  auxiliaries: the two contributions from the initial dirty state cancel over
  F₂, leaving precisely the independently checked input-to-output map.
- `independent_moment.py` reconstructs all 27 child multiplicities from the
  two tiny axis profiles. It obtains m = 575, W = 132,466,108, total rank
  76,166,165,200 and deficit 1,846,900. Independent rational logarithm and
  exponential enclosures verify the published bit saving
  408390793141/7812500000000000. The acceptance gap is at least
  3.3380659867 × 10⁻¹⁹; increasing the saving by 10⁻¹⁸ fails.
- `combined_assembly.py --check-complex` recomputes this repository's h18
  complex histogram, all label incidences, all 665,856 ordered scalar pairs
  (including a flipped-piece negative control), and scalar-prefix guard, then checks
  the candidate macro-stack inequalities exactly. These inequalities are the
  inherited enlarged-precision stack, **not** PR84's refined 47-row assembly.

The separate [conditional transport audit](interface-audit.md) checks the fixed
word, fresh row-stock constants and smaller precision exponent against those
interfaces. It found no additional exponent restriction; it does not replace
the inherited global proofs.

Passing these finite checks does not establish the all-size ordered residual
compiler, routing, prime selection, exact recovery, fixed-tape simulation or
analytic transfer. Those remain inherited dependencies. The new combination
also needs its global interface checked before it can be promoted.

## Candidate tradeoffs

The h18 complex saving is 0.000059292472. With β = 1/16, its leaf saving
exceeds the external bit saving. This reduces the complex role volume from
PR84's 537,696,432 to 32,011,680. The same product-stock arithmetic then has
coefficient 568 instead of 843; row degree 1200 has gap 1032/25 > 0 and suffix
slope 4800. This arithmetic is recorded as a candidate interface, not a
completed integration proof.

| Candidate | κ | ε | Precision power P | Axis-width power 1−ε |
|---|---:|---:|---:|---:|
| Larger saving | 0.00005227 | 0.99994 | 3 | 0.00006 |
| Wider analytic slack | 0.000037 | 0.71 | 3 | 0.29 |
| Wider slack with bit-recurrence backoff | 0.000037 | 0.71 | 3 | 0.29 |
| Lower precision with wider slack and backoff | 0.000037 | 0.71 | 2.15 | 0.29 |

The last two candidates deliberately use the smaller bit saving
522/10⁷ = 0.0000522. Its separately certified moment gap is about
3.43 × 10⁻⁸, instead of the published profile's near-root gap of about
3.34 × 10⁻¹⁹. They use δ = 0.09 and two recurrence-exponent gaps of 10⁻⁸.
Their butterfly margin is exactly 0.0000370478, exceeding κ by 4.78 × 10⁻⁸.
The last candidate's fine-field margin is 193/2000 = 0.0965.
These wider inequalities may reduce eventual thresholds; they do **not** give
a finite operational crossover or a practical multiplication runtime.

## Rounded fractional precision

The pinned stack already allows a fixed rational precision exponent:
[stack-notes.tex, Idea D](https://github.com/Swapnil-jain/integer-mult-kappa/blob/f2176bc1124821bf17eb63725bd366d7bdc020a3/notes/stack-notes.tex).
For the last candidate set P = 43/20, ε = 71/100, and

```
b' = 128 ceil(b^P),    p = 6b',    d = floor(b^ε),
α² = floor(b'/(8d)),   θ ≥ 1/(8d L^ε),   1 ≤ L ≤ b.
```

Then P − 3ε = 1/50, and the floor loss gives

```
α² θ ≥ 2 b^(P−3ε) − 1/8 ≥ (15/8) b^(1/50) > 1.
```

The same linear guard is covered because
128(d+1) ≤ 256b^ε ≤ 768b^P ≤ p for b ≥ 1. The exact upper bound
2dα² ≤ b'/4 preserves the inherited final-rounding condition. Also
128b^P ≤ b' ≤ 256b^P, so logical volume Tp remains Θ(n), and the ceiling
does not change any recorded exponent. The rational ceiling is computable
with integer arithmetic, as the smallest integer z with z²⁰ ≥ b⁴³.
The axis width remains Θ(b^0.29); the record exponent can be μ = 29/430 > 0.
These statements check precision and exponent arithmetic; they do not supply
the omitted global routing and setup proofs.

## Reproduce

From this directory, the bundled numerical inputs suffice for:

```sh
python3 independent_moment.py
python3 combined_assembly.py --check-complex
```

To check the physical words, obtain the exact external commit separately.
The approximately 132 MB archive and large word files are intentionally not
duplicated here:

```sh
git clone --filter=blob:none --no-checkout https://github.com/chafreaky/integer-mult-bounds.git /tmp/pr84-pinned
git -C /tmp/pr84-pinned checkout --detach 88ca39571907343a49e97f328971ec7bcd26fbfd
python3 independent_moment.py --source /tmp/pr84-pinned
python3 independent_literal.py --source /tmp/pr84-pinned
make -C /tmp/pr84-pinned indexed-cycle-verify
```

The `--source` checks bind the selected word/profile files and independently
hash all 580 entries of the pinned source manifest. Python assertions must
remain enabled. The saved JSON results describe the actual runs, not a promise
that any arbitrary source checkout is equivalent.

## Attribution and files

The tiny numerical profiles are copied from the pinned external source; the
complete-child summary is an explicitly identified extraction of its
certificate. Original Apache-2.0 and AI-assistance notices are preserved in
[UPSTREAM-NOTICE](UPSTREAM-NOTICE), including Chafik Boukhalfa, Thomas DiFiore,
Rohan Arun, Eumemic, Rohan Garg, Avi Eisenberg, Dominik Scholz, Alejandro
Zarzuelo Urdiales and the earlier contributor chain. All external work should
be credited through PR84 and those notices. The repository root retains the
Apache-2.0 license. The new audits and exploratory parameter combinations were
prepared with OpenAI Codex assistance at the repository owner's request.

`inputs/source.json` records URLs, commit and source hashes. `inputs/` has only
the two axis profiles and the small extracted complete bit profile.
`results/` contains the independent replay, moment, candidate-assembly and
external focused-test receipts. No main theorem, headline certificate or
existing note is changed by this directory.
