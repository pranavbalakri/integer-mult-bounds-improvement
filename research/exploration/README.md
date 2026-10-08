# Search toward kappa > 2^-10

**Status: the target has not been attained.** The main certificate and LaTeX
proof remain the previously recorded conditional result, kappa=0.000036942.
The files here include a new finite network wrapper, candidate assembly
inequalities, and explicitly scoped obstructions. Passing them does not prove
the complete inherited multiplication theorem.

## New source-role borrowing candidate

The [borrowed-source construction](compact-adjacency-target-10/README.md)
uses existing source data rows in place of separate producer input roles.
A fixed early correction cancels every arbitrary initial auxiliary value.
Retained-centre reads precede side reads, and the second stage uses the same
chronology with inverse-transposed scalar gates. This preserves the checked
frame paths. Fresh temporary storage can implement the early correction using
the existing sparse producer, without a dense-matrix expansion.

Applied to the separately credited, pinned PR84 producer below, this removes
8,146,600 bit-network roles: W falls from 132,466,108 to 124,319,508, about
6.15%. Complete scalar-basis replays at h23 and h25 include all original dirty
auxiliary inputs, both stage orientations, and a failing negative control.
Two independently constructed recursive histograms agree. Their exact moment
supports the deliberately backed-off bit saving a_b=0.00005488.

Combining that profile with the existing h18 complex word and guard gives
two [candidate assembly choices](compact-adjacency-target-10/borrowed-assembly-results.json):

| Choice | kappa | epsilon | Coefficient precision, b=ceil(log2 n) |
|---|---:|---:|---|
| Larger saving | 0.00005487 | 0.9999 | 768 b^3 |
| Lower precision and wider margins | 0.0000389 | 0.71 | 768 ceil(b^(43/20)) |

All displayed assembly inequalities are exact and strict. The fixed-tape,
ordered residual, common-basis, routing, prime and recovery interfaces remain
inherited assumptions. These candidates have not replaced the main note or
certificate. Neither reaches 1/1024, and neither is a measured speedup.
The higher candidate is about 4.97% above the pinned PR84 headline; this is an
extension of that credited construction, not an independent discovery of its
underlying producer. No claim is made about priority over later public work.

## What the search established

- [Rectangular obstruction](rectangular-target-10/README.md): if the two-stage
  network retains the specified centre losses, endpoint correction and at least
  three side-output roles per triple, its bit saving is strictly less than
  1/1300, even with ideal batching and unequal factor dimensions. The proof
  combines 9,243 exact rational checks with an analytic bound for every larger
  dimension. This rules out a whole parameter search, not other architectures.
- [Complex obstruction](complex-target-10/architecture-obstruction.md): the
  current complex data/centre/endpoint profiles force a_c<1/1028. A direct
  one-stage shortcut loses its apparent advantage when its endpoint corrections
  are charged.
- [Central-factorization obstruction](complex-target-10/central-factorization.md):
  another minimal basis of the existing central matrix saves at most one unit
  of source-label dimension. Larger improvements require changing additional
  assumptions.
- [Alternative bit fitting matrix](bit-target-10/README.md): a scalar rank-six
  candidate exists at h=7, but the exact rational transfer ranks still give a
  negative recursion deficit. It is a rejected candidate, not an improved bound.
- [Centre rectangles](centre-rectangles-target-10/README.md): some nonzero read
  frames save one local rank, but the checked point-star/exclusion dictionary
  cannot assemble a positive h14 witness. The same directory independently
  checks the new borrowing identity over complex values and its linear guard;
  that complex variant has not been used in the candidates above.
- [New motif checks](new-motif-target-10/README.md): the E8 construction has
  complementary ranks 8 and 9, but its minimal centre factorizations are too
  expensive. A proposed sparse-inner-product Leech subset has at most 425
  lines. These are scoped rejections, not impossibility results for all motifs.

The desired saving is 1/1024=0.0009765625. Under the retained assembly,
achieving it requires a bit saving greater than 1/1023. Improving the current
producer alone is insufficient under the assumptions above.

## A separately attributed finite construction and an overhead candidate

[Chafik Boukhalfa's PR84](https://github.com/CrocSwap/integer-mult-bounds/pull/84),
pinned to `88ca39571907343a49e97f328971ec7bcd26fbfd`, supplies a newer finite
bit network. Its published bit saving is approximately 0.000052274021522048.
We independently reproduced its finite word/profile checks and exact moment;
this is verification of that author's construction, not a new construction
or a claim of priority here.

The [indexed-cycle audit](indexed-cycle-target-10/README.md) also records a
proposed combination with this repository's h=18 complex instance and linear
precision guard. A conservative choice uses bit saving 0.0000522,
epsilon=0.71, delta=0.09, and kappa=0.000037. Its coefficient precision is

    p = 768 * ceil(b^(43/20)),  b=ceil(log2(n)),

rather than 768*b^3, and its axis-width exponent is 0.29. Exact parameter and
precision checks pass. An [internal transport audit](indexed-cycle-target-10/interface-audit.md)
found no additional exponent restriction. The combined global transfer still
depends on the listed inherited contracts; finite checks do not prove them.
This candidate has not replaced the main certificate. It addresses analytic
overhead and gives no measured runtime or practical crossover.

## Reproduce the independent local checks

All of these use only the Python standard library:

```sh
python3 research/exploration/rectangular-target-10/check_barrier.py
python3 research/exploration/complex-target-10/check_obstruction.py
python3 research/exploration/complex-target-10/check_central.py
python3 research/exploration/bit-target-10/check_fitting.py
python3 research/exploration/indexed-cycle-target-10/independent_moment.py
python3 research/exploration/indexed-cycle-target-10/combined_assembly.py
python3 research/exploration/compact-adjacency-target-10/borrowed_sources.py --optimized 23
python3 research/exploration/compact-adjacency-target-10/check_borrowed_rect.py
python3 research/exploration/compact-adjacency-target-10/borrowed_assembly.py
python3 research/exploration/centre-rectangles-target-10/check_rectangles.py --dictionary
python3 research/exploration/centre-rectangles-target-10/inplace_complex_audit.py 8
python3 research/exploration/centre-rectangles-target-10/inplace_complex_guard.py
python3 research/exploration/new-motif-target-10/check_e8.py
```

The independent literal-word replay needs the separately obtained pinned PR84
source archive; its directory documents the exact source and command. The
optional SAT search also needs Z3, but the saved fitting witness can be checked
without it. No search heuristic or floating-point diagnostic is a certificate
of the target.

## Remaining research questions

The obstructions leave open different fitting matrices, redundant central
factorizations with smaller transfer cost, nonzero intermediate read frames,
different endpoint transfers, and entirely different finite motifs. Any such
candidate must pay for its full scalar word, frames, recursive children and
assembly before an exponent improvement can be claimed.
