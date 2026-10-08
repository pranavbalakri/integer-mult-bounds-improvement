# Conditional transport audit: PR84 bit word with the lower-overhead stack

This is an integration audit, not a new multiplication theorem or an independent
proof of the inherited machine and analytic lemmas. It does not modify a certified
file. The selected external word is the pinned PR84 indexed-cycle construction
at commit `88ca39571907343a49e97f328971ec7bcd26fbfd`, using the constants in
`research/indexed-cycle/PROOF.md` and `scripts/experiments/indexed_cycle_compose.py`.
Some top-level historical README constants describe other constructions.

## Verdict

I found no additional exponent restriction from transporting this fixed bit word
to the proposed parameters. In particular, neither the bit transfer contract nor
the cubic-stack derivation requires `x>=2`. That previous choice supplied slack
near `epsilon=1`; the actual precision inequality is `1+x>3*epsilon` when
`y=epsilon`. The proposed smaller precision exponent is admissible under the
same inherited interfaces. The main proof boundary remains the inherited global
transfer and analytic assumptions explicitly acknowledged by PR84.

## The finite bit word fits the general recursion interface

The independently checked finite profile has `m=575`, `W=132466108`, maximum child
width 529, and complete paid children including copied centres, exterior edges,
data growth, and the endpoint correction. Each child width is positive and less
than 575. The strict moment at the backed-off saving `a_b=522/10^7` is stronger
than the moment needed for the target stack.

The generic recurrence in `notes/batched-bit-rows.tex` accepts any such fixed
finite child multiset, with each child acting on logical volume `V/W` and width
`t*floor(e/m)`. Its proof covers integer remainders, bounded outer row padding,
complete spectator fields, arbitrary scratch restoration, and fixed-tape
parking. It does not require that the producer be the older leaf-addition DAG,
that both tensor factors have the same dimension, or that roles never be reused.
The audited physical word must provide those finite invariants; changing the
compiler search heuristic does not change this interface.

All region synthesis, carry matching, coordinate selection, profile computation,
and rational matrix preparation concern the fixed dimensions 23 and 25. Their
costs belong to fixed machine preparation, not to an input-size-dependent
arithmetic row. Runtime gates are the emitted finite XOR/copy/clear word and the
paid ordered-affine address transformations. Thus the discovery compiler creates
no new `p^delta` exponent requirement. Its fixed constants can still be enormous.

## Fresh row-stock bridge

The old PR84 bridge must not be copied unchanged: its complex component is a
different circuit. For the proposed bit and complex components, exact integer
comparisons give

| Component | m | Largest child | W | Halving degree t | ceil(log2 W) |
|---|---:|---:|---:|---:|---:|
| PR84 bit | 575 | 529 | 132466108 | 9 | 27 |
| Audited h18 complex | 324 | 306 | 32011680 | 13 | 25 |

Here `m^t>2*r^t`, and the preceding integer exponent fails that test. Simultaneous
row stock includes both factors, since a complex node may invoke a bit adapter.
Therefore

```
K = 9*27 + 13*25 = 568,
1200 - (51/25)*568 = 1032/25 > 0.
```

Row degree 1200 is sufficient for the inherited bound on calls of width `e<=Cp`.
This is a polynomial stock of *existing row coordinates*, not permission to
multiply the logical data volume by `p^1200`. Outer complete-row padding remains
bounded by a factor below two. The new certificate needs to record the fresh
complex geometry and new phase guard as well as these row constants.

For `p=768*ceil(b^(43/20))`, `b>=1`,
`log2(p)<=43/20*log2(b)+11<4*(log2(b)+8)`. Hence the conservative suffix condition
`b^(1-epsilon)>4800*(log2(b)+8)` still implies enough row stock. It holds eventually
for `epsilon=.71`; its existence is not a practical crossover calculation.

## Parameter transport

Take

```
a_b = 522/10^7, a_c = 59292472/10^12,
epsilon = 71/100, beta = 1/16, delta = 9/100,
x = 23/20, P = 1+x = 43/20, y = epsilon,
lambda = 1-a_b+10^-8, lambda_prime = 1-a_b+2*10^-8,
kappa = 37/10^6.
```

The complex leaf saving `(1-beta)*a_c=22234677/400000000000` exceeds `a_b`.
The precision gap is `P-3*epsilon=1/50`; the linear-guard exponent gap is
`P-epsilon=36/25`; and `0<delta=.09<1/8`. The principal margins are

```
prefix:                    29/100
butterflies:               185239/5000000000
fine/Gaussian/chirp/axes:   193/2000
CRT reversal:              18531/500000000
packed products:           71/100.
```

The butterfly surplus over kappa is `239/5000000000>0`. The reservation
inequality has surplus `1449814761/3550000000>0`. The finite bit compiler does
not involve delta; its use at `.09` comes from the unchanged arithmetic lemma
allowing any fixed delta in `(0,1/8)`.

Digits `b'=128*ceil(b^P)` and `p=6b'` preserve `Tp=Theta(n)` and the original
coefficient-capacity and final-rounding argument. With `d=floor(b^epsilon)`,
`L<=b`, `eta=1/(4d L^epsilon)`, and the inherited prime interval giving
`theta>eta`, the Gaussian choice `alpha^2=floor(b'/(8d))` has

```
alpha^2*theta >= 4*b^(P-3*epsilon) - 1/(4*d*L^epsilon) > 1
```

in the nontrivial working regime. Thus the small positive precision gap need
not introduce a separate huge Gaussian-width threshold. For the relaxed record
condition one may fix `mu=29/430`, strictly below `(1-epsilon)/P=29/215`.

## What remains assumed, and what must not be claimed

PR84 expressly inherits the ordered-affine residual implementation, all-size
recursion transfer, finite scalar alphabet/overhead, routing, prime selection,
exact recovery, fixed-tape simulation, and analytic transfer. Its finite tests
do not independently prove those statements, and neither does this audit.
The proposed composition additionally uses the separately checked h18 complex
word and the local linear phase guard; PR84's old complex precision certificate
cannot substitute for them.

The global proofs still contain eventual record, prime, row-stock, descriptor,
small-width, and arithmetic thresholds. The improved precision and delta
parameters reduce explicit sources of asymptotic overhead, but do not establish
seconds, a feasible memory footprint, or an end-to-end crossover. In particular,
the fixed word and its tape constants remain very large. No further integration
failure was found beyond these stated inherited proof boundaries and the need
to regenerate the row-stock bridge for the actual paired components.
