# Smaller deferred entrances, with realized rational frames

This package changes the deferred entrance frames of the **pinned round-seven
bit word**. It leaves its scalar operations, physical roles, lifted producer
frames, leaf-copy frames, retained centres, and data entrances unchanged.
Its targeted candidate caps only the dimension-22 deferred frames at
dimension 21. The largest auxiliary block falls from **527/529 to 525/529**,
while only 1,760 of the 28,866 auxiliary entrances change.

These are finite frame/profile certificates. They do not establish a new
complete multiplication theorem or a practical input-size crossover. The
separate [interface audit](../deferred-interface-target-10/) assembles the
histograms, moments, and row-stock tradeoffs. The physical role count does
not decrease; the smaller maximum child can reduce the recursion-depth
constant used to bound the required row stock.
The [coordinated-cap package](../deferred-cap-target-10/) extends this idea
to caps of 20 and 19.

## Attribution and reproducibility

The original word, deferred frames, common flag family, and staircase
entrance are from [Swapnil Jain's round-seven repository](https://github.com/Swapnil-jain/integer-mult-kappa/tree/741e7aa078392553815df7926ee17ac5e25a8c38),
specifically [deferred-readout.tex](https://github.com/Swapnil-jain/integer-mult-kappa/blob/741e7aa078392553815df7926ee17ac5e25a8c38/notes/deferred-readout.tex)
and [flag-basis.tex](https://github.com/Swapnil-jain/integer-mult-kappa/blob/741e7aa078392553815df7926ee17ac5e25a8c38/notes/flag-basis.tex).
The source was imported into [CrocSwap PR #97](https://github.com/CrocSwap/integer-mult-bounds/pull/97).
Zhihao Chen's PR #97 work supplies the gauged exits, reflected scheduling,
and signed integration used by the accompanying interface audit.
This package adds the two frame modifications below and their checks; it
does not claim the original deferred construction.

Pass the unpacked pinned `swapnil-round7` directory as `--source`. The
checkers verify SHA-256 hashes of both input certificates and all imported
linear-algebra/frame modules. They require only the Python standard library.

```sh
python3 check_cap21.py --source /path/to/swapnil-round7 --output cap21-results.json
python3 check_truncation.py --source /path/to/swapnil-round7 --output exact-results.json
python3 check_truncation.py --source /path/to/swapnil-round7 --codim 22 --output zero-results.json
```

The complete unchanged geometry is independently reproducible with the
source's `check_lifted.py 23` and `check_frames.py 23`; its scalar word and
chronological frame keys are audited in the sibling interface package.
Neither this package nor a dimension table substitutes for those checks.

## Targeted modification: cap only the top deferred frames

Write `T⊥` for the dimension-22 target hyperplane in `Q^23`, under
`G = I − J/9`. An old deferred entrance σ is contained in every target
hyperplane reached by its garbage readout. If `dim σ = 22`, all those
hyperplanes must equal σ. Distinct triples have distinct target hyperplanes,
so such a slot has exactly one target. The checker also verifies this
directly on the full integer garbage coefficients, which contain the F2
support.

For each affected target T, take its largest proper deferred frame M_T.
Every proper frame on that target's old nested chain lies in M_T. In the
actual pinned witness:

| Target count | dim M_T | New cap C_T |
|---:|---:|---|
| 1,485 | 21 | Reuse M_T exactly |
| 55 | 20 | Adjoin one explicit integer vector of T⊥ |

The added vector is chosen so `G|C_T` is nondegenerate. Every affected
dimension-22 slot for T uses **the same C_T**. All other σ remain unchanged.
The receipt records the 55 extra integer vectors and their source slots.
The checker verifies over the integers that every proper frame for T lies
in M_T, and that C_T lies in T⊥. Its Gram determinants are nonzero modulo
`2^31−1`, proving exact nondegeneracy over Q.

Thus all changed Y-chain edges are nested. For a changed auxiliary slot,
the old first frame F0 has dimension 22 and contains the old σ=T⊥; hence
F0=T⊥, and `C_T ⊂ F0`. Its first internal edge gains rank one. Its exterior
edge loses rank one. Their total rank is unchanged. The scalar readouts
remain at their original position in the word, now at C_T, and every gate
still uses the same frame on both registers.

Both the integer and F2 target supports have the same preceding maximum
dimensions listed above. All changed positive factor edges have rank one
or two: the first internal step, the gauged exterior corner, the entry to
C_T, and the final step to T⊥. The checker nevertheless evaluates **every
changed edge's actual north-east pivot pattern under both original fixed
axis bases U,V**. It does not use a different basis for each edge. The
changed exterior profile is `[525, 1, 1]`; the other changed steps split
into one or two singletons. Unchanged profiles retain their original
certificate.

## Uniform alternative: intersect every σ with one ambient cut

Choose four explicit integer covectors, recorded in `exact-results.json`.
For q=1,2,3,4 let H_q be the common kernel of their first q rows and set

`σ'_u = σ_u ∩ H_q` for every deferred slot u.

The same H_q is used on **all** slots and target chains. Therefore old
inclusions `σ_u ⊂ σ_v`, `σ_u ⊂ F0(u)`, and `σ_u ⊂ T⊥` imply the required
new inclusions. No producer or leaf-copy frame moves. There are 11,565
deferred slots but only 4,157 distinct supplied σ bases.

For each supplied integer basis B of σ, the restriction of the cut is
`C_q B^T`. Its rank modulo `p=2^31−1` equals `min(q, dim σ)`, an algebraic
upper bound on its rational rank. This proves

`dim σ' = max(dim σ − q, 0)` **exactly over Q**.

An invertible pivot minor modulo p is a unit in the localization `Z_(p)`.
The kernel basis calculated modulo p is consequently the reduction of a
rational kernel basis with denominators prime to p. Multiplication by B
gives the reduction of σ∩H_q. This justifies using the reduced bases to
certify nonzero rational Gram determinants and projector minors; a mere
modular containment test would not justify an exact rational inclusion.

| q | Maximum deferred dimension | Largest auxiliary child |
|---:|---:|---:|
| 0, original | 22 | 527 |
| 1 | 21 | 525 |
| 2 | 20 | 523 |
| 3 | 19 | 521 |
| 4 | 18 | 519 |

The checker covers every changed first internal edge, exterior corner,
and consecutive Y edge for **both integer and F2 coefficient supports**.
The first internal edge gains `dim σ−dim σ'`, while the exterior loses
the same amount, preserving the rank sum. The total Y-chain rank remains
22. The unchanged X chains and scalar word are inherited.

## A deterministic endpoint: every deferred entrance is zero

The separate `--codim 22` mode uses `H=span(e_0)`, defined by the 22
coordinate covectors `e_1^T,...,e_22^T`. Every deferred σ lies in at least
one target hyperplane with normal `9·1_T−3·1`. The normal's e_0 coordinate
is either 6 or −3, so that hyperplane meets H only at zero. Consequently
every σ∩H is zero, without a generic choice of H.

`zero-results.json` verifies the restriction ranks on all 4,157 supplied
frames and **10,548 distinct high-rank changed edges under both fixed axis
bases**, with zero failures. The largest auxiliary block is 483 instead
of 527. The original scalar chronology remains valid: the deferred slots
are untouched during the early producer phase, and their garbage reads
now occur at frame zero. The retained-centre scatters also occur at zero
before these reads, exactly as in the original order.

This endpoint trades the bit saving for smaller recursive children. It is
not a uniform improvement over earlier borrowed-profile options: with the
round-seven h=24 complex choice, the sibling arithmetic table's bit saving
and row-stock degree are both worse than the earlier PR84/h=18 combination.
The purpose here is to certify the finite endpoint, not claim that it is
the preferred global profile. The common-cut negative controls involving
different cuts illustrate the codimension-one nesting requirement; they
are not counterexamples to this deterministic zero-frame construction.

## Why these are actual common-basis profiles

For nested G-nondegenerate subspaces A⊂B, orthogonal projectors satisfy
`P_A P_B = P_B P_A = P_A`. Thus `P_B−P_A` is an idempotent of rank
`r=dim B−dim A`. Complementing and reversing time gives the same difference,
so both stages are covered by the two axis conjugations.

For any h×h rank-r idempotent P, a north-east submatrix with rows 0..i and
columns j..h−1 has rank at most

`min(i+1, h−j, r, max(0,i−j+1)+h−r)`.

The first three bounds are immediate; the last follows from
`P=I−(I−P)` and `rank(I−P)=h−r`. The checkers calculate the complete
rightmost-pivot pattern of `U^T(P_B−P_A)U^−T` and
`V^T(P_B−P_A)V^−T`. Equality with that rank function modulo p proves
equality over Q: the modular minors provide lower bounds matching the
universal upper bounds. Both use exactly
`check_lifted.point_UV(23,1)`, the source's common axis point.

For r>h/2, the resulting profile consists of h−r singletons and one block
of width 2r−h. For smaller r, splitting into r singletons is always allowed.
`check_truncation.py` checks the high-rank profiles; `check_cap21.py` also
checks the small changed profiles explicitly. The source common-flag
lemma transfers these axis profiles into the tensor corner. Its existing
line factors, staircase construction, and unchanged edges are preserved.

The uniform-cut checker includes negative controls: duplicate cut rows
break the asserted dimension formula; giving different cuts to adjacent
readouts breaks a Y-chain inclusion; and replacing U,V by the identity
breaks the claimed NE pattern. These distinguish an exact frame/profile
certificate from rank bookkeeping alone.

The new rational subspaces and projector coefficients are fixed finite
data. Their denominators can enlarge the finite excluded-prime list in
the inherited address geometry. This is an eventual finite construction,
with no claimed numerical machine runtime or tested finite crossover.
