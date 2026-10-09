# Fewer references for sink unmixing

This package supplies exact scalar words and finite cleanup schedules that
reduce the reference bank of the [sink-unmixing prototype](../sink-unmixing-target-10/README.md).
In particular, the reference bank can be empty: only selected producer rows
clean early, and all other rows clean using the unchanged source data at the
full local frame. The arity-23 hierarchical variant passes a literal replay
on all 43,619 public scalar basis inputs, including arbitrary producer dirt.

These are **not new multiplication-exponent certificates**. The counts below
measure feasible scalar/frame incidences. A separately recorded full rational
profile of one zero-reference hierarchy on the pinned PR84 producer is
**worse than the existing borrowed producer**. Its source-data detours and
merged auxiliary exits are all charged. The smaller public width of the
borrowed producer is part of that comparison.

## Exact word with a selected reference bank

Use the unborrowed producer from `rtgm.build(h, retain='search')`, compiled
literally by the earlier scalar prototype. Its invertible auxiliary word is
`M`, source injection is `P`, and centre-plus-side scatter is `J`. Write
`T=MP`. The complete clean identity `JMP=I` is checked on every source column.
The data bank `X` is never overwritten in the forward word.

Choose a primary point `c(i)` in each source triple `S_i`. Choose a set of
reference pairs `(c,i)` with `c in S_i` and `c != c(i)`. Such a reference is an
arbitrary public input `B_(c,i)`; it receives `X_i` at the source-line frame.
The unchanged `X_i` itself supplies reads at its primary point, when needed.

A non-side-output producer row `s` may clean at the point-star frame `U_c`
only if both of the following hold:

1. Its **entire last producer frame** lies in `U_c`. The checker enforces the
   sufficient exact condition that every source triple in its last-node label
   contains `c`.
2. Every source `i` occurring in the actual clean coefficient row `T_s` has
   either `c(i)=c` or an allocated reference `(c,i)`.

The clean signal can be a proper subset of the last-node label, because an
operand role may retain its old scalar value after advancing to its consumer's
frame. Condition 1 cannot be replaced with the smaller signal support.

All unselected rows, including side outputs, use `X` directly at the full
local frame `D1`. **Late cleanup requires no references.** If `K_B` records
the selected reference reads and `K_X` the direct reads, including late reads,
then the exact coefficient identity is

```
T = K_X + K_B P_B.
```

The chronological word is:

1. At `D0`, apply `Y += JM Z` to the original producer dirt `Z`.
2. Create the selected reference signals; inject `PX` into `Z`; execute `M`.
3. Scatter retained centres while `Y` is at `D0`, then perform side scatters
   at their inherited target frames.
4. Clean selected rows with their available references/direct sources, then
   clean all remaining rows with `X` at `D1`.
5. At `D1`, erase the reference signals with `B += P_B X`.
6. At the common auxiliary sink, apply `Z += K_B B`, then `M^-1`.

After step 3, `Y=Y_original+X`. After step 5, the auxiliary block is exactly

```
(Z,B) = (M Z_original + K_B B_original, B_original).
```

Step 6 therefore restores every arbitrary auxiliary input. Its input offsets
are all zero and its output frames are the same global sink; the pointwise
inverse is consequently permitted there. The complete public scalar map is

```
(X,Y,Z,B) -> (X,Y+X,Z,B).
```

Transposing every elementary CNOT in the same chronological order gives the
opposite shear, because each gate and the public shear are involutions. Every
gate keeps the same two common-frame incidences. This scalar transpose check
does not remove the need to account for the second tensor stage separately;
the present early-retirement proposal is for the first stage.

If the reference set is empty, `K_B` and the reference creation, erasure, and
dirt-unmixing blocks all disappear. This is an actual arbitrary-dirt word,
not a calculation that assumes zero initial auxiliary values.

## Zero references and sparse references

`reference_plans.py` computes every row of `T` by replaying the literal
producer. With `c(i)=min(S_i)`, it selects all early rows meeting both conditions
above. It also implements improving single-source changes of the primary
assignment, and a sparse-reference heuristic that adds complete missing-
reference bundles for candidate rows. These are finite search heuristics;
neither assignment optimality nor moment optimality is claimed.

| Arity | Producer rows | Early rows, least-point rule | Early rows, locally improved rule | References |
| ---: | ---: | ---: | ---: | ---: |
| 5 | 58 | 16 | 16 | 0 |
| 6 | 186 | 43 | 60 | 0 |
| 7 | 436 | 113 | 151 | 0 |
| 8 | 824 | 262 | 356 | 0 |
| 23 | 40,077 | 12,109 | 12,109 | 0 |

The arity-23 least-point plan detours 1,769 of the 1,771 source rows through
their primary point-star frames. A source not read early follows its original
direct route. The same arity-23 producer has these sparse-plan counts:

| Reference rows | Early producer rows |
| ---: | ---: |
| 0 | 12,109 |
| 16 | 12,742 |
| 64 | 13,527 |
| 110 | 14,234 |

Every reference is used, stays at one chosen `U_c`, and has its full path
`0 -> source line -> U_c -> D1 -> global sink`. Its arbitrary original dirt
is explicitly included in the final sink unmixing. These sparse arity-23
rows have coefficient/incidence checks; the saved full dirty-basis replay
at arity 23 is for the zero-reference hierarchical variant below. Small
arity sparse plans also pass full arbitrary-dirt replay.

## Cleaning at the existing last frame, without copies

`hierarchy.py` improves some zero-reference schedules without adding roles.
Instead of taking every selected row from its last frame `A_s` to `U_c`, it
chooses a subset that can clean directly at `A_s`.

For each unchanged source `X_i`, collect all selected last frames that need
to read it. Require their source-support labels to form a nested chain, all
within its primary point star. Route `X_i` from its source line through this
chain, then through its primary `U_c` if needed, then to `D1`. Perform direct
cleanups in increasing dimension, followed by the remaining point-star
cleanups. Producer rows wait at their existing last frame until cleaned.
All producer and scatter uses finish before this cleanup sequence begins.

Support inclusion implies actual source-span inclusion. These frames are
positive definite for the inherited rational Gram form `G=I-J/9`: for triples
containing `c`, the Gram matrix is the unsigned pair-incidence Gram matrix on
the other points. Thus every nested nonzero frame difference is nondegenerate.
Equal-rank nested supports describe the same frame and require no motion.
This establishes the local common-frame schedule, but does not supply its
fixed-basis pivot profile.

At arity 23, the small-frame-first heuristic cleans 2,321 of the 12,109 early
rows at their original last frame. The remaining 9,788 clean at point-star
frames. It has 9,778 source-data local path edges from source lines through
the chosen frames and onward to `D1`. The large-frame-first heuristic cleans
1,935 directly with 7,045 source-data local path edges. These two objectives
can trade off; more direct cleanups need not mean a better moment.

The first-stage auxiliary exit of a directly cleaned row is

```
I - A_s tensor Q,
```

of rank `m-rank(A_s)`, where `Q` is the fixed other-factor line. Its **child
widths are not determined by this rank**. A complementary exit can split
into several pivot runs in the fixed basis. The source chains add new edges
whose complete profiles are also required. No single-child exit is assumed.

## Exact profile on the pinned PR84 word: no improvement

`indexed_hierarchy.py` applies the same schedule selection to the serialized
word of [Chafik Boukhalfa's PR84](https://github.com/CrocSwap/integer-mult-bounds/pull/84),
at pinned commit `88ca39571907343a49e97f328971ec7bcd26fbfd`. It uses the read-only
adapter and independent rational frame profiler in the adjacent sink-unmixing
package. The original compiler, word, and 23-by-25 rectangular basis are
credited external inputs; see the retained
[UPSTREAM-NOTICE](../indexed-cycle-target-10/UPSTREAM-NOTICE).
The checksum verifier checks the pinned input files and all 580 files in the
upstream source closure before using them.

For this producer, the least-point plan selects 6,461 of 26,834 producer rows
for early cleanup. The large-frame-first hierarchy cleans 1,472 directly at
their original frames. All 30,376 public scalar basis columns pass the
forward and same-chronology-transposed dirty-word replay, whose dense audit
realization has 2,338,001 CNOTs.

The profile audit computes 1,421 distinct exact rational merged-frame
profiles and 16,445 exact local frame differences. It replaces only the
first-stage paths and charges every selected source-data detour. The second
stage is either the original word or the separately audited borrowed word.
Both profiles retain the original deficit 1,846,900.

| Second stage | Public width W | Accepted finite bit saving | Numerical moment root |
| --- | ---: | ---: | ---: |
| Original | 132,466,108 | 0.000052143348 | 0.0000521433496954 |
| Borrowed | 128,392,808 | 0.000053451948 | 0.0000534519499464 |

The existing borrowed implementation in both stages has a larger finite bit
moment root, approximately **0.0000548890247**. Thus this fully accounted
hierarchy does not improve it. Even omitting the source-detour cost as a
diagnostic relaxation gives only approximately 0.0000536756403 for this
particular selected-row set. That diagnostic is not an implementable profile
and does not rule out other row selections or constructions.

The accepted entries in the table are exact rational finite-moment witnesses,
checked with rational logarithm and exponential enclosures. They are not
claims for the multiplication exponent `kappa`, nor a completed global
integration. The saved histogram has 46 child widths, with maximum 552;
replacing its merged exits by one child of the same total rank would give an
incorrect comparison.

## Checks and reproduction

The scripts use only standard Python and reject `-O`. Full scalar replay uses
integer bitmasks for all public basis columns, with no random field sampling.
It checks forward and same-chronology-transposed words, exact data independence
of the pre-sink auxiliary block, complete restoration of reference dirt, and
failure when a required cleanup is deleted. Sparse cases also reject omission
of the reference-dirt unmixing.

The arity-23 hierarchical word has 2,155,298 literal CNOTs and 43,619 public
scalar coordinates. Both complete basis replays pass. Dense early correction
is used for this audit; the gate count is not an optimized implementation.

```sh
python3 research/exploration/sink-reference-target-10/reference_plans.py 5 6 7 8 23
python3 research/exploration/sink-reference-target-10/reference_plans.py 5 6 7 8 --sparse --output research/exploration/sink-reference-target-10/sparse-small-results.json
python3 research/exploration/sink-reference-target-10/hierarchy.py 5 6 7 8 23
python3 research/exploration/sink-reference-target-10/hierarchy.py 23 --full-max 23 --order small-first --output research/exploration/sink-reference-target-10/h23-full-dirty-receipt.json
python3 research/exploration/sink-reference-target-10/indexed_hierarchy.py --source /path/to/pinned-pr84-checkout --order large-first
```

For sparse arity-23 plans, import `reference_plans`, build `M=model(23)`, set
`primary=[min(t) for t in M['triples']]`, then call
`sparse_greedy(M,primary,budget)` and `audit_plan(M,plan,False)` for budgets
16, 64, and 110. Their results are saved in `sparse-h23-results.json`.

The original `rtgm` schedules and the sparse-reference schedules remain at
the scalar/incidence stage. The selected PR84 hierarchy has a complete
finite rational profile and gives the negative comparison above. Further
work would require a better construction or selection, followed by the
inherited recursive, machine, and multiplication interfaces for any favorable
finite moment. The main theorem and existing exponent certificates are
unchanged.
