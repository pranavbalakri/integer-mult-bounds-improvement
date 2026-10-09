# Exact first-stage profile: the full reference variants lose

The new scalar word is valid, but its `3v` and `2v` dirty-reference versions
do not improve the checked borrowed-source construction. The reference paths
cost more than the retirement merge saves. This is a useful negative result
about this specified schedule, not an impossibility result for sink unmixing.

The comparison uses the complete pinned PR84 producer at commit
`88ca39571907343a49e97f328971ec7bcd26fbfd`, with its credited reversed `23×25`
basis and inherited baseline profile. Its producer has `R=26834` first-axis
roles, of which `21521` are eligible to retire at a common point-star frame.
The full tensor dimensions are `m=575`, `N=4073300`; the first-axis invocation
is repeated `2300` times. All rows in the following table have the same rank
deficit `1846900`.

| First stage | Second stage | Total roles W | Approximate bit-moment root |
|---|---|---:|---:|
| Original | Original | 132466108 | 0.00005227402152 |
| Borrowed source | Borrowed source | 124319508 | 0.00005488902473 |
| Sink unmixing, 3v references | Original | 144686008 | 0.00004887671772 |
| Sink unmixing, 3v references | Borrowed source | 140612708 | 0.00005002469293 |
| Sink unmixing, 2v references | Original | 140612708 | 0.00004995562769 |
| Sink unmixing, 2v references | Borrowed source | 136539408 | 0.00005115546378 |

`indexed-profile-results.json` contains every complete histogram, a rational
accepted saving just below each displayed root, and exact rational exponential
moment bounds. These are finite-profile comparisons. They do not independently
establish the all-size multiplication interfaces, and none reaches the target
`1/1023` bit saving.

## Explicit rational frame matrices

Let `G=I−J/9`, `T=I+J`, and let columns of `L` span a local source space. Its
matrix in the actual fixed basis is

```
Π_L = T L (Lᵀ G L)⁻¹ Lᵀ G T⁻¹.
```

For an envelope with core `C`, cover `V`, and `O=V\C`, put
`c=|C|`, `n=|O|`, `s=3−c`, `d=s²+(c−1)n`, and `D=3(h+1)d`.
For `c=1,2`, write `o_i=1_(i∈O)`, `w_i=3+1_(i∈C)`, and
`z_j=3(h+1)1_(j∈C)−10`. The exact entries used in the checker are

```
(Π_(C,V))_ij = δ_ij o_i
  + [s o_i z_j + 3(h+1)s w_i o_j + n w_i z_j
     − 3(h+1)(c−1)o_i o_j] / D.
```

For a source triple `S`, the line projector is

```
(P_S)_ij = (3+1_(i∈S)) [3(h+1)1_(j∈S)−10] / [6(h+1)].
```

The point-star projector is `U_c=Π_({c},[h])`. These are matrices, not just
dimension labels. The first-stage lifted local matrix is `Π⊗Q`, where `Q`
is the second-factor source-line projector, with the inherited reversed
physical coordinate permutation. The early correction is at frame zero;
the local full frame is `E=I_h⊗Q`; the common auxiliary sink is `I_m`.

`frame_profiles.py` implements these formulas independently in rational
arithmetic. `check_frames.py` also reconstructs projectors from independent
triple bases using the displayed Gram definition. It checks twenty small full
tensor matrices by direct rational north-east elimination, and compares all
48 target point-star complements with a separate exact rank-one calculation.
`frame-check-results.json` includes explicit local and global matrices for
five `5×7` examples. Omitting the point-star's missing normal direction changes
the profile and is detected.

## Every changed path is charged

For a non-side producer role at last local projector `A`, choose an allowed
common point `c`. The old exit consists of `I_h−A`, followed by the selected
exterior children `{m−2h,h}`. Replace it by

```
U_c−A, then I_m−U_c⊗Q.
```

The cleanup occurs only after growing the role to `U_c`. The added dirty
reference path is exactly

```
0 → P_S⊗Q → U_c⊗Q → E → I_m.
```

For the `2v` version, each primary `X_S` also changes from the direct edge
`P_S→I_h` to `P_S→U_min(S)→I_h`; both differences are profiled separately.
Every source/reference injection occurs before these source paths grow.
The retained-centre copy costs, side-output growth, second-stage entrances,
data exterior edges, and endpoint correction are unchanged. The auxiliary
unmixing at `I_m` is pointwise and introduces no uncharged frame edge.

There is no claimed second-stage sink merge. The hybrid comparison uses the
already checked borrowed-source second stage, deleting precisely
`N×{525,25,1,1,23}` from the original histogram and `N` roles from `W`.

## The actual merged complement

For every first-axis common point, the rank-553 complement has profile

```
I_575−U_c⊗Q : {1,21,531}.
```

It is not a single child of width 553. The corresponding second-axis matrix
has profile `{1,23,527}`, but that profile is not charged as a second-stage
retirement in this construction.

Here is the independent complementary-corner calculation. For a projector
`X` of rank `r`, let `H=I−X⊗Q` after the physical permutation. If the corner
uses the top `i` rows and columns `j,...,m−1` (zero-based column index), then

```
rank H_[0:i,j:m]
  = rank X_[prefix(i),suffix(j)]                     if j ≥ i,
  = i−j−r + rank X_[prefix(j),:] + rank X_[:,suffix(i)] if j < i.
```

Every factor of the dense line projector is nonzero. The first and last `h`
physical coordinates visit the local coordinates once in their natural
order. Consequently these local sets are initial/final intervals or the
full set, and every needed rank follows from the local exact north-east pivot
set. Second differences of the corner ranks give the large matrix's pivots.
Consecutive diagonal pivot runs give its child widths. Direct full rational
matrix checks validate both this reduction and the physical permutation.

## A stronger scoped rejection

`check_uc_relaxation.py` lets every eligible producer role independently choose
any common-point anchor or remain at the old full-stage exit. It omits **all**
reference costs and **all** extra primary-source path costs. Even this relaxed
profile has moment strictly greater than one at the previous backed-off saving
`0.00005488`; the exact excess is at least `4.07369558×10⁻⁷`.

The omission is favorable: each reference path has total rank `m`, so its
positive children add a strictly positive amount to `moment numerator−W` at
positive saving. A nontrivial primary-source detour has at least two positive
children summing `h−1`; by concavity its contribution is at least that of the
old `{1,h−2}` edge. The checker compares all possible anchors using disjoint
rational moment intervals, including the option of not retiring a role.
It leaves 1079 roles late. The profile selected at that test exponent has
approximate root `0.00005400076183`; the rigorous claim is the exact rejection
at `0.00005488`, not a universal optimized-root formula.

This rules out gaining over the existing bound through reference sparsification
alone while every early cleanup still occurs at `U_c`. Cleanup at an earlier
original last frame, with genuinely nested source-data paths, is outside this
relaxation and is explored separately in `../sink-reference-target-10/`.

## Reproduce

```
python3 research/exploration/sink-unmixing-target-10/check_frames.py
python3 research/exploration/sink-unmixing-target-10/indexed_profile.py --source /tmp/kappa-audit/pr84-full
python3 research/exploration/sink-unmixing-target-10/check_uc_relaxation.py --source /tmp/kappa-audit/pr84-full
```

The source path must contain the pinned PR84 tree; its existing hash manifest
is checked. All new checks reject `python -O`. The inherited certified files
are unchanged.
