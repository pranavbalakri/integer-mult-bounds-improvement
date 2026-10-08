# Cross-characteristic motif search

These are finite geometry checks and scoped obstructions, **not a new integer
multiplication bound**. No new recursive child histogram or global theorem is
claimed here. The current published certificate is unchanged.

Run the finite check with Python's assertions enabled:

```sh
python3 research/exploration/new-motif-target-10/check_e8.py
```

The saved output is `e8-result.json`. The checker uses only Python's standard
library. It checks all 14,400 Gram entries, all 511 nonzero binary functionals,
and the 28-block incidence pattern used below. A nonzero minor modulo 101 proves
a rational rank lower bound; the ambient eight coordinates give the matching
upper bound.

## An E8 candidate and its scope

Take one representative of each of the 120 lines consisting of the roots

- `e_i + e_j` and `e_i - e_j`, for `i < j`;
- `(1, ±1, …, ±1)/2`, with an even number of negative signs.

Their rational Gram matrix `A` has rank 8, diagonal 2, and off-diagonal entries
0 or ±1. The binary matrix `B = A + J (mod 2)` has rank 9 and diagonal 1. At each
pair of distinct indices, exactly one of `A` and `B` is zero. Thus this supplies
a concrete cross-characteristic pair with 120 indices and ranks 8 and 9.

It does not yet improve the retained-centre construction. Every nonzero row in
the row space of `B` has support spanning all eight rational coordinates. The
checker exhausts these 511 rows: 255 have weight 56, 255 have weight 64, and one
has weight 120. Conversely, the support of every nonzero rational functional on
the roots spans all nine binary coordinates. To see the converse, suppose a
nonzero binary functional vanished on that support. Every root on which the
binary functional is 1 would then lie in the rational functional's kernel,
contradicting the checked eight-dimensional support span.

Consequently a **minimal-rank** factorization of either central matrix has
support-frame cost at least `8 * 9 = 72` per axis. In the current two-stage
schedule, with `v = 120` and `N = v²`, this already gives rank deficit at most

```text
N - 2 v * 72 = 14,400 - 17,280 = -2,880.
```

This statement does not exclude a nonminimal factorization with cheaper sparse
supports, a different endpoint schedule, or another motif. Counting only the
number of indices and the two ranks would miss this obstruction.

## The high-rank label variant also has a local obstruction

Here a rank-`k` label means a `k`-dimensional nondegenerate subspace of an ambient
nondegenerate symmetric bilinear space. A required zero Gram block means the
two subspaces are orthogonal. This is the natural block analogue of the scalar
label construction; more general nonsymmetric factorizations are outside the
claim.

The rational orthogonality graph contains an eight-clique, represented by
`e_1 ± e_2, e_3 ± e_4, e_5 ± e_6, e_7 ± e_8`. Therefore its ambient dimension is
at least `8k`.

The binary orthogonality graph also contains an eight-clique: the seven roots
`e_1 + e_j` for `j = 2,…,8`, together with `(1,…,1)/2`. Call the corresponding
mutually orthogonal nondegenerate `k`-spaces `E_1,…,E_8`, and write their
orthogonal complement as `Z`, of dimension `s = d - 8k`.

There are exactly 28 further roots indexed by the pairs `{i,j}` of these eight
clique vertices. Such a root has rational inner product zero with exactly
vertices `i,j`, and is therefore binary-orthogonal to the other six. Its
hypothetical subspace satisfies

```text
U_ij ⊂ E_i ⊕ E_j ⊕ Z.
```

The checker also verifies that two of these roots have rational inner product
1 when their pairs overlap, and 0 when the pairs are disjoint.

Let `W_ij = U_ij ∩ Z⊥`. Its dimension is at least `k-s`. Because `U_ij` is
nondegenerate and `W_ij` has codimension at most `s` in it, the radical of
`W_ij` has dimension at most `s`. A complement to that radical is a
nondegenerate space `V_ij` of dimension at least `k-2s`, contained in
`E_i ⊕ E_j`.

All 28 such `V_ij` are pairwise orthogonal. For overlapping index pairs this is
a required binary orthogonality; for disjoint pairs it follows from their
support in disjoint `E` blocks. They all lie in the nondegenerate `8k`-space
`E_1 ⊕ … ⊕ E_8`. Hence

```text
28(k - 2s) ≤ 8k,
s ≥ 5k/14,
d/k ≥ 117/14.
```

Together with the rational lower bound 8, this gives rank-product ratio at
least `468/7`, and index-to-rank-product ratio at most `70/39 < 2`. This rules
out the proposed route of compressing the E8 labels to effective ranks near
7 on both sides. It is not a ceiling for all possible recursive algorithms.

## Why a sparse-inner-product Leech subset is too small

One tempting replacement is a subset of the 24-dimensional Leech minimal
vectors, of squared norm 4, with off-diagonal inner products only 0 or ±1.
It cannot provide the desired thousand-plus lines: there are at most 425.
Here is a self-contained argument that applies to any such Euclidean set.

For a unit vector `u` in `R^r`, define its trace-free cubic tensor by

```text
T(u)_ijk = u_i u_j u_k
           - (u_i δ_jk + u_j δ_ik + u_k δ_ij)/(r+2).
```

Direct contraction gives

```text
<T(u), T(v)> = <u,v>³ - 3<u,v>/(r+2).
```

Thus the entrywise cubic kernel `K_3(t) = t³ - 3t/(r+2)` is positive
semidefinite on any collection of unit vectors. If its off-diagonal inner
products are 0 or ±α, and `G` is their ordinary Gram matrix, then

```text
K_3[G] = (1-α²) I - (3/(r+2)-α²) G.
```

When `α² < 3/(r+2)`, positive semidefiniteness bounds every eigenvalue of `G`
by `(1-α²)/(3/(r+2)-α²)`. Since `trace(G)=N` and `rank(G)≤r`, some eigenvalue
is at least `N/r`. Therefore

```text
N ≤ r(r+2)(1-α²)/(3-(r+2)α²).
```

With `r=24` and `α=1/4`, this is `N≤4680/11`, hence the integer bound 425.
For comparison, the E8 parameters `r=8, α=1/2` give 120 with equality.
The Leech restriction is therefore unsuitable before any circuit cost is
counted. This argument does not apply to an indefinite bilinear form or to a
family that permits additional odd inner products such as ±3.

The relative bound is also a specialization of Theorem 2.18 in the primary
paper [Ganzhinov and Szöllősi, *Biangular lines revisited*, 2021](https://d-nb.info/1234224763/34).
The calculations above are included so the particular obstruction does not
depend on accepting a numerical optimization or an external theorem.

## A finite search in indefinite geometry

`search_lorentz.py` investigates integer vectors `x` with a height `k` satisfying

```text
sum(x_i) = 3k,
sum(x_i²) - k² = 2.
```

On this hyperplane the Gram form is `x·y - k l`, equivalently
`xᵀ(I-J/9)y`. In dimensions greater than 9 this form is indefinite, so the
Euclidean relative bound above does not apply. A compatible family may have
any odd off-diagonal inner product, including ±3, and zero; it may not have a
nonzero even one. Its binary complementary Gram is again `B = A + J mod 2`.

The initial family is all positive triples and the difference lines
`e_i-e_j`, with `i<j`. It has `binomial(h+1,3)` indices. We searched two finite
catalogs:

- `positive-level3` adds all six-set vectors of height 2 and all permutations
  of `(2,1⁷,0,…)` of height 3. It deliberately does **not** include every root
  of height at most 3.
- `complete-level4`, for the dimensions used here, additionally includes
  `(-1,1¹⁰,0,…)`, `(3,1⁹,0,…)`, and `(2³,1⁶,0,…)`. The implementation also
  contains `(2²,-1,1⁹)` when dimension 12 permits it. These exhaust the
  coordinate patterns through height 4 for dimensions at most 12, up to
  simultaneous sign reversal. This follows by using
  `sum(x_i(x_i-1))=k²-3k+2` to enumerate the possible coefficients.

The search repeatedly removes a random fraction of the incumbent, sometimes
forces another root into it, and greedily fills it with compatible roots.
Every tenth trial starts the greedy fill from an empty set. The saved runs
found no improvement:

| Dimension | Catalog | Available roots | Trials | Best family | Rational / binary ranks |
|---|---|---:|---:|---:|---:|
| 10 | positive-level3 | 735 | 75,082 | 165 | 10 / 11 |
| 11 | positive-level3 | 2,002 | 21,751 | 220 | 11 / 11 |
| 12 | positive-level3 | 5,170 | 3,191 | 286 | 12 / 13 |
| 10 | complete-level4 | 1,585 | 79,371 | 165 | 10 / 11 |
| 11 | complete-level4 | 6,743 | 8,251 | 220 | 11 / 11 |

These are **heuristic search receipts, not upper bounds or optimality
certificates**. Larger families may exist within these catalogs or outside
them. The receipts store the seeds, trial counts, actual selected indices,
exact Gram checks, field ranks, and a hash of each candidate. The positive
rank deficit needed for the current two-stage schedule is not obtained by
these candidates. No side circuit or recursive multiplication histogram was
constructed from them.

To independently validate the saved candidates without repeating the search:

```sh
python3 research/exploration/new-motif-target-10/search_lorentz.py --verify-receipts
```

To also repeat all the recorded randomized search trials with their saved
seeds (this takes several minutes):

```sh
python3 research/exploration/new-motif-target-10/search_lorentz.py --verify-receipts --rerun-search
```

A shorter new search can be run with `--dimension 10 --catalog complete-level4
--trials 1000`. Running another heuristic search cannot establish optimality,
even if it returns the same family.
