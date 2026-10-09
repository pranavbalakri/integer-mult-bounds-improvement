# One sparse extra source direction cannot lower the E8 centre cost

This package investigates overcomplete source factorizations, where source
rows need not belong to the central matrix's row space. It proves a finite,
scoped exclusion and supplies a sharp ten-source example tying the old cost.
It does not prove unrestricted overcomplete optimality or a new multiplication
bound. The target `kappa > 2^-10` remains unmet.

Let `C` be the nine-dimensional binary row space of the E8 fitting matrix,
indexed by its 120 root lines. For a binary source row `u`, charge

\[
 w(u)=\dim_{\mathbb Q}\operatorname{span}\{r_i:u_i=1\}.
\]

**Result.** Fix any binary row `f` with at most 13 nonzero entries. Every
family of arbitrary rows in `C+<f>` whose span contains `C` has total cost
at least 72. No condition puts those rows inside `C`, and no restriction
is imposed on their rational support subspaces.

The bound is sharp: replace one natural Gram row `b_i` by the singleton
`e_i` and `b_i+e_i`. With eight other code-basis rows this gives ten
independent sources of costs `1,7,8,8,8,8,8,8,8,8`, totaling 72.
The checker includes the complete dictionary as bit masks.

The next search should therefore allow at least 14 entries in the added
direction, at least two independent extra directions, or another motif.
The result does not exclude any of those possibilities.

Run from this directory:

```sh
python3 check_sparse_extension.py > exact-results.json
```

All arithmetic is integral or rational, with modular minors providing exact
rational-rank lower bounds. The standard-library checker verifies all 511
nonzero codeword support ranks, the reflection and hyperplane reductions,
the punctured code tables, 43,920 moment-tensor entries, and the sharp
source dictionary. No numerical optimization result is used.

## Fixed geometry and code

Use the standard E8 roots of norm squared two: the 56 lines `e_i +/- e_j`
and the 64 lines `(s_1,...,s_8)/2` with an even number of negative signs.
Let

\[
 B_{ij}=1+r_i\cdot r_j\pmod2,
 \qquad C=\operatorname{rowspan}_{\mathbb F_2}B.
\]

The checker proves `dim C=9`, that every nonzero codeword has rational
support rank eight, and that the minimum nonzero Hamming weight is 56.
Thus a nonzero `f` of weight at most 13 lies outside `C`.

For each root `i`, write `H_i=r_i^perp`. The natural row `b_i` is the
indicator of root `i` together with all 63 root lines in `H_i`.

## The finite puncturing gap applies to every proper support space

Every proper span of root lines is contained in a rank-seven span of root
lines: extend an independent set of roots to seven using the full-rank root
configuration. It therefore suffices to consider such hyperplanes.

The checker supplies an explicit simple-root basis and verifies two facts:
its eight reflections permute the 120 root lines, and every root's
coordinates in that basis have one sign and are integral. These facts give
the following reduction without assuming an external root-subsystem list.
Each displayed simple root is explicitly checked to be one of the root
lines. The induced line permutations preserve `C`: choosing a representative
after reflection can change a root's sign, but this does not change its
integer inner products modulo two, so the fitting matrix is permuted in both
indices.

The reflection group is finite because it permutes a spanning set of root
lines. Move a hyperplane normal into a position maximizing its inner product
with a vector having positive products with every simple root. Reflecting
in a simple root with negative product would increase that objective, so
all eight products are nonnegative. A root is then perpendicular to the
normal exactly when its simple-root support uses only zero-product simple
roots. A rank-seven root span has exactly seven such simple roots.
Consequently its normal is proportional to one of the eight dual basis
vectors, up to a root permutation preserving `C`.

For each of these eight representatives, the checker restricts all 512
codewords to the roots **outside** the hyperplane:

| Deleted simple-root index | Root lines inside | Minimum nonzero weight outside |
|---|---:|---:|
| 1 | 42 | 14 |
| 2 | 28 | 28 |
| 3 | 22 | 35 |
| 4 | 14 | 43 |
| 5 | 16 | 40 |
| 6 | 23 | 33 |
| 7 | 37 | 27 |
| 8 | 63 | 1 |

The last case is a root-perpendicular hyperplane `H_i`. Its unique
exceptional codeword below weight 24 is `b_i`, whose outside restriction
is the singleton `i`. Every other nonzero codeword has outside weight at
least 24. All outside restrictions are injective on `C`, consistent with
the full support-rank property above.

Now suppose `u=f+c`, with `c` in `C`, has rational support rank below eight.
Place its support in one of the hyperplanes just classified. Outside that
hyperplane, `c` must equal `f`, which has weight at most 13. The table forces
`c=0`, unless the hyperplane is some `H_i` and `c=b_i`. Thus every low-rank
row in `C+<f>` is either

\[
 f\quad\text{or}\quad f+b_i,
\]

where the latter possibility requires `f_i=1` and every other root in
`supp(f)` to be orthogonal to root `i`.

Let `I` be the set of those eligible roots. They are pairwise orthogonal,
so they are linearly independent and

\[
 |I|\le w(f)\le8.
\]

## The exceptional toggles still cost seven

A row `f+b_i` from the preceding list has support obtained from the 63 roots
in `H_i` by deleting at most 12 of them. These remaining roots still span
`H_i`. Here is an exact, reproducible bound stronger than needed.

For doubled integer root coordinates `R_j=2r_j`, set
`Q_i=8I-R_iR_i^T`. On the 63 roots in `H_i`, the checker verifies the full
symmetric tensor identities

\[
 \sum_j R_jR_j^{\mathsf T}=9Q_i,
 \qquad
 \sum_j(R_j\cdot z)^4=3(z^{\mathsf T}Q_i z)^2.
\]

The second identity is checked coefficient by coefficient. If `z` is a
nonzero vector in `H_i`, its second and fourth moments are respectively
`72 ||z||^2` and `192 ||z||^4`. Cauchy--Schwarz therefore implies that at
least

\[
 \frac{72^2}{192}=27
\]

of those roots have nonzero inner product with `z`. Removing fewer than
27 roots cannot put the survivors in a proper subspace of `H_i`.
In particular every eligible `f+b_i` has cost exactly seven.

## Cost proof for an arbitrary source family

A family entirely inside `C` needs at least nine independent nonzero rows,
each costing eight. Assume instead that it contains a row outside `C`.
Because it spans `C` and is contained in the ten-dimensional space
`D=C+<f>`, its span is exactly `D`. Discard dependent rows to obtain a
basis of `D`; this can only decrease total cost.

Every basis row other than `f` and the eligible `f+b_i` has cost eight.
If the basis contains `f`, let `k` be the number of its other low-rank rows.
The total cost is at least

\[
 w(f)+7k+8(9-k)=72+w(f)-k\ge72,
\]

since `k<=|I|<=w(f)`. If the basis does not contain `f`, it has at most
`|I|<=8` low-rank rows, so its cost is at least

\[
 7k+8(10-k)=80-k\ge72.
\]

This proves the stated result. It also explains why allowing the added row
to be arbitrary does not invalidate this particular exclusion.

The E8 fitting matrix and minimal-source context are documented in
[the earlier E8 package](../new-motif-target-10/README.md) and
[the exceptional quotient package](../exceptional-rectangle-target-10/README.md).
The root ingredients are standard; the finite reduction, tables, and full
argument are included here so the scope can be checked independently.
Developed with OpenAI Codex assistance; no priority claim is made for the
root-system or moment identities.
