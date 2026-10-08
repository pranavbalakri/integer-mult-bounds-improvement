# A symmetry-breaking rank-seven fit and an invariant-block exclusion

This folder records an exact finite construction found while investigating
the target `kappa > 2^-10`. The target remains unmet. The construction lowers
the rational fitting rank of the eight-point triple graph from eight to
seven, and seven is optimal even for fractional block fits. Its centre
costs prevent a gain in the current two-stage retained-centre architecture.
No multiplication exponent or complete recursive child histogram is claimed.

The accompanying theorem rules out an improvement by block fits that depend
only on intersection size. It allows arbitrary block size and noncommuting
block coefficients, but it does **not** cover general symmetry-breaking,
nonsymmetric fits. The rank-seven construction demonstrates why that scope
restriction matters.

All checks use exact integer or rational arithmetic and the Python standard
library. From this directory, run:

```sh
python3 check_fractional_motif.py > exact-result.json
```

The receipt checks every ordered pair of the 56 indices, verifies the
finite clique and centre ranks, checks the invariant-block formulas for
`7 <= h <= 100`, and checks their incidence-counting identities by direct
enumeration at `h = 7, 8, 9, 10`. The all-`h` theorem below is proved
algebraically; the finite range of the numerical receipt is not its proof.

## The eight-point construction

Let the points be the eight vectors of `F_2^3` and index the matrices by all
56 three-element subsets. Each triple lies in a unique affine plane of four
points. Color a triple by the direction of that plane, identified with its
unique nonzero normal `ell` in `F_2^3`. There are seven colors, with eight
triples of each color.

Define the rational and binary fitting matrices by

\[
 A_{S,T} = \mathbf 1\{\operatorname{color}(S)=\operatorname{color}(T)\},
 \qquad B_{S,T}=|S\cap T|\pmod 2.
\]

Both have diagonal one. Two affine planes with the same direction are equal
or disjoint. Triples in the same plane either coincide or intersect in two
points; triples in its distinct parallel plane are disjoint. Thus `A` is
zero whenever `|S intersect T| = 1`, and every off-diagonal product
`A[S,T] B[S,T]` vanishes.

The seven color-indicator columns give `rank_Q A = 7`. Binary point
incidence gives `rank_F2 B = 8`, verified independently by exact elimination.
The following seven triples form a clique for the overlap-one constraints:

```text
012, 034, 056, 135, 146, 236, 245.
```

Their pairwise intersections all have size one. In any fitting matrix with
`k`-by-`k` identity diagonal blocks and zeros for overlap one, this principal
submatrix is `I_(7k)`. Therefore its rank divided by `k` is at least seven,
over any field and without a symmetry assumption. Our rational fit attains
this bound.

### Why its centre cost does not improve the current recurrence

Each color class consists of all four triples in each of two disjoint
four-point affine planes. Its eight binary incidence rows form two
invertible `J_4 + I_4` blocks over `F_2`: indeed `(J_4 + I_4)^2 = I_4`.
Consequently every color class spans `F_2^8`.

Every nonzero binary linear incidence functional is therefore nonzero on at
least one triple of **every** color. The rational labels of its support span
all seven coordinate directions. Conversely, every nonzero rational
color-constant function contains a whole color class in its support, whose
binary incidence span has rank eight.

This proves a basis-independent statement about the **fixed matrices above**.
Any minimal scalar factorization of the rank-eight binary matrix has eight
nonzero factor columns in its incidence column space, each of support-frame
rank seven. Its total centre cost is at least `8 * 7 = 56`.
Any minimal scalar factorization of the rank-seven rational matrix has
seven nonzero color-constant factor columns, each of support-frame rank
eight, giving the same bound. The displayed factorizations attain it.
The corresponding row-space argument is identical.

With `v = 56`, the current two-stage schedule therefore has

\[
 N=v^2=3136,\qquad L=56,\qquad N-2vL=-3136.
\]

The negative deficit rules out this use of these fixed matrices before
paying any side-circuit cost. It is not a lower bound for every compiler,
every nonminimal decomposition, or every block-label implementation.

## A sharp theorem for intersection-invariant rational blocks

Let `h >= 7`, `n = binomial(h,3)`, and let `D` and `T` denote the adjacency
matrices on triples for intersection sizes zero and two, respectively.
Consider

\[
 M=I_n\otimes I_k+D\otimes U+T\otimes V,
 \qquad U,V\in\mathbb Q^{k\times k}.
\]

This is precisely the class of block fits with identity diagonal and zero
overlap-one blocks whose remaining blocks depend only on intersection
size. Neither symmetry nor commutativity of `U,V` is assumed. Then

\[
 \frac{\operatorname{rank}_{\mathbb Q}M}{k}\ge
 \begin{cases}8,&h=9,\\ h,&h\ne9.\end{cases}
\]

Both bounds are attained.

### Incidence decomposition

For `i=0,1,2,3`, let `X_i` be the inclusion matrix from triples to
`i`-subsets: its entry is one when the subset is contained in the triple.
Write `V_i` for its column space and `V_-1 = {0}`. The spaces are nested:
summing the columns indexed by `i`-sets containing a fixed `(i-1)`-set gives
`(4-i)` times that set's incidence column.

These columns are linearly independent over the rationals for each `i`.
For `i=1`, a relation says the sum of three point coefficients is zero on
every triple; comparing triples makes all point coefficients equal, hence
zero. For `i=2`, suppose edge weights satisfy
`w_ab + w_ac + w_bc = 0` for every triple. For fixed distinct `c,d`,
subtracting the equations for `abc` and `abd` gives

\[
 (w_{ac}-w_{ad})+(w_{bc}-w_{bd})=0
 \quad(a,b\notin\{c,d\}).
\]

There are at least three eligible points, so each difference is zero.
Every row of edge weights is constant; symmetry makes the constants equal,
and a triple equation makes their common value zero. The cases `i=0,3`
are immediate. Thus `dim V_i = binomial(h,i)`.

The orthogonal differences `E_i = V_i intersect V_(i-1)^perp` have
multiplicities

\[
 (d_0,d_1,d_2,d_3)
 =\left(1,h-1,\binom h2-h,\binom h3-\binom h2\right).
\]

Both `D` and `T` are symmetric and preserve this filtration. Direct counting
shows their actions on `E_i` are the scalars

\[
 a_i=(-1)^i\binom{h-3-i}{3-i},\qquad
 b_i=(3-i)(h-3-i)-i.
\]

For completeness, the column of `D X_i` indexed by `R` is
`binomial(h-3-i,3-i)` times the indicator that `R` is disjoint from the
source triple. Inclusion-exclusion gives its top incidence coefficient
`(-1)^i`, with all other terms in `V_(i-1)`.
The corresponding column of `T X_i` takes value `(3-i)(h-3)` when the source
contains `R`, value `4-i` when it contains exactly `i-1` points of `R`, and
zero otherwise. Subtracting `(4-i)` times the sum of the incidence columns
of the `(i-1)`-subsets of `R` leaves top coefficient
`(3-i)(h-3)-i(4-i)=b_i`. At `i=0` the sum is empty. Symmetry then gives the
claimed scalar action on each `E_i`.

Consequently

\[
 \operatorname{rank}M=\sum_{i=0}^3 d_i\operatorname{rank}P_i,
 \qquad P_i=I_k+a_iU+b_iV.
\]

This identity uses only the simultaneous scalar eigenspaces of `D,T`;
it requires no simultaneous diagonalization of `U,V`.

### Rank inequality

The coefficient matrix for `P_2,P_3` is

\[
 \begin{pmatrix}h-5&h-7\\-1&-3\end{pmatrix},
 \qquad\det=-2(h-4)\ne0.
\]

Eliminating `U,V` expresses each `P_j`, for `j=0,1`, as
`alpha_j I_k + beta_j P_2 + gamma_j P_3`, with

\[
 \alpha_0=\frac{(h-1)(h-2)(9-h)}{12},\qquad
 \alpha_1=\frac{(h-2)(h-3)}4.
\]

Put `r_i = rank P_i`. Whenever `alpha_j` is nonzero, rank subadditivity
gives `r_j >= k-r_2-r_3`. If `h != 9`, both coefficients are nonzero and

\[
 \operatorname{rank}M
 \ge hk+(d_2-h)r_2+(d_3-h)r_3\ge hk.
\]

Here `d_2,d_3 >= h` for all `h >= 7`. When `h=9`, keep the `j=1` bound
and omit the nonnegative `j=0` contribution to obtain

\[
 \operatorname{rank}M
 \ge8k+(d_2-8)r_2+(d_3-8)r_3\ge8k.
\]

Finally, the scalar fit `(|S intersect T|-1)/2` attains the bound. If `X`
is point incidence, its matrix is

\[
 \frac12 X\left(I_h-\frac19 J_h\right)X^{\mathsf T}.
\]

The full column rank of `X` makes its rank equal to that of the middle
matrix: `h` except at `h=9`, where it is `h-1`. Tensoring this fit with
`I_k` proves sharpness for every block size.

## Scope and context

The construction and calculations here were developed with OpenAI Codex
assistance. The affine-plane coloring, Fano clique, and triple-incidence
eigenspaces use standard finite geometry and inclusion counting; this note
does not claim priority for those ingredients.

For the broader fractional-Haemers framework and related cross-characteristic
triple constructions, see [Bukh and Cox, *On a fractional version of the
Haemers bound*, especially Lemma 12 and Proposition 17](https://arxiv.org/html/1802.00476v2).
Their disjoint-four-block clique supplies a known binary-side obstruction;
the theorem above concerns the separately specified rational invariant-block
family. We do not infer a lower bound for unrestricted rational fractional
fits by averaging a fit over permutations: averaging need not preserve rank.

Unrestricted symmetry-breaking block fits and alternative circuit
architectures remain open routes. No numerical near-fit is used as evidence
for an exact construction or a multiplication speedup in this package.
