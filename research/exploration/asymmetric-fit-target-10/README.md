# An augmented triple motif with smaller retained-centre cost

This package supplies an exact finite candidate for a new compiler: add all
coordinate-difference lines to the usual triple-incidence labels. At odd
arity `h`, it has `v = binomial(h+1,3)` indices, rational and binary ambient
dimensions both `h`, and optimal minimal-scalar centre cost `h(h-1)` on both
sides. The central rational map has a dyadic factorization attaining this
cost when `h-3` is a power of two, notably at `h=19`.

This is **not a new multiplication exponent**. A complete dirty-scratch word,
side-producer cost, child histogram, and assembly certificate remain to be
proved. The target `kappa > 2^-10` remains unmet. The existing published
precision guard and main certificate are unchanged.

Run the exact checks from this directory:

```sh
python3 check_augmented.py > exact-results.json
```

The checker uses the Python standard library. It checks every ordered pair
at `h=11,13,15,17,19`, the central identities on a rational coordinate basis,
source ranks, nondegeneracy, and every binary-functional orbit. Modular
minors modulo 101 give exact rational rank lower bounds; explicit ambient
spaces and annihilators give the matching upper bounds. No floating-point
search result is used.

| h | Indices v | Tensor ambient h² | Binary centre cost | Dyadic complex centre cost | Complex endpoint deficit |
|---|---:|---:|---:|---:|---:|
| 11 | 220 | 121 | 110 | 110 | 0 |
| 13 | 364 | 169 | 156 | 157 | 18,200 |
| 15 | 560 | 225 | 210 | 211 | 77,280 |
| 17 | 816 | 289 | 272 | 273 | 220,320 |
| 19 | 1,140 | 361 | 342 | 342 | 519,840 |

The deficit in this table is only `v²-2vL`, before any side-circuit costs.
For comparison, `h=17` has the same 816 data labels as the original
18-point triple motif, but tensor ambient dimension 289 rather than 324.
Its dyadic complex centre deficit is 220,320 rather than 164,832. These
finite changes may help a compiler; they do not establish a better recurrence.

## Fitting matrices

Fix odd `h>9`. There are two types of indices:

- A triple `S`, with rational coordinate vector `x_S=1_S` and `t_S=1`.
- A pair `i<j`, with rational coordinate vector `x_ij=e_i-e_j` and `t_ij=0`.

In both cases `sum(x)=3t`. Use the rational bilinear form

\[
 G=I_h-\frac19J_h,\qquad A_{a,b}=\frac12x_a^{\mathsf T}Gx_b
       =\frac12(x_a\cdot x_b-t_at_b).
\]

Every diagonal entry of `A` is one. The unnormalized off-diagonal Gram
entries are in `{-1,0,1}`: for two triples this is intersection size minus
one; for two distinct differences it is a signed shared-endpoint indicator;
and the mixed case is `1_S(i)-1_S(j)`.

The binary label is

\[
 b_S=1_S,\qquad b_{ij}=\mathbf1+e_i+e_j\quad\text{over }\mathbb F_2.
\]

Thus difference labels are complements of pairs. Every label has odd
weight, and

\[
 B_{a,b}=b_a\cdot b_b
       =1+(x_a\cdot x_b-t_at_b)\pmod2.
\]

Consequently `B` has diagonal one and `A[a,b] B[a,b]=0` for every distinct
pair of indices. Both fitting matrices have rank `h`: the triples alone
span the coordinate space over each field, and `G` is nondegenerate when
`h != 9`.

The augmented root family itself was already considered in
[the earlier Lorentz-root exploration](../new-motif-target-10/README.md).
The additions here are its explicit binary coordinate form, centre bases,
source-frame optimality, and dyadic central factorization. The labels and
arguments are provided in full; no priority claim is made for the root-family
construction.

## Low-cost source frames

Let `f_i(a)=b_a[i]`. This binary central source is supported on triples
containing `i` and differences avoiding `i`. Put

\[
 q_i=\mathbf1-3e_i.
\]

Every rational label in that support lies in `ker(q_i^T)`. They span this
hyperplane: differences outside `i` span the sum-zero subspace there, and
triples containing `i` add the remaining direction. Its dimension is `h-1`.
The restricted rational form is nondegenerate, since

\[
 q_i^{\mathsf T}G^{-1}q_i=\frac{36}{9-h}\ne0.
\]

The `h` sources `f_i` are the coordinate columns of the full-rank binary
label matrix, so they are independent and factor `B` directly. Their total
centre cost is `h(h-1)`.

For the rational map use

\[
 z_i(a)=t_a-x_a[i].
\]

This has value one on triples avoiding `i`, value zero on triples containing
`i`, and a signed unit value on differences incident to `i`. Its support is
exactly where `b_a[i]=0`. Those binary labels span the ordinary coordinate
hyperplane `e_i^perp`, which is nondegenerate and has dimension `h-1`.
The `h` functions `z_i` are independent: their coefficient matrix on the
rational coordinates is `J/3-I`, invertible when `h != 3`.

### Subset-indexed frames for producer graphs

These endpoint frames extend to a useful nested family. For any subset `I`
of output coordinates, let

\[
 U_I=\bigcap_{i\in I}\ker(q_i^{\mathsf T}).
\]

The normal vectors are independent, and their Gram matrix in `G^{-1}` is

\[
 9I_{|I|}+\frac{9(h-5)}{9-h}J_{|I|}.
\]

The nonconstant eigenvalues are 9. The remaining one could vanish only at
`|I|=(h-9)/(h-5)`, strictly between zero and one. Hence every `U_I` is
nondegenerate, of dimension `h-|I|`, for `h>9`. The empty subset gives the
full space. Inclusion of subsets reverses inclusion of frames. On the
binary side the analogous intersections of coordinate hyperplanes are
also nondegenerate.

This supplies valid frames for a producer whose intermediate source
supports satisfy the indicated coordinate restrictions. It does not by
itself prove a producer's role schedule or its charged child widths.

## Dyadic central identities

Using all low-dimensional sources gives the exact identity

\[
 A_{a,b}=\sum_{i=1}^h
 \left(-\frac{x_a[i]}2+\frac{t_a}{h-3}\right)z_i(b).
\]

It follows by expanding `z_i=t-x_i` and using `sum(x)=3t`. When `h-3` is a
power of two, every coefficient is dyadic. At `h=19` the scatter coefficient
is `-x_a[i]/2+t_a/16`, and its possible values are
`0, 1/2, -1/2, 1/16, -7/16`. Both sides then attain centre cost 342.

For other odd arities a flag basis removes the odd denominator at a cost
of one extra frame dimension. Use `z_1,...,z_(h-1),t` and the identity

\[
 A_{a,b}=\frac12\sum_{i=1}^{h-1}
 (x_a[h]-x_a[i])z_i(b)
 +\left(t_a+\frac{3-h}{2}x_a[h]\right)t_b.
\]

All coefficients are integers or half-integers. The source `t` is supported
on every triple and has binary frame dimension `h`, so the total is
`(h-1)^2+h=h(h-1)+1`.

These formulas preserve the Gaussian-dyadic scalar format. Once a full word
exists its scalar-prefix bound and semantic guard still require checking;
this package does not extend the previously certified `128(d+1)` guard to
an unbuilt circuit.

## Optimality for these fixed minimal scalar factorizations

The following proof concerns these fitting matrices and the inherited
support-span definition of a centre. It does not exclude nonminimal
factorizations, other matrices, or a different architecture.

Consider a nonzero binary linear functional `a·b` of weight `k=|a|`, and
split the coordinate set into `K=supp(a)` and its complement. We claim its
support has rational span dimension `h-1` if `k=1`, and `h` otherwise.

If `k` is positive and even, the selected differences are exactly the
crossing pairs between `K` and its complement. Those differences span the
whole sum-zero hyperplane. A selected triple exists and has coordinate sum
three, giving full rank `h`.

If `k` is odd and smaller than `h`, the selected differences are those
within either part. Their span has codimension two. The selected triples
have either one or three points in `K`. When `k>=3`, both kinds exist,
since `h-k` is even and at least two. Their quotient coordinates are
`(1,2)` and `(3,0)`, which are independent over the rationals. The resulting
span is full. When `k=1` only the first kind exists, giving the stated
hyperplane of dimension `h-1`. Finally, `k=h` selects all triples, again
of full rank. This proves the claim for odd `h>=5` as a statement about the
raw rational coordinates; here `h>9` makes the ambient form nondegenerate.

Every nonzero row in a minimal binary factorization therefore has rational
support rank at least `h-1`. A rank-`h` factorization needs `h` such rows,
so its total cost is at least `h(h-1)`.

Conversely, suppose a nonzero rational incidence functional had binary
support rank at most `h-2`. Its binary annihilator would contain two
independent nonzero functionals, and hence some functional of weight at
least two. The rational functional must vanish on the support of that
binary functional, but the preceding result says that support spans the
entire rational ambient space. This is a contradiction. Thus every
nonzero rational source has binary support rank at least `h-1`, proving
the same lower bound on the opposite side. The `z_i` basis attains it.

## Structure available to a compiler

The correction matrix `P=I-A` is supported on binary-orthogonal source and
target labels. With pair coordinates oriented as `e_i-e_j`, its pieces are:

- For a triple target, the old triple-to-triple correction is unchanged.
  The extra difference contribution is negative one-half of the oriented
  cut sum between the triple and its complement.
- For a difference target `(i,j)`, the triple part is negative one-half of
  the triples containing `i` but not `j`, plus positive one-half of those
  containing `j` but not `i`.
- Between differences, only other pairs sharing one endpoint contribute,
  with coefficient negative one-half of their signed Gram entry. These are
  signed leave-one-out star sums.

The dense source `z_i` needs a careful producer: summing triples avoiding
`i` independently would repeat each triple `h-3` times. A monotone interval
producer is a possible alternative. The points omitted by a sorted triple
form four intervals; partition them into segment-tree nodes, aggregate
input triples by node, and propagate ancestor sums to leaf outputs. Along
any root-to-leaf path an input appears at most once. Pair corrections add
signed incident-edge terms. The actual graph, role reuse, frame changes,
scalar prefix bound, and recurrence moment are still required before this
idea can support a multiplication claim.
