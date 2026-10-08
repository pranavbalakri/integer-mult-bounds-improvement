# The exceptional h=9 quotient in a rectangular motif

The singular nine-coordinate augmented motif has a dyadic rank-eight rational
factorization. Although its symmetric centre deficit is negative, pairing it
with a larger augmented motif gives a positive rectangular endpoint deficit.
This is a finite candidate for reducing dimensions, **not a certified new
multiplication exponent**: no complete side circuit or child moment is proved
here, and the target `kappa > 2^-10` remains unmet.

The package also proves that adding higher-norm vectors to all 120 fixed E8
root lines cannot preserve the same integral-parity fitting rule. A separate
small result excludes one particular overcomplete-source strategy while
allowing the added source to lie outside the original row space.

Run:

```sh
python3 check_exceptional.py > exact-results.json
```

The standard-library checker verifies every fitting pair, the quotient
identity, source-frame costs, and 7,560 orthogonal-pair rank conditions used
in the overcomplete-source argument. Rational rank lower bounds use exact
nonzero minors modulo 101, with matching explicit upper bounds.

## An explicit dyadic quotient

Use the 120 labels from the nine-coordinate augmented motif: triples have
`x=1_S,t=1`; differences have `x=e_i-e_j,t=0`. Thus `sum(x)=3t`. The form
`G=I_9-J_9/9` has radical spanned by the all-ones vector. For `i=1,...,8`, set

\[
 r_i(x)=2x_i-t+x_9.
\]

A direct expansion gives

\[
 \sum_{i=1}^8r_i(x)r_i(y)
 =4(x\cdot y-t_xt_y),\qquad
 A_{x,y}=\frac18\sum_{i=1}^8r_i(x)r_i(y).
\]

The `r_i` values are integral (`0, +/-1, +/-2`). This realizes the rational
quotient in eight Euclidean coordinates, with all central coefficients
Gaussian-dyadic. The resulting 120 lines are E8 root lines in a choice of
sign convention.

The binary labels remain nine-dimensional: triples use `1_S`, and
differences use the complement of their pair. Every nonzero `r_i` source
has binary support span nine. Every binary coordinate source has rational
support span eight after quotienting. Both full ambient forms are
nondegenerate, so each side has centre cost `8*9=72`.

For the more general minimal-source statement and the original E8 finite
checks, see [the earlier E8 package](../new-motif-target-10/README.md).
The new fact used here is the explicit dyadic quotient and its possible
rectangular use, not a new claim that the symmetric E8 deficit is positive.

## Rectangular endpoint bookkeeping

Pair the exceptional factor with the odd `h` augmented motif from
[the augmented-centre package](../asymmetric-fit-target-10/README.md).
It has `v_h=binomial(h+1,3)`, both ambient dimensions `h`, binary centre
cost `L_b=h(h-1)`, and dyadic complex centre cost `L_c=L_b` when `h-3` is a
power of two, or `L_b+1` using the flag basis otherwise.

The rectangular data count is `N=120 v_h`. Its bit and complex ambient
dimensions are `8h` and `9h`. The respective endpoint deficits are

\[
 \Delta_b=N-72v_h-120L_b,\qquad
 \Delta_c=N-72v_h-120L_c.
\]

| Other h | N | Bit ambient | Complex ambient | Bit deficit | Complex deficit |
|---|---:|---:|---:|---:|---:|
| 13 | 43,680 | 104 | 117 | -1,248 | -1,368 |
| 15 | 67,200 | 120 | 135 | 1,680 | 1,560 |
| 17 | 97,920 | 136 | 153 | 6,528 | 6,408 |
| 19 | 136,800 | 152 | 171 | 13,680 | 13,680 |

Positive entries leave room before side-circuit charges. They do not bound
those charges, establish a moment inequality, or establish a practical
crossover. The inherited precision certificate is unchanged; a newly built
word would still need its own scalar-prefix and assembly audit.

## Saturation against adding higher-norm vectors

Fix the standard E8 roots: the `D_8` lines `e_i +/- e_j` and the half-root
lines `(s_1,...,s_8)/2` with an even number of negative signs. There is no
nonzero real vector `x` whose inner product with **every** fixed root is
zero or an odd integer. In particular, no new higher-norm root type can be
added while keeping these roots and the same rule `B=1+Gram mod 2`.

To prove this, all products `x_i+x_j` and `x_i-x_j` must be integers. Hence
all coordinates are integral or all are half-integral. The all-plus
half-root also makes `sum(x)/2` integral, so the coordinate sum is even.

In the integral case, two coordinates of the same parity have even sum and
even difference. Both must therefore be zero. At most one coordinate could
be odd, and every even coordinate is zero because there are at least seven
of them. The even-sum condition excludes a single odd coordinate. Thus
`x=0`.

In the half-integral case, exactly one of each pair's sum and difference is
even. Its required vanishing forces all absolute coordinates to equal the
same positive half-integer `a`. Write `x=a s`, where each sign is `+1` or
`-1`. The even-sum condition forces an even number of negative signs, so
`s/2` is an existing half-root. But

\[
 x\cdot(s/2)=4a
\]

is nonzero and even, a contradiction.

This is a saturation statement about the **fixed integral-parity rule and
all the old roots**. Deleting roots, rescaling old labels, changing the binary
fitting matrix, using block labels, or a different architecture is outside
its scope. The proof is elementary and included in full; no priority claim
is made for this property of the E8 root configuration.

## One arbitrary extra source cannot improve a natural toggle scheme

Let `b_i` be row `i` of the binary fitting matrix on the 120 roots. It is
one on root `i` and on the 63 lines orthogonal to it. Its rational support
rank is eight. Suppose a source factorization adds one arbitrary binary
row `f`, not necessarily in the fitting matrix's row space, and replaces
some distinct natural sources `b_i` by `g_i=b_i+f`. Require each new source
`g_i` to be supported in the hyperplane orthogonal to root `i`.

This particular recoding cannot reduce the sum of rational source-span
ranks. If only one source is replaced, rank subadditivity applied to
`b_i=f+g_i` already proves this.

For two or more replacements, `g_i(i)=0` forces `f(i)=1`. If roots `i,j`
being replaced were nonorthogonal, the support condition for `g_i` at
coordinate `j` would force `f(j)=b_i(j)=0`, a contradiction. Thus all `k`
replaced roots are pairwise orthogonal. Their presence in `supp(f)` implies
`rank_Q supp(f) >= k`.

For any distinct replaced pair `i,j`, every root other than `j` that is
orthogonal to `i` but nonorthogonal to `j` is forced into `supp(g_i)`: the support condition
for `g_j` forces `f` to be zero there, while `b_i` is one. These roots span
all of `i^perp`, of dimension seven. The exact checker verifies that finite
rank assertion for all `120*63=7560` ordered orthogonal pairs. Therefore
all `g_i` have rank exactly seven, and the new cost is at least

\[
 k+7k=8k,
\]

matching the cost of the replaced sources.

This does not resolve arbitrary overcomplete factorizations. In particular,
it does not cover several added rows, different support hyperplanes,
general codeword generators, or block-valued sources. It avoids the invalid
assumption that all added source rows must belong to the original row space.
