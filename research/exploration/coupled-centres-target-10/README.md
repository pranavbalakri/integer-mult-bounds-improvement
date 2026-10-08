# Coupled tensor centres: exact scalar ingredients, unresolved geometry

This directory investigates a change to the centre-transfer architecture. It
contains exact scalar identities, a shared arbitrary-input auxiliary word, and
an exact *hypothetical* moment budget large enough for the target. It does not
contain a valid frame schedule attaining that budget and does not improve the
published exponent. The certified construction is unchanged.

## The two separate corrections cannot simply be merged

Let `B_h` be the binary point-by-triple incidence matrix,
`C_h=B_h^T B_h`, and `S_h=I+C_h`, the intersection-one side matrix. On the tensor
index set of size `N=v_a v_b`, define

```
C1 = C_a tensor I,       C2 = I tensor C_b,
S1 = I+C1,               S2 = I+C2.
```

The original stages apply `Y += X` and then `X += Y`. Keeping only their side
parts applies `Y += S1 X` and then `X += S2 Y`. Relative to the desired result,
the error on arbitrary original data `(x,y)` is

```
E(x,y) = ((C1+C2+C1 C2)x + C2 y, C1 x).
```

Because the two centre matrices commute, its image is exactly

```
{ (u+z,u) : u in im(C1), z in im(C2) }.
```

Indeed set `u=C1x` and `z=C2(x+y+u)`; after choosing `x`, the arbitrary `y`
realizes every `z` in `im(C2)`. Thus

```
rank(E) = a v_b + b v_a,
rank(C1 C2) = ab.
```

The shared mixed totals alone do not span the complete omitted correction.
This is a scalar rank statement, not a lower bound on the cost of every possible
geometric implementation. In particular, it does not rule out transporting the
same number of independent values along a cheaper joint schedule.

## A complete arbitrary-dirty mixed-centre word

For `h=2 mod 4`, `B_h B_h^T=0`: its diagonal is `binomial(h-1,2)` and its
off-diagonal is `h-2`, both even. Consequently, if both arities have this form,
`C1^2=C2^2=0`. Write `L(A)` for `Y += A X` and `U(B)` for `X += B Y`.
The four chronological operations

```
L(C1), U(C2), L(C1), U(C2)
```

have total map

```
diag(I+C1 C2, I+C1 C2).
```

Thus interchanging the two centre corrections changes the word only by a mixed
correction of rank `2ab` on the two data banks.

Factor `C1 C2 = H V` with

```
H = B_a^T tensor B_b^T,     V = B_a tensor B_b,     V H = 0.
```

A single shared bank `Z` of `ab` arbitrary auxiliary values implements the
mixed update on one data bank by the literal chronological block additions

```
X += H Z;
Z += V X;
X += H Z;
Z += V X.
```

The final state is `(X+HVX,Z)`, with `Z` restored for every possible initial
value. The same bank can then serve `Y`. No zero initialization or discarded
public dirt is assumed. `scalar_coupling.py` checks this on the complete input
basis at `6 x 6` and `6 x 10`, checks the correction rank also on other small
rectangles, and rejects a word missing its final cleanup gather.

This is a genuine scalar mechanism for sharing a centre bank. It does **not**
yet establish common frames for its gathers and scatters, nor remove the two
original centre transfers. A commutation identity alone does not make either
individual correction free.

## Why product centres would be large enough

There is a valid binary fitting matrix on pairs of triples:

```
C_a tensor C_b.
```

For intersection sizes `s,t`, its coefficient is `(s mod 2)(t mod 2)`, while
the normalized rational tensor Gram entry is `(s-1)(t-1)/4`. Every nonzero
off-diagonal central coefficient therefore has rational Gram zero. The full
side matrix is

```
I + C_a tensor C_b
  = S_a tensor I + I tensor S_b + S_a tensor S_b.
```

The point-pair source span is `U_c tensor U_d`, of dimension
`(a-1)(b-1)`. Each local point-star span is nondegenerate under `I-J/9` for
integer arities other than 9, so its tensor product is nondegenerate as well.
There are `ab` such spans. Their total copied rank would be

```
L_product = ab(a-1)(b-1),
L_product/N = 36/((a-2)(b-2)),
```

instead of the two separate-stage loss

```
L_old/N = 6/(a-2) + 6/(b-2).
```

At `23 x 25`, this replaces `2226400` copied ranks by `303600`, about an 86%
reduction. However, using this global fitting motif in the unchanged two-stage
construction would square its ambient dimension again, to `(ab)^2`. It is
incorrect to claim this loss while automatically retaining the old `m=ab`
paths and child profile.

To measure the potential value of a successful new routing, the diagnostic
`tensor_relaxation.py` nevertheless grants precisely that unproved combination:
it retains the fully batched data paths and `rho=2` flag-shaped auxiliary
profiles from the borrowing relaxation, then replaces their centre cost by
`L_product`. Exact moment enclosures give:

| Arity pair | Root of the hypothetical bit moment |
| --- | ---: |
| `10 x 14` | between `0.001301334732` and `0.001301334733` |
| `12 x 12` | between `0.001314198811` and `0.001314198812` |
| `23 x 25` | between `0.000779644462` and `0.000779644463` |

The first two hypothetical profiles pass the clean witness `a_bit=1/800`,
above the required `1/1023`. These are exact budgets for an explicitly
unrealized profile, **not** certified savings. The `10 x 14` pair also satisfies
the nilpotence conditions for the shared dirty scalar mechanism above.

## The direct shared-bank frame route fails

A data row indexed by `(T,U)` belongs to nine point-pair stars,
`{U_c tensor U_d : c in T, d in U}`. Directly gathering and scattering through
these frames requires visiting a `3 x 3` grid of subspaces. A switch between
stars differing only in their second coordinate has rank cost `2(a-1)`;
differing only in the first costs `2(b-1)`; differing in both costs
`2(a+b-3)`. These costs grant a route through the sum of the two frames, which
can only be cheaper than requiring the full ambient frame.

Any path visiting all nine has at least eight switches and at least two changes
in the more expensive coordinate. A snake attains the resulting bound

```
minimum switching cost
  = 12 min(a-1,b-1) + 4 max(a-1,b-1).
```

The initial source-line entrance and final full-frame exit already total
`ab-1`, exactly the direct data path's rank. All the switching cost is extra;
in the symmetric case it is `16(h-1)` per data row. This overwhelms the
available credit, which is less than one rank per data pair. The bound applies
to this literal visit-to-all-stars implementation, not to a different DAG,
shared retained-source construction, or a changed endpoint contract.

The promising remaining question is therefore specific: can the centre
commutation identity relocate the original two corrections so that existing
retained-source streams share useful **nested** frame paths? The scalar
identity, the product fitting matrix, and the hypothetical moments do not
answer that question by themselves.

## Additional interface audits

[The delayed-inverse audit](delayed-inverse-obstruction.md) gives an exact
zero-mixed-totals witness with a nonzero cleanup residue. It also explains why
a fresh retained-output copy generally includes public dirt and cannot be
treated as a clean point total.

[The endpoint-merging audit](borrowed-endpoint-merging.md) checks the separate
PR96 proposal against the complete borrowed word. All remaining auxiliary
columns participate in the early correction, and all participate in the late
full-frame inverse, so merging only the positive producer endpoints would miss
actual required gate incidences.

## Reproduce

From the repository root:

```sh
python3 research/exploration/coupled-centres-target-10/scalar_coupling.py
python3 research/exploration/coupled-centres-target-10/tensor_relaxation.py
```

Both scripts use standard Python, write exact JSON receipts, and reject
execution with assertions disabled. The moment checker is the repository's
independent rational exponential enclosure; no external producer is imported.
