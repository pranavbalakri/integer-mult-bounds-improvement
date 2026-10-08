# One delayed cleanup cannot use only the mixed point totals

This is a scoped obstruction to moving the first producer's final source
cleanup past a second-bank shear. It is not an obstruction to all coupled
centre words. In particular, it assumes the second data bank is arbitrary;
it does not apply unchanged to a bank promised to start at zero.

Let `B_a` and `B_b` be the binary point-incidence matrices on all triples,
with `v_a=binomial(a,3)` and `v_b=binomial(b,3)` columns. Both have full row
rank. Set

\[
 C_2=I_{v_a}\otimes(B_b^{\mathsf T}B_b),\qquad
 T=B_a\otimes B_b.
\]

The mixed point centres record `T Y`, which has rank `ab`. If the second
shear changes the first bank from `X` to `X+C_2Y`, postponing an ordinary
source-slot erase until after that change leaves residue `C_2Y` in the
first producer's source slots. This statement assumes the usual invertible
producer has already undone its mixing, so those slots contain the old
copied `X` plus their original arbitrary contents. Erasing with the changed
bank adds exactly the displayed residue.

The residue map has rank `v_a b`, greater than `ab` for `a>=5`. More directly,
let `r` be the vector supported on the four first-coordinate triples

```text
012, 013, 024, 034.
```

Every point occurs evenly, so `B_a r=0`. Let `z` be the coordinate vector for
any second-coordinate triple. Then

\[
 T(r\otimes z)=0,
 \qquad C_2(r\otimes z)
   =r\otimes(B_b^{\mathsf T}B_bz)\ne0.
\]

The second assertion follows because the chosen triple's diagonal Gram
entry is one. Thus the same mixed-centre observation, zero, is compatible
with both zero residue and a nonzero residue. No deterministic correction
that sees only these centres can restore this delayed source bank.

The loss of information is exact. Writing
`W=(I_(v_a) tensor B_b)Y`, the residue is an injective embedding of `W`,
whereas the mixed centres reveal only `(B_a tensor I_b)W`. At least
`b(v_a-a)` additional independent binary linear observations are needed to
recover arbitrary `W`. This counts information, not recursive rank cost:
access to `Y` itself or a new factored producer may supply those observations
with a separately proved schedule.

## Why a zero-bank variant is not excluded

If `Y` is known to equal `C_1 X_0`, where
`C_1=(B_a^T B_a) tensor I_(v_b)`, then the late residue becomes

\[
 C_2Y=(B_a^{\mathsf T}B_a\otimes B_b^{\mathsf T}B_b)X_0
     =(B_a^{\mathsf T}\otimes B_b^{\mathsf T})T X_0.
\]

It is then recoverable from mixed point totals. This is a different data
promise, and any such architecture must account for the clean bank and
prove the desired multiplication transform on its actual input space.
The arbitrary-bank rank obstruction makes no claim about it.

## Dirty retained snapshots require cancellation

A retained total in an arbitrary-scratch producer generally has the form
`d_(c,j) + sum_i B_a[c,i] X_(i,j)`, where `d` is the retained image of the
original scratch. A fresh private copy preserves both terms. Gathering
these copies along the second coordinate produces the desired mixed signal
**plus** `sum_j B_b[d,j] d_(c,j)`.

Since the producer is invertible and its retained rows are independent,
those original retained dirt values can be arbitrary. Their gathered image
has full rank `ab`; it cannot be silently discarded from the scalar identity.
One may gather a matching pre-copy snapshot to cancel it, or prove a clean
signal-extraction word. Either route must be included in the scalar word,
its frame chronology, and its child-cost accounting before a private-copy
producer is a certificate.

The exact checker `check_delayed_inverse.py` verifies the explicit witness
and rank ingredients at `(a,b)=(6,6),(10,14),(12,12),(18,18),(23,25)`.
The argument itself does not require the nilpotent congruence condition
used by the separate coupled-shear identity. It applies for `a>=5,b>=4`.
