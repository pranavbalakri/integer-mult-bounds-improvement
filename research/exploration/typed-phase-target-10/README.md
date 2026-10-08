# Recursive phase types: exact closure and its costs

This explores a four-channel representation of the fourth roots of unity and
the two-bank endpoint encoding from the retained complex construction. There
are constructive finite identities, but **no improved recursive moment or
multiplication bound is claimed**. The main quantified result rules out a
specific proposed implementation: making the batched endpoint encoding free
by pointwise changes of a fixed number of channels around a common phase
transform. Nonlocal encoded mixers and redesigned recursive contracts remain
outside that result.

Run the standard-library check with assertions enabled:

```sh
python3 research/exploration/typed-phase-target-10/check_typed.py
```

`results.json` records exact operator checks over Q(i), checked channel-rank
examples for one through five active directions, and conditional rank budgets.
No search or numerical rounding is used. Modular nonzero minors give rank
lower bounds over Q(i); the explicit factorizations below give matching upper
bounds.

## Two-bank relative mixers have an explicit word

Let P have order four, with P^2 an address translation, and write

    T(P) = [[0,-I],[P^-2,P^-1]].

For commuting phases P0 and P1, direct block multiplication gives

    T(P1)^-1 = [[P1,P1²],[-I,0]],
    T(P0)T(P1)^-1
      = [[I,0],[P0²P1-P0^-1, P0²P1²]].                (1)

Thus an unequal-type source pair (u,v) can be converted to the other type by

    u' = u,
    v' = P0²(P1 u) - P0^-1 u + P0²P1² v.

This is a constructive word: copy u as needed, perform one P1 phase and one
inverse P0 phase, perform the indicated translations, add, and erase the
copies. If each P has f independent active directions, it uses two phase
children of rank f, in addition to the translation and scalar work. These
are upper costs of this literal word, not a lower bound on every possible
encoded mixer. Translations have not been declared free in the tape model.

An ordinary scalar gate adding one pair into another can use (1) on a copied
source pair before adding it. Keeping the old scalar gate without this
conversion is incorrect, as checked in
[the preceding endpoint note](../block-labels-target-10/README.md).

## The four-channel regular representation really does close

Let W4 be the unitary four-point Fourier matrix on the channel index, with
entries i^(-jk)/2. Define

    E(P) = W4^-1 diag(I,P,P²,P³) W4.                   (2)

For commuting P0,P1, this representation has exact closure:

    E(P0) E(P1)^-1 = E(P0 P1^-1).

On an address eigenvector where P has eigenvalue i^s, E(P) cyclically shifts
the four channels by -s with this Fourier convention. Consequently this is a genuine constant-dimensional
representation of Z4, not merely a formal naming convention.

The direct implementation of (2) is also explicit. Apply the fixed W4
channel transform pointwise; apply P and P^-1 to modes 1 and 3 and P² to
mode 2; apply W4^-1 pointwise. Modes 0 and 2 require identity and translation,
respectively. When P0 and P1 have orthogonal rank-f active spaces, the relative
phase Q=P0 P1^-1 has rank 2f up to a translation. The relative four-channel
word therefore uses two rank-2f phase children, with total phase rank 4f.
This costs more than the literal two-bank conversion (1).

This mode decomposition also prevents a misleading width credit. E(P) alone
does not apply a full phase to four arbitrary channels: its nontrivial phase
modes are P and P^-1. If F is the full e-dimensional phase and P is its
rank-f factor, the modes of F E(P) have active phase dimensions

    e, e-f, e, e-f,

after allowing translations and inverse phases. Their phase-rank volume is
4e-2f, rather than 4e. This observation is an accounting test; it is not a
general lower bound on algorithms for these operators.

There is also a boundary to a Clifford-only type tracker. In the address
Walsh basis, E(C)^-1 contains controlled increment modulo four. With channel
bits r0,r1 and control address bit z, its binary action is

    z' = z,  r0' = r0 xor z,  r1' = r1 xor (z and r0).

This is a non-affine permutation, verified by its nonzero second difference
in the checker. A computational-basis permutation that is Clifford must be
affine: conjugating every Pauli translation must again yield a fixed Pauli
translation. Hence even this one-direction regular lift is not a Clifford
operation on the address and channel bits. Quadratic reversible types remain
possible, but their implementation cannot be justified solely by closure of
the Clifford group.

## Free pointwise channel changes require exponentially many channels

Here is the precise model being tested. A q-channel implementation consists
of a pointwise input map B(y), q identical copies of F, and a pointwise output
map A(x). Rectangular B(y) and A(x) are allowed, so this model grants extra
zero channels and permits their discard. Such an implementation has kernel

    F(x,y) A(x) B(y).                                  (3)

Since every entry of F is nonzero, divide the desired kernel entrywise by
F(x,y). Its matrix, with output-bank/address and input-bank/address indices,
must then have rank at most q. This argument even grants arbitrary exact
pointwise scalar coefficients.

Take P=C^tensor f on f active coordinates, with C=aI+bX,
a=(1+i)/2 and b=(1-i)/2, and first restrict to those coordinates. Spectator
coordinates can be fixed. Put n=2^f, s_x=(-1)^wt(x), c=(-i)^f and d=a^-f.
For the uncorrected two-bank target F T(P), the normalized kernel is

    K(x,y) = [[0,-1],[c s_x s_y, d 1[x=y]]].            (4)

The lower-right n-by-n submatrix is d times the identity, so its rank is n.
In fact the complete 2n-by-2n matrix also has rank exactly n. One explicit
factorization is

    A(x) = [ -d^-1 1^T ; e_x^T ],
    B(y) = [ c s_y s | d e_y ],
    K(x,y) = A(x) B(y).

The upper-left entry vanishes because sum_x s_x=0. Hence the exact minimum
number of common-F channels in model (3) is

    q_min = 2^f.                                      (5)

This includes an actual construction attaining the bound. At f=1, its two
input and output channel maps are square and invertible; it recovers the
exceptional free two-bank change of coordinates. At f=2, four channels
suffice. At f=3, at least eight are necessary, so a fixed four-channel
extension cannot solve the batched problem in this model.

The witness is compatible with the inherited phase-frame setting. For an
orthonormal active space and orthonormalizable complement, an adapted binary
orthogonal address basis makes P a product of active coordinate phases and
F a product of those phases and signed spectator phases. The latter have
nonzero kernels and cancel in the same normalized-kernel calculation.
Changing the address labels in this argument is only a relabeling of rows
and columns of the rank witness; it is not an assertion that the basis
change is a free algorithmic operation.

The literal four-channel target F E(P) has an even larger rank in model (3).
After its constant channel Fourier transform, its four normalized kernel
blocks have ranks

    1, 2^f, 1, 2^f.

Its complete normalized kernel consequently has rank 2^(f+1)+2. The checker
verifies both rank formulas through f=5.

In the existing batching regime f grows proportionally to e/m. The exact
construction attaining (5) then needs 2^(e/m+O(1)) common-F channels or that
many streamed common-F calls. Streaming can reduce its space requirement,
but does not remove this time factor. It is not a fixed-width recursive
primitive and does not preserve the desired polynomial-in-e overhead.

This is not a lower bound for arbitrary algorithms using arbitrary phase
calls: (3) deliberately restricts the intervening transform to identical F
channels and the surrounding maps to pointwise operations. Nonlocal Clifford
mixers or a new typed primitive could evade it, but require their own child
costs and all-depth contract.

## A one-node normalization test for free phase charts

There is another exact identity. For D(P)=diag(P,I) and
R=[[0,-I],[I,I]], all phases commuting with F,

    F T(P) = D(P) [(F P^-1) R] D(P)^-1.                (6)

If a proposed interface declares the direction-dependent D(P) input and
output charts free, its target in that chart is two complement phases of
dimension e-f. With the unchanged full-phase child contract, the uncorrected
data paths each have exactly that total rank: the X path goes f to e and the
Y path goes 0 to e-f. Every original auxiliary path still totals e.

For any 0<a<1, splitting a monotone path cannot make its child potential
smaller than the potential of its total rank:

    sum_j r_j^(1-a) >= (sum_j r_j)^(1-a).

Summing over these paths shows that their full-phase children already
dominate the properly normalized two-complement-phase and auxiliary target.
Copied centre transfers add positive work. Thus the combination of **free
D(P) charts and unchanged full-phase children** has no positive recursive
moment saving. Counting the pair as two full e-dimensional transforms would
instead award phase work absent from the target.

At the checked h=18 instance, normalized by the active block count f, the
copied-centre rank is 501,024. The original corrected rank deficit is 164,832.
Deleting the endpoint while retaining the full-phase denominator would
appear to increase it to 830,688; using the free-chart target instead gives
-501,024. Separately, even granting that optimistic full-phase denominator,
the literal two-bank relative word (1) could tolerate fewer than 415,344
unequal-type mixers before its zeroth-moment deficit becomes nonpositive.
No claim is made here that this is the actual number of mixers in a closed
encoded algorithm.

This normalization test is intentionally scoped. If the children also return
further encoded types, their actual input/output maps and a conserved type
potential must be specified. One cannot repeatedly discard their missing
directions by applying the free-chart convention independently at each depth.

## Precision remains linear; recursive work is the unresolved cost

For any Walsh-diagonal fourth-root phase Q on e bits, its physical matrix
entries lie in 2^-e Z[i]. Therefore E(Q) from (2) has entries in
2^(-e-2) Z[i]. It is unitary on four channels, so its row l1 norm is at most
sqrt(4*2^e)=2^(e/2+1). Relative types retain these same bounds, because they
are again E(Q), rather than an uncontrolled product of unrelated matrices.

Consequently an actual finite encoded word with the common-frame prefix
invariant can use the same semantic-reset argument as the existing guard:
completed children have their exact endpoint formats, and the one active
child supplies the recursive term. The extra channel Fourier coefficients
only contribute a fixed dyadic allowance. A guard of the form
g(e)<=g(rho e)+O(e)+J is still linear in e for fixed rho<1 and fixed word J.
No new numerical layer-guard constant is asserted without a complete word
and its scalar-prefix audit. Neither precision closure nor the algebraic
closure (2) supplies the missing favorable child moment.

Written with OpenAI Codex assistance. Starting phase conventions and the
two-stage endpoint are inherited from the
[primary endpoint note](https://github.com/Paureel/integer-mult-bounds/blob/c82d09eb4781b68435e3fa3eb6beea72fa900fab/notes/two-stage-phase-transfer.tex).
