# Block labels and deferred complex endpoints

This records exact algebra for two possible changes to the recursive interface.
The transported product below is a positive identity, but **these checks do not
certify a better multiplication exponent**. In particular they do not remove
the endpoint child from any published histogram.

Run the independent standard-library check with assertions enabled:

```sh
python3 research/exploration/block-labels-target-10/check_interfaces.py
```

`results.json` contains the checked receipt. All endpoint arithmetic is over
the exact field Q(i), with rational real and imaginary parts. Two address
directions suffice for the noncommutation witness; the statements below hold
at arbitrary address dimension.

## An exact product in the uncorrected endpoint encoding

The retained two-stage complex word has outputs

    A = -F b,
    B = F P^-2 a + F P^-1 b,

where P^4 = I and P^2 is an address translation. Factoring out the common F
leaves the invertible two-bank encoding

    T(P)(a,b) = (u,v) = (-b, X a + Q b),
    Q = P^-1,  X = P^-2 = P^2.

Here X is a permutation, so X(f g) = (X f)(X g) for pointwise products.
If (u,v)=T(P)(a,b) and (s,t)=T(P)(c,d), then

    T(P)(ac,bd)
      = (-us, (v+Qu)(t+Qs) + Q(us)).                     (1)

Indeed, v+Qu=Xa and t+Qs=Xc, while us=bd. Thus the second component
on the right is X(ac)+Q(bd), as required. Formula (1) uses three applications
of Q, two ordinary pointwise products, and additions. A naive componentwise
product of the two encoded banks is incorrect; the checker includes a
nonzero counterexample.

If P is a tensor product of b rank-one phase factors, applying Q directly by
b butterflies costs O(b V) scalar work on V entries. Accordingly, a **single
globally deferred** encoding can be handled with linear work per active
direction at the pointwise-product boundary. This identity alone does not
establish a full FFT or multiplication algorithm with fewer recursive
children. Its bit cost also includes the coefficient precision and exact
phase guard; O(b V) is a scalar-operation count, not a runtime prediction.

## Why an order-three role type does not yet repair the recursion

Let

    R = [[0,-I],[I,I]],          D(P) = diag(P,I).

Then the exact factorization is

    T(P) = P^-1 D(P) R D(P)^-1,
    R^3 = -I,
    T(P)^3 = -P I,  T(P)^6 = P^2 I,  T(P)^12 = I.       (2)

The outside P^-1 in (2) acts on both banks. Reducing the signs in R modulo
two gives the familiar order-three matrix for multiplication by a primitive
element of F4. Over Q(i), R has order six. Crucially, neither the prefactor
P^-1 nor the conjugation D(P) is independent of the phase direction.
The order-three observation therefore does not supply a uniform scalar role
type for the complex recursion.

There is a minimal explicit obstruction to keeping the old parent mixers.
Take two role pairs with encodings T(P0),T(P1), where P0 and P1 are phases on
two independent address directions. Let E=diag(T(P0),T(P1)), and let G add
the second role pair to the first. Block multiplication gives

    E G - G E = [[0, T(P0)-T(P1)], [0,0]],              (3)

which is nonzero. The checker evaluates (3) on the entire relevant basis
until it finds and records a nonzero witness. As a positive control, replacing
P1 by P0 makes the same parent mixer commute with the encoding. Two independent
rank-one directions are valid for this test: the two-triple labels formed
from two orthogonal weight-three lines and a shared weight-three line are
orthogonal, have weight nine, and restrict to exactly this two-axis phase
algebra after an orthonormal change of basis.

The existing algorithm uses different triple-conditioned directions and
ordinary scalar role mixers. If a recursive child returns its encoded output,
those mixers no longer implement the old logical word. To preserve that word
one must implement an encoded mixer E G E^-1; its off-diagonal block contains
T(P0)T(P1)^-1, which is itself address-dependent. Formula (1) at the final
product does not implement these internal encoded mixers.

Also, removing correction calls only at the root of a fixed recursive
algorithm changes a constant factor, not the exponent from its child moment.
Deleting the endpoint term from that moment requires a recursive contract
that accepts and propagates the encoding at **every depth**. Such a typed or
matrix-valued recursive contract remains open here. Equations (2) and (3)
rule out only the proposed drop-in use of the unchanged scalar mixers;
they are not an impossibility theorem for redesigned algorithms.

Address-dependent twiddle multipliers add another obligation: even for a
single phase C_w=aI+bX_w, conjugation by a diagonal twiddle gives
aI+bD X_w D^-1. The latter is a weighted translation. With several independent
directions the resulting invariant blocks can grow, so a constant-size
encoding cannot be assumed without proving a bound on that growth.

## A useful single-direction simplification, with a precise boundary

For one active direction write P=aI+bX, with a=(1+i)/2, b=(1-i)/2 and
X=P^2. Put S=[[I,0],[-bI,I]]. In this case there is the stronger identity

    T(P) = S diag(I,X) S^-1 R.                          (4)

It follows that two such encodings differ by

    T(P0) T(P1)^-1 = S diag(I,X0 X1) S^-1.             (5)

Thus, in this restricted one-direction case, implementing an encoded parent
mixer only calls for fixed scalar role changes and an address translation.
This is an actual simplification, although an address translation still
needs an implementation and cannot simply be declared free in the tape
model. The exact checker verifies (4).

The active endpoint in the recursive algorithm is generally a product of f
phase factors, not a single factor. Formula (4) uses the special identity
P^-1=bI+aP^2. That identity fails already at f=2: on the joint eigenvector
with both phase eigenvalues i, P=-I, so P^-1=-I while bI+aP^2=I. The checker
includes a two-address-direction counterexample. Hence applying (5) directly
to a batched endpoint would be incorrect. Two factors give eigenvalues
1,i,-1; three or more can give all four fourth roots of unity.

Keeping a separate two-bank encoding for each active factor is a different
proposal and needs its own width and cost analysis. Executing the old finite
word separately for every factor incurs work linear in the number of active
directions, so it cannot be assumed to retain the recursive exponent saving.

## A block-label compression limit for the full triple graph

Replace each triple line by a nondegenerate rank-r subspace U_S in an ambient
nondegenerate bilinear space of dimension d. Impose the same required
orthogonality as the original triple labels: U_S is orthogonal to U_T when
the two triples have even intersection. This is the natural projective
orthogonal representation; nonsymmetric factorizations are outside this
claim.

The following clique is the p=2 construction in
[Bukh and Cox, Lemma 12](https://arxiv.org/html/1802.00476v2),
with its elementary block-rank consequence included here.
Partition 4 floor(h/4) points into blocks of four. In each block take all four
triples. Two distinct triples in one block intersect in two points, and
triples from different blocks are disjoint. These 4 floor(h/4) labels must
therefore be pairwise orthogonal. Their subspaces have a direct sum: a vector
in one U_S and the sum of the others would be in the radical of U_S and hence
zero. Consequently

    d >= 4 r floor(h/4),      d/r >= 4 floor(h/4).       (6)

For h divisible by four, the existing scalar representation has d/r=h and
already attains (6). For other h this argument allows a reduction of less
than four in effective dimension; it provides no construction achieving one.
The checker verifies these exact orthonormal triple cliques for h=4,...,64.

A plain tensor blowup U_S=<t_S> tensor F_2^r has ambient dimension rh and
effective dimension h, so it does not beat this limit. A promising block-label
construction would need a changed or restricted motif, rather than simply
compressing the full triple side graph. This scoped obstruction says nothing
about the existence of such a different motif.

Written with OpenAI Codex assistance. The starting endpoint identity and
phase convention are inherited from the primary
[two-stage phase-transfer note](https://github.com/Paureel/integer-mult-bounds/blob/c82d09eb4781b68435e3fa3eb6beea72fa900fab/notes/two-stage-phase-transfer.tex).
