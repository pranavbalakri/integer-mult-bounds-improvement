# Centre rectangles: checked local constructions and remaining gap

This is exploratory work toward a saving above `2^-10`. It proves no new
multiplication exponent and changes none of the published certificates.
The rectangle checks below concern scalar identities and binary frame feasibility;
no improved full network has been built from those rectangles. A separate
[audit of the in-place producer](inplace-producer.md) checks the borrowing
construction over complex values, including its local frames, conservative
histogram and linear guard. That variant remains separate from the candidate
multiplication assembly.

## Valid nonzero read frames

A centre containing a scalar sum of source triples is held at a nondegenerate
binary label U containing every source triple with a nonzero coefficient.
Let V span its target triples. A common read frame E is legal only if
E is contained in U intersect V-perp and both E and the residual U intersect
E-perp admit the inherited phase kernels. In particular a nonzero alternating
residual is not an inherited whole-residual child. The paid return rank is
`dim(U)-dim(E)`.

For even h, let P_i be the sum of triples containing point i, and E_i the sum
of triples avoiding i. Their exact support spans are

    U(P_i) = (1+e_i)-perp,       U(E_i) = e_i-perp,

both nondegenerate of dimension h-1. Here 1 denotes the all-ones vector.
There are genuine read-frame savings for individual rectangles:

- A source E_i read only by targets avoiding j, with i != j, can use
  the read line <e_j>. Its return rank is h-2, with a valid nonalternating
  residual.
- A source P_i read only by targets containing j, with i != j, can use
  the read line <1+e_j>. Its return rank is also h-2 and valid.

These are concrete legal building blocks. The unresolved issue is assembling
few enough more general rectangles into a fitting matrix. The exact
point-star/exclusion dictionary is now ruled out at the critical h=14 cost
budget; see `point-star-dictionary.md` for the proof and finite certificate.

A closely related tempting rectangle fails the frame test. A source P_i read
only by targets avoiding i has cross-Gram rank h-2 and candidate read line
<e_i>. But the residual consists of the even-weight vectors outside i: it
is nonzero, nondegenerate and alternating. The current phase-child interface
does not implement it as a single child of rank h-2. Counting cross-Gram rank
alone would incorrectly accept this rectangle.

## Two complete scalar fitting identities

Write r=|S intersect T| for source and target triples. The inherited central
matrix A[S,T]=(r-1)/2 has the point-star factorization

    A[S,T] = sum_i ((1_{i in S}-1/3)/2) 1_{i in T}.

It uses h source stars. Each source has label dimension h-1, but its target
coefficient is nonzero on every triple. Consequently its common read is zero
and the total centre rank is h(h-1), one unit below the inherited
h^2-h+1. This identity alone does not show that the changed producer has a
competitive role count. With the current copied two-stage deficit formula it
only reaches zero deficit at h=14; it gives no positive deficit at h<=14.

The same matrix can be written

    A = (1/3) sum_i P_i P_i^T - (1/6) sum_i E_i P_i^T.

The apparently cheaper cross terms are precisely the alternating-residual
case described above, and there are more terms. This is not an improvement.

An alternative fitting matrix is

    B[S,T] = binom(r,2)/3
           = (1/3) sum_{pairs p} 1_{p subset S} 1_{p subset T}.

Both A and B have diagonal one and vanish whenever r=1. Thus their off-diagonal
corrections occur only between orthogonal binary triple lines. For B, each
pair star consists of h-2 mutually orthonormal triples, with an admissible
complement and admissible proper side residuals. This gives a particularly
simple exact fitting construction, but its centre rank is
`binom(h,2)(h-2)=3 binom(h,3)`, far too large for the copied two-stage deficit.

## A rank constraint that a candidate must respect

Suppose a centre has nondegenerate source label U and contributes nonzero
values to diagonal entries indexed by a set of triples spanning D. Then
D is contained in U. Nondegeneracy identifies U with its dual, so restriction
of the pairing to D is onto D*. Hence

    dim(U intersect D-perp) = dim(U)-dim(D).

Every legal common read frame is contained in this kernel. Its paid rank is
therefore at least dim(D). For example, a family of Fano-type triples can have
a very low internal Gram rank, but its degenerate source span cannot itself be
used as U. Completing the source label restores this diagonal-rank constraint.
This is a necessary local test, not a lower bound on all decompositions.

## Fixed rational coefficients need not destroy the linear guard

The point-star identities contain thirds and sixths. They are outside the
literal Gaussian-dyadic scalar word currently certified, but this is not an
intrinsic obstacle to a future exact implementation.

For any fixed rational scalar circuit, choose a fixed integer D divisible by
every denominator of its scalar-prefix matrices and scalar operations. Along
one recursive invocation, a completed framed prefix has a denominator dividing
D times 2^e times the incoming denominator. The common-Walsh frame argument
still bounds its magnitude by a fixed scalar envelope times 2^(e/2).
A completed child has its exact semantic phase map, with dyadic denominator
at most 2^t; the child's internal odd denominator therefore disappears at its
return. No factors D are multiplied across completed siblings.

Only the active stack accumulates these fixed odd factors. With a strict
child-width contraction rho<1, its depth is O(log e), so one can use an exact
common scale D^H, where H bounds the remaining depth. This adds only O(log e)
bits. Division or multiplication by a fixed integer is a linear scan of a
fixed-width integer record; the scale ensures divisibility at the needed
operations. Returning to the parent removes known exact denominator factors.
The guard recurrence remains

    G(e) <= 2e + J + max_child G(t),    t <= rho e,

for a fixed J depending on the new finite circuit. Therefore G(e)=O(e).
This is a guard extension for a *valid completed phase word*, not a proof that
an arbitrary scalar fitting decomposition supplies that word. Its constant
must be recomputed if an actual new circuit is adopted.

## Reproduction and scope

Run `python3 check_rectangles.py` in this directory. The checked dimensions
are h=6,8,10,12,14. Every ordered source/target pair is checked with rational
arithmetic. The script also checks the valid nonzero reads, rejects the
alternating residual, and verifies the pair-star source and side labels.
The receipt is `results.json`. Also run
`python3 check_rectangles.py --dictionary` for the exact h=14 dictionary
certificate in `dictionary-results.json`.

The formulas, local binary calculations and rational-guard argument above are
new exploratory analysis, written with OpenAI Codex assistance. They use the
phase-frame and copied-centre interfaces already credited in this repository's
main note and in `research/independent/guard-improvement`. No new exponent,
practical runtime or complete global integration is claimed.
