# What another basis of the same central matrix can and cannot improve

This is a restricted obstruction, not a lower bound on all phase networks.
It concerns the current matrix

  A[S,T] = (|S intersection T|-1)/2

on all triples of [h], a rank-minimal scalar factorization of that matrix, and the inherited convention that a retained sum's source label contains the binary span of every triple with a nonzero coefficient. It keeps the copied read's low frame at zero. Alternative fitting matrices, additional redundant centres, nonzero intermediate read frames, and nonmonotone phase schedules are not ruled out.

## Source-support lemma

For h>=8, let f(S)=sum(c_i:i in S), where the c_i are rational or complex numbers, not all zero. Let U_f be the span over F_2 of the weight-three indicator vectors t_S with f(S)!=0. Then dim U_f >= h-1.

Proof. If U_f is proper, choose a nonzero binary annihilator with support Q. Every triple meeting Q oddly must have coefficient sum zero. Write k=|Q|.

If both Q and its complement have at least three points, triples with one point in Q and two outside force every outside coefficient to equal b and every inside coefficient to equal -2b. A triple wholly inside Q then forces b=0, a contradiction. If Q is all points, all triple sums vanish and all coefficients are zero. Thus only k=1,2,h-2,h-1 are possible.

For k=1 or 2, the same one-inside/two-outside equations give c_i=-2b inside and c_i=b outside, with b nonzero. For k=h-1, all coefficients inside Q vanish and the one outside coefficient is nonzero. For k=h-2, all coefficients inside Q vanish and the two outside coefficients are opposite nonzero values.

In each of these four cases the support triples span exactly the indicated hyperplane:

- k=1: all triples avoiding the distinguished point span its coordinate hyperplane.
- k=2: triples wholly outside Q span all outside coordinates; a triple containing both points of Q supplies their sum.
- k=h-1: triples containing the distinguished outside point span a space of dimension h-1. Their differences span the even-weight subspace on the other points, and one odd triple adds the remaining dimension.
- k=h-2: the support consists of triples containing exactly one of the two outside points. Their differences give the sum of those two points and the even-weight subspace on Q; one odd triple again gives total dimension h-1.

The facts used here follow directly by subtracting pairs of triples. In particular, all weight-three vectors on at least four points span the full coordinate space. This proves the lemma, including the case U_f is already full. Notice that an h-1 dimensional span can be degenerate for the binary dot product; requiring a valid nondegenerate phase label can only increase the required label dimension.

## Consequence for a rank-minimal central factorization

Let B be the h-by-C(h,3) incidence matrix. Then

  A = (1/2) B^T (I-J/9) B.

The incidence matrix has full row rank over Q. For h!=9, I-J/9 is invertible, so rank A=h and its row space is exactly the code f(S)=sum(c_i:i in S). In a factorization A=G F with exactly h nonzero centre rows, G has full column rank, and a left inverse proves rowspace F=rowspace A. Every nonzero centre therefore has source-span dimension at least h-1. The sum of these dimensions is at least h(h-1).

The existing retained basis has summed label dimension h^2-h+1. Changing only to a different rank-minimal basis of this same central matrix can therefore save at most one unit of this quantity. It cannot supply the order-of-magnitude change needed for the current target.

The exceptional h=9 matrix has rank h-1; the same source-support argument gives a lower bound (h-1)^2. That is still larger than v/2 (64>84/2), so this exceptional rank drop does not make the copied-centre deficit positive.

All h point stars P_i=sum(x_S:i in S) have source-span dimension h-1 and span the incidence code. They meet this dimension lower bound. Expressing the current central map through them uses coefficients (1_{i in S}-1/3)/2, hence thirds/sixths; this is outside the current Gaussian-dyadic scalar format. In addition, a suitable producer and its frame costs would need to be proved. This is a possible tiny refinement, not a certified improved multiplication witness.

## Directions left open

A genuinely new central construction could use:

1. A different fitting matrix with diagonal one and zero off-diagonal entries on pairs of triples meeting in exactly one point.
2. More than rank(A) centres whose source supports span smaller binary subspaces, with cancellations outside the incidence code.
3. Nonzero common read frames. For a source span U and target span V, a compatible intermediate frame inside U intersection V-perp could charge the rank of their cross Gram form instead of dim U; nesting, nonalternation and the complete phase schedule must all be checked.
4. A different endpoint or tensor motif, outside the complete architecture obstruction in architecture-obstruction.md.

None of these is proved here to attain a_c>2^-10.
