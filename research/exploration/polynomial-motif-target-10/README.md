# The full five-subset polynomial motif

This is a rejected candidate for the retained two-stage transfer, not a lower
bound on arbitrary multiplication algorithms. It tests a different family
from the earlier search for five-subsets avoiding intersection three.
The use of complementary polynomial fitting matrices follows the standard
characteristic-separation approach; see [Bukh and Cox, Section 4](https://arxiv.org/html/1802.00476v2).

For all five-subsets S,T of [h], put s=|S intersect T| and take

    A(S,T)=(s-1)(s-3)/8 over Q,
    B(S,T)=s mod 2 over F2.

Both diagonals are one. For distinct sets, A vanishes at intersections 1 and
3, while B vanishes at 0, 2 and 4. Thus this pair has exactly the complementary
zero patterns needed for the two scalar motifs; no subfamily restriction is
needed. The decomposition A(s)=binom(s,2)/4-3s/8+3/8 uses pair, point and
constant incidence terms with dyadic coefficients. Thus an uncompressed
central implementation does not require odd-denominator scalar arithmetic
in the complex network; no assertion about a smaller complex circuit is made.

## Rational labels and central cost

Let phi(T) be the indicator of the ten pairs contained in T. Write r=binom(h,2),
and let U be the pair-versus-point incidence matrix. On Q^r define

    H = I/4 - 3 U U^T/128 + 3 J/800.

Then phi(S)^T H phi(T)=A(S,T), and each label has norm one. The eigenvalues of
H on the constant, point-difference, and remaining pair subspaces are

    (3h^2-78h+475)/1600, (38-3h)/128, 1/4,

with multiplicities 1, h-1, and binom(h,2)-h. They are all nonzero for integer
h>=7: the first numerator has discriminant 384, not a square, and the second
has no integer zero. In particular the rational label dimension is r.

Keep the binary point-star central factorization of B. The span of labels of
five-subsets containing a fixed point i has dimension binom(h-1,2). To see
this, project to pairs avoiding i. These are the incidence vectors of all
four-subsets of h-1 points. They span the full pair-coordinate space for
h-1>=6. Indeed, if coefficients c_jk sum to zero on every four-subset,
comparing two four-subsets with three common elements gives

    sum_{j in K}(c_jp-c_jq)=0

for every three-subset K avoiding p,q. Differences of such equations make
all c_jp-c_jq equal; the displayed equation makes this common value zero.
Hence every c_jk is equal, and the four-subset equation makes it zero.
The coordinates on pairs containing i are recovered from the others by
dividing each incident-coordinate sum by three. This proves the span claim.

A retained-centre copy to the zero frame must therefore pay at least
binom(h-1,2) rank units for each of the h centres. Possible degeneracy of a
star span could require a larger nondegenerate completion; it cannot lower
this bound. In fact the star restriction is singular at h=11 and h=21:
its eigenvalues on outside-pair coordinates are 1/4, (21-h)/72 and
(11-h)/36. We grant the smaller span dimension even at these exceptional
values; we do not claim that all the literal centre frames are valid.

## Why the rank gain still misses the target

Use the unchanged symmetric two-stage data endpoints and their rank-one
endpoint correction, at m=r^2, with v=binom(h,5) and N=v^2. Retaining these
centre copies bounds the deficit per data index by

    d <= 1 - 2h binom(h-1,2)/v
       = 1 - 120/((h-3)(h-4)).

For 7<=h<=14 this is nonpositive. The smaller cases h=5,6 have respectively
v=1,6 and point-star span dimensions 1,5, so their deficits are both -9.
For h>=15, r>=105. Grant the construction
zero auxiliary roles and optimally batch the unchanged data paths; this
relaxation only makes its moment better. Their entropy alone is at least

    E >= (4r-2) log r > 4*(4*105-2) = 1672,

because log(105)>4. This is the data-path relaxation also used in
[the borrowed-source ceiling](../borrowed-ceiling-target-10/README.md): per
data index, it has two blocks of width (r-1)^2, four of width r-1, and one
singleton. The last five blocks contribute

    (4r-2) log r + 4(r-1) log(r/(r-1)).

Dropping the second term and the positive corner-block entropy is safe.
If the exact moment at a positive saving a is below one, exp(x)>=1+x
requires a E<d<1. Consequently

    a_bit < 1/1672 < 1/1023.

The final inequality excludes kappa>1/1024 under the retained assembly.
This family improves index density by a higher-degree fitting polynomial,
but its rational address space grows too much. Unequal factor dimensions,
other factorization costs, asymmetric block fits, and other transfer
contracts are not claimed to be covered here.

## Check

`python3 research/exploration/polynomial-motif-target-10/check_polynomial.py`

The checker uses exact rational arithmetic and modular ranks, verifies the
coefficient identities and complementary zeros, checks the displayed form
and eigenvalue formulas on small instances, and checks the universal
inequalities. It does not assert that a full framed circuit for this
rejected candidate has been constructed. Written with OpenAI Codex assistance.
