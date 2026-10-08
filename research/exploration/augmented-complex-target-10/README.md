# An actual complex side producer for the augmented motif

This is a checked finite scalar word, frame schedule, child histogram and
precision-envelope calculation for the augmented triple/difference motif.
At h=19 its complex moment root lies between **0.0000548119 and
0.0000548120**. This is below the existing h=18 complex certificate
0.000059292472, so **no improvement to the published multiplication bound is
claimed**. The extra side-production work outweighs this motif's larger
rank deficit in the present construction.

Reproduce the saved receipt from the repository root:

```sh
python3 research/exploration/augmented-complex-target-10/check_word.py 7 11 19 --dirty-small --output research/exploration/augmented-complex-target-10/results.json
```

Only Python's standard library is required. Assertions must remain enabled.
The producer imports the existing complex triple-exclusion DAG, binary frame
checker, augmented motif definitions, and independent rational moment
enclosure. It constructs its own new signed pieces, role compiler, scalar
word, chronological trace and histogram.

## Labels and the exact scalar identity

For odd h, the indices comprise ordinary triples S and differences (i,j),
i<j. Their rational labels x and binary labels b are

    triple S:       x=1_S,       b=1_S,               t=1;
    difference ij:  x=e_i-e_j,   b=1+e_i+e_j,         t=0.

There are v=binomial(h+1,3) indices. All binary labels have odd norm. The
rational fitting matrix is

    A(s,t) = (x_s dot x_t - t_s t_t)/2.

It has diagonal one and is zero whenever the corresponding distinct binary
labels are not orthogonal. Its complementary side matrix is I-A. The exact
motif and central-factorization proof is in
[the augmented fitting note](../asymmetric-fit-target-10/README.md).

The retained centres are

    z_i = sum_t (t_t-x_t[i]) input_t,

and their coefficient at target s is

    -x_s[i]/2 + t_s/(h-3).

Each z_i uses the coordinate hyperplane i=0 as its binary frame. Its paid
copy rank is h-1, giving local centre cost h(h-1). For h=7,11,19 the central
coefficients are dyadic with denominators 4,8,16, respectively.

The producer constructs I-A from four kinds of contributions:

1. The ordinary triple-to-triple side uses the existing triple-exclusion
   and pair-star DAG. Build that DAG at even arity h+1, set every input
   containing the extra point to zero, and retain only actual triple targets.
2. For a triple target S, add -1/2 times each signed incident-edge sum at
   i in S whose other endpoint is outside S. Internal edges never enter
   these sums.
3. For a difference target (i,j), add -1/2 times the triple sum containing
   i but not j, and +1/2 times the triple sum containing j but not i.
4. At that difference target, also add -1/2 times the signed incident-edge
   sum at i excluding j, and +1/2 times the corresponding sum at j excluding i.

An incident edge at i has value input_(i,j) for i<j and -input_(j,i) for j<i.
Interval sums provide the cut and deleted-edge pieces. For the deleted-triple
pieces, fix i and regard triples containing i as pairs on the other points.
For each omitted point j, add the pair total below j, the pair total above j,
and the crossing terms. Shared prefixes, suffixes and row suffixes give an
O(h^3) construction of this new part over all i.

Every DAG value is represented exactly by disjoint positive and negative
support masks. The checker verifies the compiled clean producer on every
source column, every retained centre, and every entry of the resulting
identity matrix. At h=19 this checks all 1,299,600 scalar matrix entries.

## Frames and the odd-arity repair

Incident-edge labels at a fixed i are orthonormal when h is odd. Therefore
an interval sum can use the span of its incident labels, with orthonormal
residuals along its growth. A side piece whose final residual is alternating
is split into its proper summands before being read; the checker verifies
the resulting pieces rather than assuming that every complete cut is valid.

The deleted-triple DAG needs a different family. Let I be the set of omitted
target points served by a node at fixed i. Use

    K_I = (span{ b_(i,j) : j in I })^perp.

The incident labels are orthonormal, so K_I is nondegenerate, has dimension
h-|I|, and nested changes K_I -> K_J for J contained in I have orthonormal
residuals. K_I is also nonalternating: its canonical vector is
1+sum_(j in I) b_(i,j), which is nonzero. The node's triple inputs belong to
this frame because every served omitted point is absent from each input.
Single-input nodes keep their original line. All entry, exit, argument and
side-read residuals are checked explicitly.

Filtering the old even-arity triple DAG introduces another issue: an odd
pair-star span can have an alternating residual into its exact coordinate
cover. Such a transition is not accepted as a whole child. At that consumer,
the compiler instead reads the star's proper summands separately, preserving
the scalar value and destination frame. There are 672 such argument
expansions at h=19. All final transitions pass the exact binary test for
nondegeneracy and nonalternation. The h=19 chronological construction checks
99,449 distinct frame transitions before the full scalar trace is added;
the saved receipt records its resulting checked counts.

## A complete word with arbitrary original scratch

Let M be the compiled invertible producer on R role slots. Its elementary
operations are signed pivot additions and adding a pivot into its other
output slots. Those original slots may contain arbitrary dirty values.
Let S inject the v input values into its source pivots and J
be the combined central and side read. The clean scalar checks prove

    J M S = I.

For arbitrary original scratch Z, the forward chronological word is

    Z <- M Z;             Y <- Y-JZ;      Z <- M^-1 Z;
    Z <- Z+S X;
    Z <- M Z;             Y <- Y+JZ;      Z <- M^-1 Z;
    Z <- Z-S X.

It returns exactly (X,Y+X,Z). The first three operations occur at the low
frame. Source injection occurs at each source line. The second M uses the
checked node labels; central reads occur from copied values moved to the low
frame while Y stays low; side reads occur after Y reaches its target
perpendicular. Cleanup and the final source subtraction occur at the full
frame. Centre copies are fresh temporaries under the inherited copy/erase
interface, rather than initialized original scratch.

At h=7, an additional check evaluates the complete forward word and its
inverse-transpose on every source and dirty column simultaneously. It uses
balanced positional encoding of the basis into large integers, with an
independently computed coefficient bound smaller than half the radix; hence
zero equality cannot conceal a carry between basis columns. The basis
dimension is 1,283. Both shears and scratch restoration pass, and omitting
the first cancellation leaves a nonzero dirty contribution. Larger instances
use the proved elementary-word identity, the exact clean column checks and
the checked compiled inverse, rather than materializing a dense dirty matrix.

## The second stage uses the same chronological frame path

There is no stage-two bank relabeling and no reversal of the local frame
sequence. Replace every chronological elementary scalar gate G by G^-T in
the **same order**. If G1,...,GL compose F0=GL...G1, these new gates compose
F0^-T=GL^-T...G1^-T. Thus the second local word performs the opposite shear
X<-X-Y and restores its scratch. Each ordinary gate has the same physical
role incidence and common label as before; only its scalar direction and
coefficient change.

A central read requires its explicit transposed bridge. The forward read
copies a centre at U, moves the copy U->0, scatters into Y, and erases it.
The inverse-transpose gathers the relevant Y combination into a fresh
temporary at 0, moves that temporary 0->U, adds it into the original centre,
and erases it. Its paid phase rank is the same dim(U), once per centre,
not once per target. This step inherits the same copy/erase interface.

For the tensor lift, write a and b for the two odd-norm input lines in
F_2^h. Stage one, repeated for each fixed b, uses U tensor b. Its data
endpoints are F_h tensor b on X and a_perp tensor b on Y. Stage two,
repeated for each fixed a, uses

    D_U = (a_perp tensor F_h) orthogonal-sum (a tensor U).

Before stage-two source-line gates, X grows from F_h tensor b to D_b;
before its low-frame operations, Y grows from a_perp tensor b to D_0.
Both joins have residual a_perp tensor b_perp, rank (h-1)^2. Then both data
roles follow the increasing local path of rank h-1. Their final frames are
the whole h^2-dimensional space on X and (a tensor b)^perp on Y.

Starting with logical coordinates (P^-1 x,y), the two shears give
(-y,P^-1 x+y). Hence the physical outputs are exactly

    A=-F y,       B=F P^-2 x + F P^-1 y,

the retained endpoint identity. Apply its one rank-one inverse phase child
to a copy of A, add it into B, and use the diagonal signs and final bank
exchange. Every augmented label is odd, so P^2 is still the required
translation. Some products of label weights are 3 rather than 1 modulo four;
the signed rank-one phase and fourth-root wrappers account for that sign.
No additional address permutation is declared free.

The checker records every scalar gate at its chronological frame and
requires that this trace reproduce the exact elementary word. It verifies
the original-role growth histogram and the centre bridges from that trace.
The inverse-transpose has identical ordinary incidences and the explicitly
reversed temporary bridges above. This is why the two stages have equal
local histograms; equality is not inferred from dimensions alone.

## Counts and the finite moment

Let H be the traced ordinary auxiliary histogram for one local invocation,
excluding centre copies. Its weighted sum is hR. Put m=h^2 and N=v^2.
The two stages have 2v invocations, with one rank-(m-h) outer auxiliary
edge per original auxiliary per invocation. The complete child histogram is

    2v H + 2vR [m-h]
       + 2N [(h-1)^2] + 4N [h-1] + N [1]
       + 2vh [h-1].

The last term is the copied centres and the singleton is the endpoint
correction. With W=2N+2vR, its total rank is

    Wm-N+2vh(h-1).

| h | v | R | W | Rank deficit |
|---|---:|---:|---:|---:|
| 7 | 56 | 1,171 | 137,424 | -1,568 |
| 11 | 220 | 6,684 | 3,037,760 | 0 |
| 19 | 1,140 | 44,672 | 104,451,360 | 519,840 |

The first two instances supply scalar and frame validation, not positive
savings. At h=19, all child ranks are at most 342, and the exact rational
log/exp enclosure accepts a=548119/10^10 and rejects the next 10^-10 grid
point. The exact gaps and complete histogram appear in `results.json`.

## Semantic guard and remaining integration

The full scalar prefix, including low-frame cancellation and both cleanup
words, has a forward envelope 1,631,376 and inverse-transpose envelope
5,045,272 at h=19. These already allow forming integer numerators before
division by 16. Including both stages and the endpoint addition gives

    G = 16,461,471,308,544 < 2^44,

with four scalar denominator bits per stage. The prior loose prefix charge
2e+64 remains sufficient for this finite word. With contraction rho=18/19,
the returned-format argument gives recursive guard coefficient 42, since

    42(1-rho) >= 2 + 64/361.

This fits the coefficient budget of the inherited 128(d+1) layer allowance.
That numerical check remains conditional on the inherited integer-width,
copy, outer-phase and layer-assembly interfaces. The new finite word has not
been integrated into a complete all-size compiler or a new multiplication
theorem. In particular, this package does not claim executable integer
multiplication timings, a practical algorithm, or the requested 2^-10 bound.

The existing triple DAG and phase interfaces retain their original authors'
credits. The added signed pieces, frame repairs, role-word audit and exact
certificates were written with OpenAI Codex assistance.
