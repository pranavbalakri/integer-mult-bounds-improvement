# In-place producer: independent scalar, frame and guard audit

This audits the source-role borrowing construction proposed in
`../compact-adjacency-target-10/borrowed_sources.py`. It is positive local
progress, not a new certified multiplication exponent. The conservative child
histogram and local h18 frame word are checked below. All-size compiler
integration and inherited global interfaces remain separate obligations.

## Exact scalar identity

Let M be the existing invertible producer on its R physical slots. Partition
its input coordinates as (X,Z), where X are its v source pivots and Z are
its other R-v arbitrary dirty slots. The source pivots can be identified with
the v actual source data roles. The combined retained-centre and side read
map J satisfies

    JM(X,Z) = X + DZ.

This is checked on every source and dirty column in the supplied h=8 audit.
The chronological scalar word is

    Y <- Y-DZ;
    (X,Z) <- M(X,Z);
    Y <- Y+J(X,Z), with retained-centre reads before side reads;
    (X,Z) <- M^{-1}(X,Z).

It returns exactly (X,Y+X,Z). In characteristic two minus equals plus. Over
complex values the minus sign in the early cancellation is essential.
This identity holds for arbitrary dirty Z; it assumes no initialized original
scratch. The temporary retained-centre copies use the separately inherited
copy/erase interface.

## Common-frame chronology

At the low frame D0, all target Y roles and all Z roles can perform the dense
fixed cancellation DZ. X is excluded from this early map. Its initial frame
is its source triple line. The old producer assigns a common label to both
roles of every elementary mixer update. Replacing a source-pivot auxiliary
by its matching X role preserves those histories: both begin the labelled
part at the same triple line, and every later inclusion is unchanged.

After M, retained-centre values are read from copies moved down to D0 while
all Y roles remain there. Each Y then grows to its target-line complement for
the side reads. All original X/Z roles finally grow to the high frame D1,
where M^{-1} restores their logical values. The scalar action does not depend
on the contents of these roles while they carry the intermediate labels.

For stage two, use the **same chronological order**, replacing each elementary
scalar gate G by G^{-T}. If G1,...,GL compose F=GL...G1, this order composes
F^{-T}=GL^{-T}...G1^{-T}. Thus the completed second word is the opposite shear
X<-X-Y. Every ordinary gate still involves the same two roles at the same
common frame, so all its phase-frame incidences remain valid. New logical
supports do not need to equal the old producer supports; those supports were
used to find valid incidences, and common-frame conjugation then applies to
any scalar map on the aligned roles.

A copied centre read transposes by gathering the relevant Y combination into
one temporary at D0, applying the increasing child D0->Uc, and adding to the
original centre at Uc. This is the inverse-transpose of the logical read and
uses the same paid phase rank as the original decreasing copy. In particular
one should not charge a separate phase child for every target in the gather.

## Two-stage inclusions

Let a and b be the two triple lines, and let F_h be the one-factor space.
Stage one uses local U lifted to U tensor <b>. Its X/Y outgoing frames are
F_h tensor <b> and a-perp tensor <b>.

For stage two put

    D0 = a-perp tensor F_h,
    D_U = D0 orthogonal-sum (<a> tensor U).

Before any middle gate, explicitly place every borrowed X at D_<b>.
The growth from F_h tensor <b> has the inherited residual
a-perp tensor b-perp, of rank (h-1)^2. Subsequent growth from <b> to
the first actual local label U is exactly the old producer edge. We retain
this intermediate frame even if the first gate already has a larger label;
no merged-edge improvement is claimed. The first Y frame is D0, which
contains a-perp tensor <b>.
The final frames are the whole h^2-dimensional space on X and
(a tensor b)-perp on Y. The two scalar shears still give (-Y,X+Y), so the
inherited one-child endpoint correction applies unchanged.

Explicitly grow **every** remaining Z role to D0 at the beginning of stage
two, even if its D column is zero. This guarantees that every initial outer
child has width m-h; omitting that harmless intermediate frame can instead
create a single larger child and requires a different histogram and guard
constant. Include these explicit edges in the actual histogram.

## Exact local checks

`inplace_complex_audit.py 8` checks the entire 957-dimensional dirty scalar
basis with exact integer numerator arithmetic, including all half-coefficients.
It checks both orientations and confirms that all 845 remaining auxiliaries
restore. It separately checks every actual middle label inclusion, every side
residual, every retained-copy residual and final growth against the existing
binary nondegenerate/nonalternating checker. The receipt is
`inplace-complex-results.json`. Omitting a nonzero early cancellation leaves
an explicit nonzero dirty coefficient.

The larger local audit is also efficient: run
`python3 inplace_complex_audit.py --frames-only 18`. It checks 46,096 distinct
exact binary frame transitions, 11,437 middle common-frame gates, all 8,160
side incidences, 18 retained-copy frames, 18,799 final growths, and all 816
triple perpendicular spaces. The receipt is `inplace-h18-frame-results.json`.
No full 30-million-role tensor matrix is materialized. The tensor entry
and outer residuals have the exact product forms
`a-perp tensor b-perp`, `a-perp tensor F_h`, and `F_h tensor b-perp`.
Each factor has a checked orthonormal basis, so the tensor products do too.
All local residuals tensor with an odd-norm line and preserve their Gram forms.

## Conservative chronological histogram

`inplace_complex_profile.py` reconstructs every local role's growth directly
from the mixer chronology, followed by retained copies, side reads and final
cleanup. It includes the explicit source-line and D0 vertices above, both
interstage data-entry edges, all remaining auxiliary outer edges, and the
unchanged endpoint correction. It then independently verifies agreement with
subtracting 2N children of each width m-h, h-1 and 1 from the old histogram.
The width changes from `2N+2vR` to `2vR`; the raw deficit remains unchanged.

| h | New width W | Raw rank deficit | Exact certified complex saving |
|---|---:|---:|---:|
| 16 | 13,942,880 | 43,680 | 4113971/100000000000 |
| 18 | 30,679,968 | 164,832 | 15363143/250000000000 |

The positive rational moment gaps and complete histograms are in
`inplace-profile-results.json`. At h18 the local saving is 0.000061452572.
This is a separately checked new complex profile; it is not automatically a
new multiplication exponent or a replacement for the unchanged assembly.

## Explicit guard survives

`inplace_complex_guard.py` bounds every scalar prefix, including the dense
cancellation, without storing its dense matrix. The producer M uses positive
additions. Therefore |D| is entrywise bounded by |J|M restricted to the dirty
columns. Row sums of this bound control the forward early cancellation;
column sums control the inverse-transpose cancellation. Subsequent absolute
updates follow the exact chronological mixer and read order. A factor two
covers forming integer numerators before halving.

The resulting bounds, including both stages and the endpoint addition, are:

| h | Forward envelope | Inverse-transpose envelope | Combined envelope | Scalar denominator bits |
|---|---:|---:|---:|---:|
| 8 | 948 | 3,168 | 6,006,528 | 2 |
| 16 | 27,370 | 68,264 | 3,736,771,360 | 2 |
| 18 | 46,256 | 123,054 | 11,383,971,648 | 2 |

In particular the h18 combined envelope is below 2^34. With the explicit
D0 growth and inherited whole-residual integer widths, the previous loose
prefix charge 2e+64, recursive coefficient 40, and layer guard 128(d+1)
remain sufficient. At h16 the recursive coefficient remains 36. These are
conditional guard checks for the proposed new word; they do not replace its
inherited recursive-width or global-transfer hypotheses.

Written with OpenAI Codex assistance. Existing producer, phase and copy
interfaces retain the contributor credits in the main repository. The
source-role borrowing proposal is coordinated with the bit construction
linked above; this note records its independent complex audit.
