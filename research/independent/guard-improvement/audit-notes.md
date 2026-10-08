# Independent audit of the semantic guard refinement

This audit concerns the exact arithmetic inside the inherited two-stage complex
phase primitive. It does not independently validate the full multiplication
algorithm, its tape layouts, prime selection, or Gaussian inversion interface.

## Common-frame prefix argument

Fix a recursive node with `e` active binary directions. At a boundary between
completed recursive children, each block of the prefix operator has the form

`M[j,i] D_current[j] D_initial[i]^-1`,

where `M` is the prefix of the scalar role program. This follows by induction:
scalar gates meet equal frames, a frame edge changes only the corresponding
`D_current`, and a copied zero temporary adds a duplicate scalar row in the
source's frame. Erasing a temporary replaces that row by zero. Neither copy nor
erasure needs to be invertible for this identity to hold.

All the `D` operators are diagonal in the same Walsh basis, with eigenvalues in
`{1,i,-1,-i}`. This is also true partway through the list of *completed* children
implementing one edge. Consequently `D_current[j] D_initial[i]^-1` is again one
such operator, rather than a product that needs two independent worst-case
precision allowances. Its matrix entries are Gaussian integers divided by
`2^e`. Its infinity norm is at most `2^(e/2)` by unitarity and Cauchy--Schwarz,
and in particular at most `2^e`.

If every scalar-prefix row sum is at most `G` and its scalar coefficients have
at most `B` denominator bits, every node-boundary value has at most `e+B`
additional denominator bits and magnitude amplification at most `2^e G`.
Only one recursive child is active in the depth-first schedule. The guard thus
charges the boundary allowance plus the peak within that one child; it does
not add the implementation depths of all already-completed children.

This is an exact algebraic statement. With one sufficiently wide fixed-point
format, the low bits that cancel at a completed-child boundary are already
mathematically zero. No approximate truncation, detection of cancellation, or
resizing operation is needed. Input/output signs and the inverse-child sign
wrappers preserve magnitude and denominator bounds. Bit-routing subroutines
are treated as exact permutations of encodings at their interfaces, as in the
inherited proof; their temporary bit strings are not used as coefficients.

## Independent finite checks

From the `research` directory, run `python3 independent/guard-improvement/audit.py`.

* At `h=16`, the script checks all **313,600 ordered triple pairs**. The retained
  scatter plus the signed half-weighted pieces equals the identity entry by
  entry. Flipping one emitted piece's sign makes the check fail.
* An independent compiler reconstructs the signed elementary updates in the
  actual ten-part dirty-scratch word. At `h=16` it has 12,449 scratch roles and
  93,724 rational updates per invocation. Positive-envelope propagation bounds
  the forward and inverse rational-update boundaries by 135,386 and 256,057,
  respectively. These particular small bounds describe rational-update
  boundaries; the production `semantic_guard.py` additionally charges the
  integer products used while realizing coefficients such as `11/2` and uses
  the conservative common allowance `J=64`.
* At `h=8`, the script runs **all 112 invocations** of both scalar stages over a
  prime field, with 3,136 data pairs and **100,912 arbitrary dirty scratch
  values**. Every scratch value is restored, and the exact data map is
  `(-y,x+y)`. The scalar part of the charged endpoint correction restores the
  original role order. Copying the retained totals too early fails a negative
  control.

The scalar check uses literal copied retained values for the second scatter.
The copy is a zero-initialized temporary; it does not consume an arbitrary
scratch input. Its recursive phase edge remains separately charged in the
existing child histogram. Excluding it from the logical role count `W` is
therefore consistent with the fixed-tape, irreversible machine model.

## Endpoint correction pinned to its primary source

The original construction is Aurel Prosz's (Paureel) two-stage extension, pinned
here to commit `c82d09eb4781b68435e3fa3eb6beea72fa900fab`:

* [Two-stage construction](https://github.com/Paureel/integer-mult-bounds/blob/c82d09eb4781b68435e3fa3eb6beea72fa900fab/notes/two-stage-construction.tex)
* [Two-stage complex endpoint correction](https://github.com/Paureel/integer-mult-bounds/blob/c82d09eb4781b68435e3fa3eb6beea72fa900fab/notes/two-stage-phase-transfer.tex)

Write `P=C_<w>`, `F=C_full`, and `E=F P^-1`. The two framed stages yield

`A=-F y`, `B=F P^-2 x + E y`.

The correction copies `A`, applies the one charged inverse phase child `P^-1`
to that copy, and adds it into `B`. Because `P^-1 F=E`, the `y` term cancels.
This uses one scalar addition and no scalar halving. The remaining operations
are input/output diagonal signs, fourth-root phases, a negation, and bank
reassembly. They preserve absolute values and Gaussian denominator exponents.
The scalar addition can be allowed one extra magnitude bit, well within the
production allowance.

The diagonal signs are not Walsh-diagonal operators. They are factored as
integral norm-one input/output wrappers in the guard argument, rather than
incorrectly described as members of the common Walsh-diagonal frame family.
