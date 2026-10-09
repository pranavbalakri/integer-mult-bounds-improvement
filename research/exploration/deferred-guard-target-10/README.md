# Linear guard for the actual signed PR97 complex word

The actual signed h24 word in PR97 supports the existing conditional guard

```
Δ = 128(d+1),       C1 = 1.
```

This conclusion is based on its complete event chronology and copied-centre
operations, not just on its sharing the old producer source. It permits the
same numerical precision choice `p=768 ceil(b^(43/20))` at `ε=71/100`, subject
to the inherited exact compiler, routing, normalization, and analytic
interfaces. This package does not prove those all-size interfaces or claim a
new multiplication exponent.

The external source is Zhihao Chen's PR97 integration, commit
`f5f9c56e637463cac1e300d1589ccf42838f688a`, including Swapnil Jain's signed
round-six complex producer. Its manifest still names four absent historical
logs. The audit records that discrepancy and does not report the upstream
top-level verifier as passing. The actual complex ledger is regenerated
independently of those missing logs in a disposable output directory.

## Exact event audit

`audit.py` invokes the pinned, credited `round6_complex_literal_ledger.py`,
then independently reads the resulting binary events and frame/copy metadata.
It checks ordinary common-frame gates, reflected frame continuity, and each
copy's actual live source and destination. The generated ledger has exactly
535744 events, including 24 copied-centre macros. Its summary is identical to
the pinned result after removing elapsed time.

Every ordinary scalar update and every update expanded from a copied-centre
macro is included. For each prefix, the audit tracks nonnegative upper bounds
on scalar coefficient row sums, starting every arbitrary data and scratch
input at norm one. It also tracks dyadic denominator exponents. A temporary
copy inherits its source row **at the actual copy time**; its paid frame
conversion preserves logical coefficients. The temporary is discarded after
its scatter, and the original arbitrary scratch input remains present.

The inverse event order reverses chronology, swaps data banks, and negates
each scalar coefficient. Reflected copied-centre macros create a fresh copy
from the reflected source, perform the inverse signed scatter, then discard
the private copy. They do not try to invert an erasure operation. Each
orientation has 374340 scalar updates, including the copied scatter.

The exact bounds are:

| Quantity | Forward | Inverse opposite shear |
|---|---:|---:|
| Maximum scalar-prefix row-sum bound | 1006577 | 1925603 |
| Maximum literal integer numerator product | 520904 | 1003352 |
| Maximum scalar denominator bits | 1 | 1 |
| Copied temporary streams | 24 | 24 |

The numerator is charged **before** division, including nonreduced literal
coefficients such as `2/2`. The largest coefficient is `19/2`. Exact dyadic
arithmetic allows cancelling a numerator's factor of two when determining
the denominator of the value after division; no approximate rounding is used.

Invocations within either tensor stage have disjoint banks. Uniform bounds
therefore compose once between stages. The copied endpoint addition costs
at most one further factor two. Consequently

```
G = 2 × 1006577 × 1925603 = 3876535381862 < 2^48,
B = 1 + 1 = 2 ≤ 8.
```

These equal the independently retained `semantic_guard.certify(24)` values.
The producer files are byte-identical, but the audit does not use that fact
as a substitute for event replay. It also reruns the pinned exact support
identity on all 4096576 ordered source/target coefficients.

## Signed phases and endpoint wrappers

For an e-axis address space, write

```
F_φ = H_e diag(i^φ) H_e,       φ : F₂^e → Z/4Z.
```

The inherited signed residual frames have this form. At a boundary between
completed recursive children, a block of a scalar prefix is

```
M[j,i] F_current[j] F_initial[i]⁻¹.
```

This identity is preserved by a common-frame scalar addition, a frame move,
and a copied temporary. A copy duplicates the logical row; a frame conversion
changes only its frame. Erasing a private temporary deletes that row. The
argument allows unequal initial frames and arbitrary scratch inputs. It also
allows descending copied-centre moves, since their inverse phase children
remain explicitly charged.

Every frame ratio is again Walsh-diagonal with fourth-root eigenvalues.
Its entries lie in `2^(-e) Z[i]`, and unitarity gives row norm at most
`2^(e/2)`. Thus a paused parent prefix consumes at most `e+B` fractional
bits and `e/2+log₂G` magnitude bits. The looser common charge `2e+64`
covers both quantities and the literal scalar intermediates.

The signed inverse wrapper does not invalidate this argument. On one axis,
with `Z=diag(1,−1)` and `X` the address swap,

```
C⁻¹ = −i Z C Z = X C.
```

The audit checks these identities exactly over Gaussian rationals. `Z` is
not being called a Walsh-diagonal frame: it is a norm-one, integral wrapper.
At a boundary after a completed inverse child, its complete returned operator
is `C⁻¹`, which is Walsh-diagonal. Inside the currently active child, its
precision is covered by that child's recursive guard. The input/output signs,
fourth-root factors, and address permutations preserve both numerical bounds.

The full signed endpoint needs both the inverse correction and the input
translation. Let `u=q_U(z)` and `w=wt(z)`, modulo four. Before correction,

```
X′ = −i^w y,
Y′ = i^(w−2u) x + i^(w−u) y.
```

Translate the input x by `P_U 1`. Its Fourier multiplier is `(-1)^u` because
`q_U(z) mod 2 = z·P_U 1`. Copy `X′`, apply the one charged inverse rank-one
phase child, and add that copy to `Y′`. The y term cancels; negate `X′` and
exchange banks. Every returned data role then has the same completed
`C_full` contract as each restored auxiliary. The translation, signs, and bank
exchange introduce no denominator loss. The endpoint addition is bounded
separately by the factor two above; no unsupported common-frame assertion is
needed for it.

Exact tests cover all sixteen pairs `(u,w)`. Replacing the inverse correction
by the forward phase fails eight cases, and omitting the input translation
fails eight. The finite tests accompany this general algebra; they are not a
substitute for the inherited signed residual compiler.

## Recursive and layer bounds

The regenerated histogram has largest child 552 at parent dimension 576.
For a divisible internal call on e axes, depth-first evaluation therefore gives

```
D(e) ≤ 2e + 64 + max_j D(t_j),       t_j ≤ (23/24)e.
```

Only the currently active child adds internal excess: completed children
already have their exact semantic returned format. All unused low positions
are mathematically zero; no return rounding, cancellation detection, or
dynamic format change is introduced.

With `C=51`,

```
C/24 ≥ 2 + 64/576,       C ≥ 8,
```

so the standard induction gives `D(e)≤51e`. Fewer than 576 remainder axes
are handled individually first. Their returned cost is at most their width,
and the elementary implementation uses the inherited eight-bit-per-axis
allowance, so the same induction covers all integer widths.

For one normalized layer, the inherited disjoint active-axis schedule charges
at most d bits for earlier completed pieces, `51d` for the active piece,
`8d` for individually handled preprocessing/tails, and `d+9` for outer
normalization. The total `61d+9` is less than `128(d+1)`.

The newly deferred **bit** word does not change this numerical calculation if
it supplies the inherited exact routing/permutation interface. Its internal
bit strings are not being treated as intermediate complex coefficients.
Its routing, product-row reservation, and finite-tape validity remain separate
proof obligations.

## Precision 43/20 with a finite numerical margin

Set `P=43/20`, `ε=71/100`, `b′=128 ceil(b^P)`, `p=6b′`, and
`d=floor(b^ε)`, for integer `b≥1`. Then

```
128(d+1) ≤ 256 b^ε ≤ 768 ceil(b^P) = p,
P−3ε = 1/50 > 0.
```

Under the inherited prime-selection band
`θ ≥ 1/(8 d L^ε)` with `L≤b`, let `α²=floor(b′/(8d))`.
Since `b^P/d≥1`, the floor estimate gives

```
α² ≥ 15 b^P/d,
α² θ ≥ (15/8) b^(P−3ε) ≥ 15/8 > 1,
2d α² ≤ b′/4.
```

Thus this guard and this particular Gaussian inequality require no new huge
numerical cutoff. The same argument covers cubic precision `P=3` with any
fixed `ε<1`. It does not supply the common cutoff for prime intervals, record
sizes, routing, or the remaining error/recovery estimates, and it is not a
practical running-time measurement.

## Reproduce and scope

```
python3 research/exploration/deferred-guard-target-10/audit.py --source /tmp/kappa-audit/pr97-full
```

The source root must contain the pinned PR97 integration. The checker rejects
`python -O`, does not alter downloaded source files, and writes `results.json`
in this directory. The receipt binds the regenerated event bytes and producer
hash, along with exact prefix bounds and the missing-log report.

This transfers the **conditional guard proof** to the actual signed word.
Whole-residual compiler correctness, exact fixed-tape wrappers, the active-axis
schedule, prime selection, normalization, and Gaussian recovery are inherited
interfaces. They have not been proved anew here. No historical certified file
is modified and no new multiplication theorem is claimed by this package.
