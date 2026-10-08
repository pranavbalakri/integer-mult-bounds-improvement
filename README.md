# A conditional improvement beyond 2^-33

Starting from [CrocSwap/integer-mult-bounds at commit
6e564879f51ae16f23d392e9e196c605f36d90df](https://github.com/CrocSwap/integer-mult-bounds/tree/6e564879f51ae16f23d392e9e196c605f36d90df),
this proposed conditional research extension derives

$$
T(n)=O\!\left(n(\log n)^{1-\kappa}\right),\qquad
\kappa=\frac{125}{10^{12}}=1.25\times10^{-10}>2^{-33}.
$$

The previous displayed saving is `83/10^12`. The new saving is larger by
`125/83 ≈ 1.506`, or about **50.6%**. This concerns the asymptotic exponent,
not a measured runtime improvement. The original theorem and the repository's
compact-control, paired-bit, Gaussian and other retained extensions remain
assumptions. The additional argument was prepared with OpenAI Codex and has
not received independent mathematical review.

## What changes

Share the first and third stages' **complex** auxiliary roles. Two explicit
matchings handle side and central roles. Both preserve dirty-scratch
restoration, and their new binary frame transitions are nested.

The extra issue for complex frames is that nondegenerate binary residuals
must also be nonalternating. The proof exhibits a norm-one vector in every
new residual, obtaining the required orthonormal bases and translation-kernel
phase factorizations. Data endpoints and signed corrections stay unchanged.

At `h=25`, the new network has

| Quantity | Value |
| --- | ---: |
| Roles removed | 19,540,339,540,000 |
| Roles W | 39,105,013,080,000 |
| Total residual dimension s | 611,015,825,672,000,000 |
| Preserved absolute deficit Wm-s | 3,703,000,000 |
| Relative deficit | 1 / 165,005,625 |
| Certified complex exponent saving | 627 / 10^12 |

Every deleted role removes exactly `m=15625` residual dimensions. The
absolute deficit is preserved while its denominator decreases.

Keep the existing bit saving `296/10^11`. Choose

```
epsilon = 1999/10000     c = 1/4
beta = 1/1000           zeta = 1/10000
delta = 1/10^6          C1 = 49961/10000
lambda = 1-6265/10^13   lambda' = 1-626/10^12
```

All constraints hold strictly. The minimum assembly margin is
`G = 625687/(5*10^15) = 1.251374e-10`, leaving the positive gap
`G-kappa = 687/(5*10^15)`. The old unshared guard constants are retained
as conservative upper bounds; the time recurrence uses the new role counts.

## Included artifacts

This repository contains the three requested research artifacts:

1. This summary (`README.md`).
2. [Complete standalone proof note](notes/complex-reuse-note.tex).
3. [Independent manuscript patch](patches/complex-reuse-33.patch).

The original project's license and attribution are retained in
[LICENSE](LICENSE) and [NOTICE](NOTICE).

## Applying the patch

The patch applies directly to the unmodified pinned manuscript bundled in
the baseline repository. It includes the necessary prior refinements and is
an alternative to the earlier patches, not an additional patch to layer on top.
From this repository's root:

```sh
git clone https://github.com/CrocSwap/integer-mult-bounds.git source
git -C source checkout 6e564879f51ae16f23d392e9e196c605f36d90df
git -C source apply --check --directory=upstream ../patches/complex-reuse-33.patch
git -C source apply --directory=upstream ../patches/complex-reuse-33.patch
```

The note is a self-contained LaTeX file. The manuscript patch modifies the
baseline's `upstream/build` sources. This three-artifact bundle does not include
the additional certificate generators or test files mentioned in the note;
those were used in the development checkout. Running the baseline's test suite
alone does not reproduce the new extension's verification.

## Verification and scope

In the development checkout on October 7, 2026, all 168 tests and 18 independent
patch-application checks passed through `make verify`. Earlier tracked
certificates, patches, and upstream sources regenerated without changes.
The standalone proof note compiled successfully in the desktop LaTeX editor.
These are arithmetic, finite-identity, and integration checks, not formal
verification of the multiplication theorem or independent review of this proof.

The baseline research is by Douglas Colkitt, building on OpenAI's pinned
manuscript. This additional complex stage-sharing extension was prepared with
OpenAI Codex at the repository owner's request. It is not an official OpenAI
release or endorsement. See the proof note for assumptions and attribution.
