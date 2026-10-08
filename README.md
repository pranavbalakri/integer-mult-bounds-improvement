# A conditional improvement beyond 2^-31

Starting from [CrocSwap/integer-mult-bounds at commit
6e564879f51ae16f23d392e9e196c605f36d90df](https://github.com/CrocSwap/integer-mult-bounds/tree/6e564879f51ae16f23d392e9e196c605f36d90df),
this proposed conditional research extension derives

$$
T(n)=O\!\left(n(\log n)^{1-\kappa}\right),\qquad
\kappa=\frac{59}{10^{11}}=5.9\times10^{-10}>2^{-31}.
$$

The previous extension gave `125/10^12 = 1.25e-10 > 2^-33`.
The new exponent saving is **4.72 times larger**, and approximately **7.108
times** the original repository's `83/10^12` witness. These are asymptotic
exponent comparisons, not measured runtime improvements.

The original theorem and the repository's compact-control, paired-bit,
Gaussian and other retained extensions remain assumptions. The additional
argument was prepared with OpenAI Codex and has not received independent
mathematical review.

## What changes

Replace individual complex side corrections with **signed rectangle sums**.
Disjoint triple pairs receive coefficient `+1/2`; intersection-two pairs
receive `-1/2`. Reversible mixers share each rectangle's intermediate sum
while restoring arbitrary initial scratch values.

Four explicit binary frame choices handle singleton source families,
singleton target families, disjoint families, and intersection-two families.
Their residuals are nondegenerate and nonalternating in both circuit
directions, supplying the orthonormal bases needed for the phase interface.
First/third-stage auxiliary banks are then shared as in the preceding proof.

At `h=25`, the construction has:

| Quantity | Value |
| --- | ---: |
| Rectangles per invocation | 33,951 |
| Side roles per invocation | 455,434 |
| Total roles W | 4,843,100,800,000 |
| Total residual dimension s | 75,673,446,297,000,000 |
| Preserved absolute deficit Wm-s | 3,703,000,000 |
| Relative deficit | 7 / 143,050,000 |
| Certified complex exponent saving | 5 / 10^9 |

Keep the bit saving `296/10^11` and choose

```
epsilon = 1999/10000     c = 1
beta = 1/10             zeta = 1/10000
delta = 1/10^6          C1 = 46001/10000
lambda = 1-2959/10^12   lambda' = 1-2958/10^12
```

Now `sigma < tau`, so the internal layer exponent is `tau`. All constraints
hold strictly. The minimum assembly margin is
`G = 2956521/(5*10^15) = 5.913042e-10`, leaving the positive gap
`G-kappa = 6521/(5*10^15)`. Conservative unshared guard constants still
cover the new scalar gates; the time recurrence uses the new role counts.

The fixed bit exponent and retained Gaussian/assembly inequalities give the
scoped ceiling `kappa < (1-tau)/5 = 5.92e-10 < 2^-30`. The new witness
reaches about 99.66% of it. This is not a ceiling for other bit primitives
or multiplication algorithms.

## Included artifacts

This repository contains the three requested research artifacts:

1. This summary (`README.md`).
2. [Complete standalone proof note](notes/complex-reuse-note.tex).
3. [Independent manuscript patch](patches/complex-rectangles-31.patch).

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
git -C source apply --check --directory=upstream ../patches/complex-rectangles-31.patch
git -C source apply --directory=upstream ../patches/complex-rectangles-31.patch
```

The note is a self-contained LaTeX file. The manuscript patch modifies the
baseline's `upstream/build` sources. This three-artifact bundle does not include
the additional certificate generators or test files mentioned in the note;
those were used in the development checkout. Running the baseline's test suite
alone does not reproduce the new extension's verification.

## Verification and scope

In the development checkout on October 7, 2026, the certificate generators, all 175 tests,
and all 19 independent patch-application checks passed. Earlier tracked
certificates, patches, and upstream sources regenerated without changes.
The standalone proof note compiled successfully in the desktop LaTeX editor.
These are arithmetic, finite-identity, and integration checks, not formal
verification of the multiplication theorem or independent review of this proof.

The baseline research is by Douglas Colkitt, building on OpenAI's pinned
manuscript. The additional complex stage-sharing and signed-rectangle extensions were prepared with
OpenAI Codex at the repository owner's request. It is not an official OpenAI
release or endorsement. See the proof note for assumptions and attribution.
