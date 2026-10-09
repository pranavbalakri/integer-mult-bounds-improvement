# Dirty references and auxiliary unmixing at the sink

This is a checked scalar prototype and a proposed common-frame schedule. It
changes the arbitrary-source word, rather than merely changing the endpoints
of the existing borrowed word. The complete scalar identity, including all
reference dirt, passes finite checks. The exact fixed-basis first-stage
profiles are now in [PROFILE.md](PROFILE.md): both full reference versions
lose to the checked borrowed-source construction. A stronger relaxed check
also rejects improving that bound with point-star cleanup alone. Earlier
cleanup at nested original frames remains a separate possibility.

## Exact scalar identity

Use an **unborrowed** positive producer. Its middle word `A` acts invertibly
only on its `R` auxiliary roles; it does not overwrite source data `X`. Let
`P` inject the data into its leaf pivots, and `J` include retained-centre and
side scatters. The inherited scalar identity is

```
J A P = I.
```

Introduce reference roles indexed by `(c,S)` with `c` in the source triple
`S`. A reference starts with arbitrary public value `B_(c,S)`, receives `X_S`
at the source-line frame, and is preserved until cleanup. With all three
references per triple, the additional bank has `3v` roles.

A smaller variant uses `2v` reference roles: choose one primary point `c(S)`
for every triple (the checker uses its least point), and use the unchanged
source data `X_S` itself as its reference for that one point. The other two
references remain arbitrary auxiliary inputs. Every source injection and
reference creation is completed before any primary `X` leaves its source line.

Write `P_B X` for the vector copied into the auxiliary references. For each
producer role, let its ideal clean-source row be the corresponding row of
`T=A P`. Choose a common point `c` in its last producer frame. Its support is
a sum of triples containing `c`, so its clean source contribution is the sum
of their references. In matrix form this gives

```
A P = K_B P_B + K_X,
```

where `K_X` is zero in the `3v` version and contains the direct primary-`X`
terms in the `2v` version.

The chronological scalar word is:

1. At the base frame, apply the early correction `Y += J A R` using the original
   arbitrary producer dirt. Reference dirt does not occur in this correction.
2. Copy `P_B X` into `B`, inject `P X` into the producer leaves, and run `A`.
3. Scatter retained centres and then sides, with the inherited centre-copy
   transfers. This makes `Y = Y_original + X` because the two copies of
   `J A R_original` cancel.
4. After every producer role's final use, add its reference combination:
   `R += K_B B + K_X X`.
5. At the full local-stage frame, clear the reference signals with
   `B += P_B X`. The source data is still unchanged.

After step 4,

```
R = A R_original + K_B B_original;
B = B_original + P_B X.
```

After step 5, the complete auxiliary map is therefore

```
(B,R) -> (B, A R + K_B B).
```

It is invertible and independent of **both** data banks. In particular, the
reference dirt is not assumed zero. At the common auxiliary sink frame,
apply `R += K_B B` and then `A^-1`. This restores every arbitrary producer and
reference input. Omitting the `K_B B` unmixing fails the check.

The public scalar map is exactly

```
(X,Y,R,B) -> (X,Y+X,R,B).
```

Transposing every elementary CNOT in the **same chronological order** yields
the opposite shear and again restores all auxiliaries. Both orientations are
checked on the full public input basis.

## Why the sink unmixing can be pointwise

Use auxiliary input offset `0` and output frame `I` for every role in the
block being unmixed. They therefore all represent their own full transform
at the same sink frame. The map being undone acts only on original auxiliary
inputs with this same offset, and has no remaining `X` or `Y` component.
Consequently its inverse is a pointwise role operation there.

This argument does not permit free mixing between roles with different input
frame offsets. No shift by minus an exterior projector is assumed here. The
inverse `A^-1` has moved to the common global sink; it is not performed at the
local full-stage frame for every producer role.

## Actual cleanup frames and required paths

In the first tensor stage, with fixed other-factor line `Q`, a producer's
local frame `U` lifts to `U tensor Q`. A point-`c` reference is created at
`P_S tensor Q`, then grows to `U_c tensor Q`, where it serves only cleanups
anchored at this same point. It subsequently grows to the local full frame
`F_h tensor Q` for signal erasure, and then to the global sink. Its complete
path must be charged. The `2v` version also grows each primary source `X_S`
through its **one** chosen `U_c tensor Q` before the local full frame.

A nonterminal producer role with last label `A_last` can be cleaned early only
if that label is contained in a common point-star span `U_c`. The checker
verifies the support inclusion explicitly. A role without such an anchor
must remain for full-stage cleanup. The cleanup itself occurs at **`U_c`**,
of dimension `h-1`: it generally requires growing from `A_last` to `U_c`.
It is incorrect to use the original `A_last` as the cleanup frame merely
because its source support consists of triples containing `c`.

After cleanup, that role can go directly from `U_c tensor Q` to the global
sink. Its exit residual is

```
I - U_c tensor Q,
```

of rank `m-h+1`. This can be profiled as one complementary residual before
the pointwise sink unmixing. It is only one rank larger than the exterior
rank `m-h`; the extra reference paths and the growth to `U_c` must also be
counted. A gain cannot be inferred from the size of this exit alone.

Side-output roles have already reached their target hyperplane when read, so
the prototype keeps their cleanup at the local full-stage frame. Retained
centre roles may clean at their point-star frame after their centre copy has
been read. Every cleanup is after the role's final producer and scatter use.

The second tensor stage requires a separate account. With input offset `0`,
its early dirty correction is at the nonzero base frame
`E_i=P_i^perp tensor F_b`; that entrance already pays the exterior. The
first-stage retirement merger does **not** automatically give a corresponding
second-stage saving. A valid next comparison can use the new word only in the
first stage, retaining the old second-stage implementation. A full two-stage
variant needs its actual matrix residuals profiled, including reference
creation and the primary-source data paths.

## Finite checks and remaining work

`scalar_word.py` checks `h=5,6,7,8`, both reference counts, and:

- the literal unborrowed producer and scatter identity `JAP=I`;
- every producer role's ideal source row and its containment in its last label;
- a common-point anchor for every early cleanup (with a full-stage fallback);
- all public input basis vectors, including arbitrary references;
- the same-chronology transposed word;
- data independence and invertibility of the intermediate auxiliary block;
- rejection when the final reference-dirt unmixing is omitted.

The first-stage frame matrices, all changed reference/data paths, and complete
moments have now been checked against the pinned PR84 baseline; see
[PROFILE.md](PROFILE.md). The result is negative for the full reference
versions, so no updated multiplication exponent is claimed. A different
cleanup schedule would still require its own full path account, literal-word
integration, pointwise costs, and inherited global interfaces.

```sh
python3 research/exploration/sink-unmixing-target-10/scalar_word.py
```

The checker uses standard Python and rejects `-O`. The historical/certified
files are unchanged.
