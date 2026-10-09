# Coordinated caps of 20 and 19 on deferred entrances

This package constructs smaller entrance frames for the pinned deferred
round-seven word, using one shared replacement for every connected group of
readouts that share targets. It gives exact finite frame certificates for
caps of 20 and 19. A concrete pair of unchanged predecessor frames obstructs
a cap of 18 in the same fixed-order, fixed-lower-frame construction.

The scalar operations, public roles, lifted producer frames, leaf-copy
frames, retained centres, and data entrances are unchanged. These certificates
do not establish a new multiplication exponent or a practical crossover.
They extend the [cap-21 construction](../deferred-truncation-target-10/README.md);
the separate interface audit accounts for recursion and numerical tradeoffs.

## Credited source and exact scope

The source is [Swapnil Jain's round-seven word](https://github.com/Swapnil-jain/integer-mult-kappa/tree/741e7aa078392553815df7926ee17ac5e25a8c38),
imported into [CrocSwap PR97](https://github.com/CrocSwap/integer-mult-bounds/pull/97).
The original deferred frame construction and common axis bases are inherited
from that source; Zhihao Chen's PR97 contribution supplies its signed and
reflected integration. This package supplies the componentwise cap choices,
their exact certificates, and the scoped cap-18 obstruction.

Both input certificates and the imported reader, linear algebra, and frame
modules are checked against the SHA-256 pins in the adjacent truncation
package. The original scalar word and its original frame inclusions remain
inherited obligations. No upstream file is edited.

All supplied frame bases are integer matrices representing rational row
spaces in dimension 23, with Gram form `G=I-J/9`. Containments are verified
over the integers against exact rational annihilators. A nonzero Gram
determinant modulo `2^31-1` proves exact nondegeneracy over the rationals;
small modular ranks are never used as rational rank upper bounds.

## Why the caps must be coordinated

For a chosen cap `k`, keep every old readout frame of dimension at most `k`
unchanged. Connect two old higher-dimensional readouts when they contribute
to a common target under the **integer** garbage-coefficient matrix. Use a
single cap `C` of dimension `k` throughout each connected component.

Equal-dimensional frames on a nested target chain must coincide. Treating
each high readout independently would therefore be insufficient. The
integer target set contains the parity target set, so the common cap serves
both the integer and binary versions.

For each component the checker chooses:

- `H`, a least-dimensional original high frame, and verifies **H lies in
  every original high frame in the component**;
- `B`, a largest unchanged predecessor frame, and verifies **every unchanged
  predecessor frame for every target in the component lies in B**;
- an explicit nondegenerate `k`-plane `C` with `B ⊆ C ⊆ H`.

These facts are checked on the actual matrices, not inferred from dimension
order. When `dim B=k`, use `C=B`. Otherwise append explicit integer linear
combinations of the rows of `H`, accepting only candidates with nonzero Gram
determinant. Every selected combination is stored in the receipt. Each cap
is also checked against every corresponding target-hyperplane normal.

Consequently every new frame is contained in its original frame, every
required unchanged predecessor is contained in the cap, and all shared
target chains remain nested in the original readout order.

## The two constructions

| Cap | Changed readout slots | Connected components | Existing caps reused | Components extended |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 3,356 | 1,651 | 1,188 | 463 |
| 19 | 4,663 | 1,661 | 1,660 | 1 |

For cap 20, 462 components start from a 19-dimensional predecessor and need
one added vector. One component has no nonzero predecessor and receives an
explicit 20-plane inside its old high frame. All old dimension-21 and
dimension-22 readouts happen to have a single integer target in this witness.

For cap 19, 100 changed readouts have two integer targets, so the connected-
component construction matters. Of the 1,661 components, 1,660 reuse an
existing 19-dimensional predecessor frame. The one component with no
predecessor receives an explicit 19-plane inside its old high frame.

The largest auxiliary block changes from `527/529` in the original word to
`523/529` for cap 20 and `521/529` for cap 19. Public width does not decrease.
Smaller maximum children can reduce a recursive row-stock bound, but their
effect on the finite moment and global parameters must also be accounted for.

## Explicit first-edge implementation

Let `F0(s)` be an auxiliary role's unchanged first producer frame. The
certificate retains the original readout frame as an **identity vertex**:

```
new cap C -> old sigma_s -> unchanged F0(s).
```

No scalar operation or new role occurs at the identity vertex. The new first
piece is checked explicitly; the original `sigma_s -> F0(s)` edge retains its
original certificate. This avoids pretending that the merged projector
`F0(s)-C` was independently checked.

For both caps every piece, and their total rank, is at most `23-k ≤ 4`, below
half the factor dimension. Their inherited inner profiles are therefore
singleton lists. Rank telescoping gives exactly the same multiplicity as
the direct dimension count `dim F0(s)-k`. Retaining the identity vertex does
not alter that histogram.

The new auxiliary corners `C -> F`, changed target-chain transitions, and
target endpoints are also emitted as actual matrix differences. The checker
verifies the required north-east pivot pattern for **every distinct changed
edge under both original fixed axis bases U and V**. It does not choose a
different basis for each edge. The component with no predecessor includes
a rank-20 or rank-19 target entry and is checked explicitly as well.

## A concrete obstruction to cap 18

All slot and target indices here are zero-based indices of the pinned data.
Original readout slot **23806**, whose frame has dimension 19, contributes
with coefficient 1 to both targets:

| Target triple | Unchanged preceding slot | Predecessor dimension |
| --- | ---: | ---: |
| `(0,3,4)` | 24415 | 18 |
| `(0,3,5)` | 25074 | 18 |

Their two 18-dimensional frames span a **19-dimensional** space. The saved
certificate selects 18 basis rows from slot 24415 and the first basis row
from slot 25074. A specified 19-by-19 integer minor has determinant **1**.
Both predecessor spaces lie in the original 19-plane of slot 23806, checked
exactly against its annihilator, so the union rank is exactly 19 and their
intersection has dimension 17.

Any replacement for slot 23806 that occurs after those unchanged predecessor
reads must contain both spaces. A frame of dimension at most 18 cannot do so.
Both coefficients are odd, so this obstruction applies to the integer and
binary supports. Fourteen individual slots give such local obstructions in
the saved scan.

This rules out cap 18 **while retaining those lower frames and the readout
order**. It does not rule out modifying the lower frames, changing the
schedule, or using a different construction.

## Reproduction

Only standard Python is required. Both scripts reject `-O`.

```sh
python3 research/exploration/deferred-cap-target-10/check_caps.py --source /path/to/swapnil-round7 --cap 20
python3 research/exploration/deferred-cap-target-10/check_caps.py --source /path/to/swapnil-round7 --cap 19
python3 research/exploration/deferred-cap-target-10/check_cap18_obstruction.py --source /path/to/swapnil-round7
```

The cap receipts record all component memberships, source-frame slots,
integer extension rows, exact containment counts, and both-basis profile
checks. `cap18-obstruction.json` includes the actual determinant-one minor,
its source rows and columns, and the containing plane's exact annihilator.
The main theorem and published certificates are unchanged.
