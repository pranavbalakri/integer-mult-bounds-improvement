# A paired-tree and reassociation improvement to the retained-total bit producer

At `h=23`, the exact exponential-moment certificate gives

```
a_b = 36943733 / 10^12 = 0.000036943733
R = 40077; W = 148225616; m = 529.
```

The original round-six bit witness is `36667/10^9`. The new witness is about
0.755% larger. This is a conditional improvement to the bit primitive. A complete
integer-multiplication claim additionally needs the complex primitive and the
inherited assembly, tape-layout, arithmetic, precision, and phase contracts.
It is not a measured runtime improvement.

## Construction

`paired_loo.py` changes only the leave-one-out subroutine of the paired exclusion
producer. Delete zero inputs. For a nontrivial list, rotate its last entry to the
front and group adjacent entries in pairs. Form each pair's total, recursively
compute all pair exclusions, and add back the surviving member of the omitted
pair. For each original zero input the omitted sum is simply the total. The base
case contains at most two nonzero entries.

Induction on the list length proves that output `i` is the sum of all inputs
except input `i`. Every internal addition has disjoint positive supports.
The original pair-exclusion recursion and global-matching order are unchanged.
The tree alone certifies `36881/10^9` with the original conservative moment bound.

`reassociate.py` then revisits an existing sum with support `S` and finds existing
nodes with disjoint supports `A,B` satisfying `A union B = S`. It replaces the
operands only when an integer score strictly improves. It creates no new support.
Strict inclusion of operand supports prevents cycles, and a final sort by support
cardinality produces a topological order. At `h=23` this changes 368 sums and
retains the same 34,741 addition nodes and 40,077 roles.

The integer score concentrates recursive child widths. Write `I(r)` for the
inherited inner profile: `r` singleton children when `2r<=h`, and `h-r` singleton
children plus one width-`2r-h` child otherwise. For a sum of dimension `d` with
operand dimensions `a,b`, the terms changed by reassociation are

```
I(a), I(d-a), I(b), I(d-b).
```

The two `I(a),I(b)` terms account for changed operand fan-out slots; the other
two account for input transitions at the sum. We maximize the sum of squared
widths in these four profiles. This is a deterministic search heuristic, not an
assertion of global moment optimality. The final emitted network's full histogram
and exact moment inequality certify the actual improvement.

## Why the inherited frames remain valid

Every active support has a common point `c`. For triple indicator vectors
`t_(cij)` and the inherited Gram form `G=I-J/9`,

```
<t_(cij), t_(ckl)>_G = |{i,j} intersection {k,l}|.
```

Consequently the source-span Gram matrix is `B^T B`, where `B` is the unsigned
incidence matrix of the corresponding pairs on the other `h-1` points. It is
positive definite on the source span: `B x=0` also forces the common coordinate
of the triple combination to vanish. Its dimension is the incidence rank,
`number of incident vertices - number of bipartite connected components`.
Every source span therefore remains nondegenerate. Disjoint unions give nested
source spans along each physical slot. Each target output still contains precisely
the source triples that meet its target triple in the single common point, so its
source span lies in the target-line complement. Each retained total is the full
point star and has dimension `h-1`. Reassociation preserves every one of these
supports; hence it changes neither these subspaces nor the admissible frame
family. The existing generic common flag-basis proposition applies to the finite
family of inherited nested frame differences.

The original selected auxiliary and data-edge profiles and copied-retained-total
schedule remain in force. The final rank identity is

```
s = W*m - N + 2*v*h*(h-1),   v=binom(h,3), N=v^2.
```

## Reproduce

From the `research` directory:

```
python3 independent/bit-improvement/check_paired_loo.py --reassociate 23
python3 independent/bit-improvement/certificate.py
```

The first command checks 2,047 leave-one-out cases with all zero patterns of lengths
0 through 10, every active global node's exact support and common point, all 5,313
side outputs and 23 retained totals, all 40,077 nested slot chains, every source
rank by the incidence formula, and all 1,970 local active-node ranks independently
by Gaussian elimination modulo an odd prime. The identity of ranks over the
rationals follows from the incidence-matrix argument above, rather than from the
finite-field checks alone.

The second command uses the existing `complex-twostage/cert.py` exact logarithm and
exponential upper bounds, and writes `certificate-h23.json`. Its rigorous moment
slack is positive; the next `10^-12` grid point does not pass that upper-bound
certificate. This does not claim mathematical optimality beyond that certificate.
The arithmetic checks do not independently verify the full inherited algorithm.
