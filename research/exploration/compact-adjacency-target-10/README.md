# Compact full-adjacency exploration

This directory now contains a constructive borrowed-source candidate, earlier
small producer experiments, and an exact bound on one proposed change of set
family. It does not modify the certified construction or claim the requested
2^-10 multiplication witness.

## Borrowing the source data roles

`borrowed_sources.py` replaces each leaf-pivot auxiliary by its existing X data
row. A fixed early correction cancels the arbitrary auxiliary contribution;
the middle producer, centre-first scatter, side scatter and full-frame cleanup
then realize the same shear and restore every arbitrary public auxiliary.
`borrowed-sources.tex` gives the scalar identity, frame chronology, intended
same-order inverse-transpose second stage and conditional profile derivation.

The exact finite replay passes h5..8 and the complete optimized h23 producer,
including every public input basis vector (41848 at h23), source-span inclusion,
and both shear orientations. Meaningful controls reject a missing dirty-value
correction and side scattering before the centre reads. The sparse early
correction uses fresh zero temporary source pivots; its transpose explicitly
clears their private output contents. Those temporaries are not arbitrary
public input roles counted in the volume denominator.

For our h23 producer, the reconstructed histogram removes 1,771 auxiliary roles per
invocation, keeps rank deficit 1,344,189, and certifies candidate bit saving
38238149/10^12 versus 36943733/10^12 (about 3.5%). See
`borrowed-profile-h23.json` and `borrowed_profile.py`. This is a finite-network
candidate with an exact moment; the full global transfer and multiplication
interfaces remain separate inherited proof obligations.

## Borrowed wrapper on the credited PR84 rectangle

`borrowed_indexed.py` applies the wrapper to the pinned external h23/h25
producer, with full arbitrary-input replay in both intended orientations.
`check_borrowed_rect.py` independently reconstructs the new histogram from the
pinned local profiles. Its result agrees exactly with subtracting, for each
axis, N copies of `{m-2h,h,1,1,h-2}`. The checked finite counts are
W=124,319,508, total rank 71,481,870,200, and unchanged deficit 1,846,900. The bit moment
certifies 5.488e-5. The moment at 1/1023 is strictly above one, so this candidate
does not attain the requested multiplication target 2^-10.

`borrowed_assembly.py` checks two separate parameter choices under the same
inherited global interfaces: conditional kappa=5.487e-5 with epsilon=0.9999, or
kappa=3.89e-5 with epsilon=0.71 and precision exponent 2.15. These are different
tradeoffs, not simultaneously available parameter settings. No end-to-end
runtime or crossover is established. Attribution and the pinned external
source are retained in `../indexed-cycle-target-10/UPSTREAM-NOTICE`.

The geometry explicitly places every surviving auxiliary at D0, and inserts
each X source-line entrance before its first middle gate. This preserves the
old selected and data-corner idempotents even when a source has no actual gate
at its initial line. The finite moment and the new local frame argument do not
independently prove the inherited global machine, routing, prime, recursion,
or analytic interfaces.

## Full intersection-one outputs

For triples, write `S[T,U]=1` when `|T intersect U|=1`. If `B` is point incidence,
then over F2, `S=I+B^T B`. When `h=2 mod4`, `BB^T=0`, hence `S^2=I`.
This scalar identity does not provide a free address-frame implementation.

`search_positive.py` constructs shared positive addition DAGs for the entire
output `sum_{U:|U intersect T|=1} x_U`, rather than three separate common-point
pieces. The two saved circuits have these independently checked counts:

| h | v | Additions | Full outputs | Side roles | W with direct copied centres | Charged rank sum | Wm-s |
|---|---:|---:|---:|---:|---:|---:|---:|
| 5 | 10 | 20 | 10 | 30 | 900 | 22900 | -400 |
| 6 | 20 | 93 | 20 | 113 | 5560 | 201200 | -1040 |

`check_positive.py` recomputes every support and its rational span rank, emits and
checks the source-injection/CNOT coefficient propagation, checks each physical
role's monotone source-frame chain and output orthogonality, and reconstructs
the complete inherited child histogram. Here `G=I-J/9` is positive definite,
so every source span and every nested orthogonal difference is nondegenerate.
The direct copied centres are charged: their loss is `2vh^2`. Neither example
has a positive rank deficit, so neither can certify a positive saving. The
checker proves the listed finite producer properties; the global transfer
interface remains inherited.

Simply appending two additions to merge the old three pieces removes two
terminal uses but adds two DAG nodes. Under the old `R=additions+outputs`
compiler this changes no role count. The affected small frame increments also
split into the same singleton children at the dimensions of interest. A real
improvement needs different synthesis or role reuse, not just merged names.

### A restricted temporal role bound

Consider a producer with input injections at their source lines, binary
common-frame CNOTs, monotonically growing labels, and full adjacency outputs
whose source spans equal the distinct hyperplanes `t_T^perp` (in particular,
ordinary triples at h6). At the last CNOT producing a full output, both roles
must reach that hyperplane. Its donor cannot have provided, or later provide,
a different full output: monotonicity would force its label to contain two
distinct hyperplanes. The same donor cannot serve a different full target's
last CNOT. Consequently there is at least one distinct non-output donor for
each output, and `R>=2v` in this restricted model.

This is a temporal argument, not a claim about undirected connected components:
a wire can leave a computation before later mixing, and its history need not
contain every input in that component. The bound does not cover paid descents,
newly erased temporaries, or arbitrary coupled address transformations, and it
does not establish `R>=3v` in general.

A stronger bound for the strict source-history model is proved in
`neighborhood-bound.md`: no target neighborhood is covered by two other target
neighborhoods, forcing three working roles per complete output. This counts
borrowed data rows as working roles, so it allows an auxiliary floor of 2v when
v data rows participate. It is not a universal interchange lower bound.

## Five-element sets with intersection three forbidden

Using every five-element set with `G=I-J/25` and the old point-incidence centres
is invalid. Two distinct sets meeting in three points have scalar point-centre
coefficient one, but Gram inner product two; an intersection-one side circuit
cannot cancel that coefficient. Restricting the family to forbid intersection
three repairs that particular scalar issue, but reduces its size.

`k5_lp.py` computes exact Delsarte bounds for this restricted family. With
Johnson distance `i=5-|S intersect T|`, its normalized distance distribution
satisfies

```
A_0=1, A_2=0, A_i>=0,
sum_i A_i P_i(j)/v_i >=0       (j=1,...,5),
v_i=C(5,i) C(h-5,i),
P_i(j)=sum_{t=0}^i (-1)^(i-t) C(5-t,i-t) C(5-j,t) C(h-5+t-j,t).
```

The objective `sum_i A_i` equals the family size. These are the usual
positive-semidefinite Johnson-scheme constraints; the same constraints apply
when one distance is forbidden rather than all distances below a threshold.
See the association-scheme LP discussion in
[Chailloux and Debris-Alazard, *New Solutions to Delsarte's Dual Linear Programs*](https://arxiv.org/html/2405.07666v1#S2.SS6)
for the underlying LP principle.

There are four variables, at distances 1,3,4,5. The program enumerates rational
vertices, then separately checks a nonnegative rational dual combination giving
the same upper bound. It additionally verifies all eigenmatrix orthogonality
relations and the Johnson adjacency recurrence exactly. Thus the bounds do not
rely on a floating-point solver. Both primal and dual certificates for h10..40
are in `k5-lp.json`.

| h | Exact LP upper bound | Integer family upper bound | 2h(h-1) |
|---|---:|---:|---:|
| 10 | 56/3 | 18 | 180 |
| 15 | 231/2 | 115 | 420 |
| 20 | 1083/2 | 541 | 760 |
| 22 | 31878/37 | 861 | 924 |
| 23 | 168245/158 | 1064 | 1012 |
| 24 | 9108/7 | 1301 | 1104 |

The proposed size condition `|F|>2h(h-1)` is impossible through h22; h23 is the
first dimension this LP leaves open. An LP-feasible distribution is not a
family construction. Nor does this comparison rule out every smaller family:
a special family might reduce its actual sum of point-star ranks below
`h(h-1)`, which would require a separate full rank and producer audit. The form
`I-J/25` is singular at h25, another condition a later construction must check.

## Reproduction

Run from the repository root:

```
python3 research/exploration/compact-adjacency-target-10/check_positive.py
python3 research/exploration/compact-adjacency-target-10/k5_lp.py
python3 research/exploration/compact-adjacency-target-10/borrowed_sources.py
python3 research/exploration/compact-adjacency-target-10/borrowed_sources.py --optimized 23
python3 research/exploration/compact-adjacency-target-10/borrowed_profile.py
python3 research/exploration/compact-adjacency-target-10/check_borrowed_rect.py
```

All use only Python's standard library and the bundled producer/certificate modules. Repeating the exploratory producer
search is optional; the saved candidates are checked directly.
