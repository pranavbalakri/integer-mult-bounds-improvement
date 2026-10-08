# Research toward a saving above 2^-10

**No new exponent witness is claimed here.** These files record failed shortcuts,
a rigorously scoped architectural obstruction, and one exact scalar fitting
candidate that does not survive the retained rank budget. They do not change the
certified construction.

## What was established

- `architecture-barrier.md` proves `a_b<1/1080` for the stated balanced retained-total
  architecture, even granting zero-cost additions and perfect batching. The sibling
  `../rectangular-target-10/` gives the stronger `1/1300` result for unequal arities.
  Neither result applies to every possible motif or multiplication algorithm.
- `fitting-h7-rank6.json` contains a GF(2) matrix on the 35 triples of seven points.
  It has diagonal one, vanishes off the diagonal when intersection size is zero
  or two, and factors into six rectangles. Its scalar rank is six, compared with
  seven for the usual point-incidence Gram matrix.
- The six rectangles have exact rational cross-Gram ranks `4,6,7,5,7,6`, summing
  to 35. Even granting that sum as the complete copied-centre transfer loss, the
  retained two-stage deficit would be `35²-2*35*35=-1225`. Thus this scalar
  improvement does not yield a positive recursion saving. A valid frame compiler
  for the candidate has not been built.
- At `h=6`, the intersection-one adjacency matrix squares to the identity over
  GF(2). However, it mixes roles with different rank-one address projectors and
  fails to commute with those role-dependent frames. Its scalar inverse therefore
  cannot simply be applied as a free endpoint repair. The exact checker includes
  this negative control.

Run the standalone, standard-library checks:

```
python3 research/exploration/bit-target-10/check_fitting.py
python3 research/exploration/rectangular-target-10/check_barrier.py
```

## Exploratory numerical diagnostics

`barriers.py` evaluates deliberately optimistic, generally unrealized child
histograms. `hypotheses.py` uses the existing paired producer but artificially
removes one or both centre-return costs. These scripts are parameter diagnostics,
not algorithm certificates. The latter shows why small changes are insufficient:
halving the cost peaks around `9.04e-5`; removing it entirely would pass the target
at `h=6`, but no valid zero-cost schedule is supplied.

## Optional scalar fitting search

`fitting_search.py` requires the external `z3-solver` Python package; it is not
vendored here. It searches factorizations `B=UV^T` over GF(2), with diagonal one
and the required forbidden entries zero. An independent family of triples
containing a fixed pair is normalized to dual coordinate vectors. Such a
normalization is valid for scalar rank searches, but it changes rectangle costs,
so the resulting factorization must not be called a minimum-transfer-cost one.

Example, after installing that optional dependency in an environment of choice:

```
python3 research/exploration/bit-target-10/fitting_search.py 7 6
```

The script writes a candidate in the current directory. The saved candidate is
independently validated without Z3. During this exploration Z3 reported
unsatisfiability for scalar ranks at most four at `h=5`, at most five at `h=6`,
at most five at `h=7`, and at most seven at `h=8` (using the maximal tested rank
in each case; smaller ranks are contained in it). These search reports are not
used in the architectural proof or in any exponent certificate, and solver proof
certificates are not included.
