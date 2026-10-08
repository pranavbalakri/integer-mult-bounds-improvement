# How far can source borrowing go?

These are architecture-specific necessary bounds and exploratory budgets, not a
new multiplication certificate. The certified files are unchanged. None of the
relaxations below is asserted to be a realizable circuit, and none reaches the
requested multiplication saving `kappa > 2^-10` with a checked construction.
The corresponding necessary bit saving is `a_bit > 1/1023`.

The main result is that reducing the producer to two auxiliary roles per target
is insufficient while retaining the current data endpoints and centre copies.
A proof valid for every pair of arities gives:

| Optimistic borrowed architecture | Necessary bit-saving ceiling |
| --- | ---: |
| Current two-block auxiliary selected edges | `a_bit < 1/2500` |
| Even ideal one-block auxiliary selected edges | `a_bit < 1/1700` |

A stronger target-specific check also shows failure whenever **both** axes
use at least `v_h/5` auxiliary roles with flag edges, or at least `2v_h/5`
with perfect selected edges. These are all-arity statements as well.

Both ceilings are below `1/1023`. The hypothesis is **at least two additional
roles per target on each axis**, with the other explicitly stated profile
assumptions. The earlier strict source-history model with separate full-output
roles motivates this role count; it is not a lower bound for arbitrary borrowed,
resetting, encoded, or coupled-data circuits.

## The exact relaxed profile

Write `v_h = binomial(h,3)`, `N=v_a v_b`, `m=ab`, and let `rho` be the auxiliary
role count per target on each axis. The number of arbitrary-input roles is
`W=(2+2rho)N`. Brackets below denote a child rank, and coefficients are per `N`.
The optimistic data histogram is

```
2[(a-1)(b-1)] + 2[a-1] + 2[b-1] + [1].
```

Its total rank is `2m-1`. The last singleton is the inherited paid endpoint
copy. Each data path has been granted complete batching within its three
inherited tensor increments. In particular, the indexed `23 x 25` data profile
is a refinement of this one: its large rank `528` is granted as the merger
`481+21+17+9*1`, its rank `22` as `21+1`, and its rank `24` as `23+1`, with two
copies of each path and one remaining endpoint singleton.

For each axis `h` in `{a,b}`, the current flag-shaped auxiliary profile is

```
rho * ([m-2h] + [h] + [h]).
```

The first two terms are its selected outer edge. The last grants all local
monotone growth a single rank-`h` child. The even more favorable selected-edge
relaxation replaces this by `rho * ([m-h] + [h])`.
The retained centre copies additionally contribute

```
(h / v_h) * [h-1]
```

on each axis. Consequently the available rank deficit is exactly

```
D/N = 1 - 6/(a-2) - 6/(b-2),
```

independent of `rho`. Actual finer flags cannot improve these relaxed moments:
for `0 < alpha < 1`, merging positive ranks lowers
`sum r^(1-alpha)`, and also lowers `sum r log(m/r)`. Source roles that have been
borrowed are included in the fully batched data path, not counted again as
auxiliaries. Increasing either axis's auxiliary ratio above two adds positive
moment excess while leaving `D` unchanged.

## Proof of the all-arity ceilings

Define

```
M(alpha) = sum_r n_r r (m/r)^alpha / (W m),
E        = sum_r n_r r log(m/r).
```

Since `exp(t) >= 1+t`, `M(alpha)<1` requires `alpha < D/E`. If `D<=0`, no positive
saving passes. Otherwise `a,b>=9`; by symmetry assume `a<=b`.
For `d=(a-1)(b-1)`, the relaxed data entropy per `N` is at least

```
2 d log(m/d)
+ (2a-1) log b + (2b-1) log a
+ 4 - 2/a - 2/b.
```

This follows by applying `log(h/(h-1)) >= 1/h` to the two smaller data ranks.
For `rho=2`, add, for each `h` and the other axis `o`:

```
flag:     2(m-2h) log(m/(m-2h)) + 4h log o;
perfect:  2(m-h)  log(m/(m-h))  + 2h log o.
```

The positive centre-copy entropy may be omitted for a weaker bound.
`check_ceiling.py` evaluates these lower bounds using rational arithmetic only:
24 positive atanh-series terms give downward log bounds, and eight positive
terms of `-log(1-t)` bound the large-block terms.

For the flag case, every `9<=a<=b<209` with `D>0` has `E>2500D`. In the remaining
tail, the displayed log terms include `(6b-1) log a`; since `log a>2` and
`b>=209`, they already exceed `2500`, while `D/N<1`.
For perfect selected edges, the finite check is `9<=a<=b<213` and `E>1700D`.
Its tail uses `(4b-1) log a > 1700` for `b>=213`. These cover all arities, not
merely a finite numerical search.

For the stronger role-count statements, the same exact lower bound is used
with `rho=1/5` for flag edges and `rho=2/5` for perfect edges. In each case every
`9<=a<=b<214` with positive deficit satisfies `E>1023D`. Both tails have log
coefficient `(12b/5-1) log a > 1023` for `b>=214`. The smallest finite slack
is positive at `24 x 24`. This proves target failure if both axes are above
the stated respective thresholds; it does not require each axis individually
to be below the threshold in an asymmetric possible construction.

## Exact target budgets

`requirements.py` brackets the unique positive roots using the independent
rational log/exponential enclosure. The following rounded numbers summarize
those exact intervals:

| Relaxation | Bit moment root |
| --- | ---: |
| Fully batched data, two auxiliary roles, flag edges, `24 x 25` | `0.000382335198…` |
| Same, even perfect selected edges | `0.000572500654…` |
| Fixed indexed `23 x 25` data, zero auxiliary roles, retained centres | `0.000523469425…` |
| Fixed indexed data, two auxiliary roles, retained centres | `0.000273221711…` |
| Fully batched `23 x 25` data, two auxiliary roles, no centre cost | `0.000842373283…` |
| Fully batched `24 x 25` data, zero auxiliary roles, retained centres | `0.001150728778…` |

The first two entries maximize a floating exploratory scan through arity 100;
their displayed individual roots are checked exactly. Global maxima are not
claimed. The last row is important: the zero-role, fully batched relaxation has
room for the target, so the preceding ceilings must not be presented as a
universal barrier caused by the triple geometry alone.

At `23 x 25`, fully batched data with retained centres requires
`rho < 0.1740685` to reach the target; the exact boundary lies between
`0.1740684` and `0.1740685`. Even with centre costs entirely removed, the fixed
indexed data requires `rho < 0.3976827` (boundary above `0.3976826`).
Deleting the paid endpoint singleton *hypothetically* would allow a relaxed
`13 x 13` profile with `rho=2` to have root `0.001686166214…`. No construction
that removes that singleton is supplied.

For a paid-reset route, the receipt also removes all old centre costs and then
counts the maximum total number of singleton charges per invocation that the
ideal `rho=2` profile could absorb at `1/1023`. For example the budgets are only
37 at `h=10`, 57 at `h=12`, and 85 at `h=16`; they must pay **all** new resets and
replacement centre transfers. The corresponding old retained-centre ranks are
90, 132, and 240. These are generous budgets, not a verified way to spend them.

## Zero-initialized auxiliaries that can be discarded

The dirty-role ceiling does not apply to private zero registers. Accordingly,
`zero_discard.py` gives them every advantage: a register starts directly at its
first used frame and disappears immediately after its last use. There are no
initial auxiliary edges, final full-frame edges, dirty cancellation, or inverse
cleanup. It replays the complete clean-input scalar map, counts only intervening
frame growth, and retains the necessary final side-frame extension and copied
centre transfer. Dead internal outputs incur no extra exit cost.

Let `q_h` be that surviving local growth rank per invocation, excluding centre
copies. Even with the ideal data profile, the resulting rank deficit is

```
D/N = 1 - 6/(a-2) - 6/(b-2) - q_a/v_a - q_b/v_b.
```

All tested existing producers fail already at this zeroth-moment rank budget.
The ordinary `h=23` producer has `q=329125`, or `q/v=185.8413…`, and
`D/N=-657491/1771`. The paired/reassociated producer has `q=332441`.
The pinned indexed PR84 producer has `q_23=355959`, `q_25=511939`, and its
rectangular profile has `D/N=-74934903/177100`.
A new symmetric zero-register producer retaining the same endpoints and
centres would need `q_h/v_h < (1-12/(h-2))/2` even to have positive savings;
at `h=23`, this is `<3/14`, over 867 times smaller than the ordinary producer.
Because the optimistic rank budget already fails, no unproved global transport
rule for discarded auxiliaries is needed to reject these particular attempts.
The exact scalar and histogram checks do not establish such a rule for future
positive candidates.

## In-place point-star transvections with paid resets

When `h=2 mod 4`, write `u_i` for the point-star incidence vector. Then
`u_i^T u_j=0` over `F_2` for every pair, including equal indices. Thus the
commuting transvections `I+u_i u_i^T` have product `I+B^T B`, precisely the
intersection-one side matrix, and preserve every original point total.
`point_stars.py` checks this complete scalar identity and involutivity for
`h=6,10,14,18,22,26`.

Grant a rank-two switch from the point-star frame `E_i` to `E_j` through the
full frame. Every data role `X_T` still visits three distinct star frames, so
it has two interior switches whatever the point ordering. The last star frame
never equals any target frame `K_U=t_U^perp`. Indeed, their Euclidean normal
vectors are respectively `1-3e_i` and `3t_U-1` for the Gram `I-J/9`. Choose one
index in `U` and one outside `U`, both different from `i`; the first normal has
the same coordinate at both, and the second does not. Hence these normals are
not proportional for `h>=5`. Permuting target roles cannot avoid the final
switch either.

The resulting source path pays

```
(h-2) + 2*2 + 2 + 1 = h+5
```

ranks: initial line to first star, two internal switches, final target switch,
and target to full frame for restoration. The inherited direct source path
pays `h-1`. This literal construction therefore adds `6N` per stage, or `12N`
in total, while the original rank deficit is less than `N`. It fails even if
all its scratch registers are free. This rules out that particular star-visit
implementation, not other factorizations or altered data endpoints.

## Reproduce

From the repository root, using standard Python with assertions enabled:

```sh
python3 research/exploration/borrowed-ceiling-target-10/check_ceiling.py
python3 research/exploration/borrowed-ceiling-target-10/requirements.py
python3 research/exploration/borrowed-ceiling-target-10/point_stars.py
python3 research/exploration/borrowed-ceiling-target-10/zero_discard.py --source /path/to/pinned/PR84
```

The final command verifies the source revision before reading its literal
words. The external producer is credited and pinned in
[the upstream notice](../indexed-cycle-target-10/UPSTREAM-NOTICE).
Omitting `--source` runs only the repository's ordinary and paired/reassociated
producers. All verification entry points reject `python -O`.
