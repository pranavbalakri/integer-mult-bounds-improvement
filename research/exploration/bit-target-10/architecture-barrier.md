# A barrier for improving the retained-total two-stage architecture

This is a necessary-condition bound for the existing architecture, not a bound
on all possible integer-multiplication algorithms.

Let `v=C(h,3), m=h², N=v²`. Retaining the same dirty-scratch copied-total schedule
gives rank deficit

```
D = Wm-s = N - 2vh(h-1) = N (h-14)/(h-2).
```

Consequently a positive saving requires `h>=15`. The existing scalar-role
compiler has `R>=3v+h`: its outputs alone require `3v` side outputs and `h`
retained uses, before paying for a single addition. Hence `R>=3v`.

A necessary condition for the normalized recursion moment to be below one at
saving `a>0` follows from `exp(z)>=1+z`:

```
a < D / E,
E = sum_children width * log(m/width).
```

Consider only the `2vR` auxiliary roles. Each has a selected edge of rank `m-h`
and a remaining chain of total rank `h`. Grant an unrealistically favorable
batching in which the selected edge is one child of width `m-h`, and the entire
remaining chain is one child of width `h`. Splitting a width only increases its
contribution to `E`, so this is a valid optimistic lower bound:

```
E >= 2vR [(m-h) log(m/(m-h)) + h log h]
  >= 6N [(h-1) + h log h].
```

Here `-log(1-1/h)>=1/h`. For `h>=15`, `log h>8/3` and `1-1/h>=14/15`, giving

```
E > (108/5) N h,
a < 5(h-14)/(108 h(h-2)) < 1/1080.
```

The last inequality is exactly

```
h(h-2)-50(h-14) = (h-26)^2 + 24 > 0.
```

For completeness, the elementary logarithm bound can itself be made rational:
`log 15 = 3 log 2 + log(15/8)` and positive atanh-series terms give
`log 15 > 3*2*(1/3+1/81)+2*(7/23) = 1666/621 > 8/3`.

Therefore **even zero-cost additions and perfect batching cannot reach
`a_b>1/1023` while preserving these centre losses and role structure**.
The target `kappa>2^-10` needs `a_b>1/1023` under the retained assembly's
reservation constraint. Improving the leave-one-out circuit or its ordering
cannot suffice. One must change the central-return schedule, the number/type
of output roles, the motif, or the assembly itself.

The `barriers.py` exploration also keeps the actual data-edge profiles and the
current two-block selected auxiliary profile. Its optimistic finite scan
`15<=h<=200` peaks near `0.0002114` even with no addition roles, versus the target
`0.0009775`. Those floating-point figures illustrate the gap; the all-dimension
rational argument above is the rigorous obstruction.
