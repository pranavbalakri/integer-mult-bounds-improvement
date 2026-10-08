# A rectangular-arity obstruction for the 2^-10 search

This is a necessary-condition result for the family specified below. It is not
an impossibility result for integer multiplication, and is not a new exponent
saving. It does not modify the previously certified construction.

The balanced-arity obstructions do not by themselves exclude rectangular
choices. This certificate closes that particular gap, under explicit role and
centre assumptions. It even grants a better child partition than our compiler
can generally implement.

Assume two tensor factors of integer dimensions a,b, with v_h=C(h,3),
N=v_a*v_b, m=a*b. In each dimension h assume at least 3*v_h auxiliary roles per
invocation, one exterior edge of rank m-h per role, a remaining chain of total
rank h, copied-centre loss h*(h-1) per invocation, and one endpoint rank per
data pair. These are assumptions to check for any new compiler; a future
construction need not satisfy them.

The normalized rank deficit is

    D/N = d = 1 - 6/(a-2) - 6/(b-2).

Nonpositive d gives no positive moment saving. Positive d forces a,b>=9.
For a candidate saving s, exp(z)>=1+z makes s < D/E necessary, where

    E = sum(child_width * log(m/child_width)).

Merging children reduces E. Give every auxiliary role just two ideal children,
m-h and h, and ignore all data and centre contributions to E. The role lower
bounds and -log(1-x)>=x then give

    E/N >= 3*[a*(log(b)+1-1/b) + b*(log(a)+1-1/a)].

The accompanying exact rational program proves that this lower bound exceeds
1300*d for every a,b>=9 with d>0. By symmetry take a<=b. It checks the finite
range b<=144 using positive atanh-series lower bounds for logarithms. For the
infinite tail, log(9)>3*(2/3)+2/17>19/9, so the second summand alone gives
E/N>9*b>=1305 when b>=145, while 1300*d<1300.

Consequently every member satisfying the assumptions has

    s < 1/1300 < 2^-10.

The retained assembly is slightly stricter: its reservation and butterfly
constraints require s>1/1023 to obtain kappa>1/1024. Thus varying two unequal
arities, reordering a producer, or merely improving batching cannot reach the
target while preserving the stated lower bounds and losses.

Run `python3 research/exploration/rectangular-target-10/check_barrier.py`.
The certificate covers all dimensions; it does not infer an infinite claim
from a finite numerical parameter sweep. Escaping it requires changing at
least one stated hypothesis, such as the centre loss, role budget, endpoints,
or transfer/assembly contract.
