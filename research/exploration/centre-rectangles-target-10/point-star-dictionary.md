# Exact h=14 obstruction for the full point-star/exclusion rectangle dictionary

This proves a restriction on one explicit dictionary, not on arbitrary centre
functions or phase circuits. A positive construction must leave this dictionary.
Coefficients may be arbitrary rational or complex numbers, and the number of
centres need not be rank-minimal.

Let P_i(T)=1 when i belongs to T and E_i(T)=1-P_i(T), for triples T.
Allowed centres are nonzero scalar multiples of any outer product of these
2h source functions and 2h target functions. Use the cheapest admissible
common read frame from the inherited nonalternating whole-residual interface.

For even h, every source support label has dimension h-1. A rectangle has
return rank at least h-2. The only rectangles attaining h-2 are
P_i P_j and E_i E_j with i != j. The same-index mixed rectangle has an
alternating residual and therefore cannot attain h-2 in this interface.
All other rectangles have rank at least h-1.

## The fitting constraints fix the scalar matrix inside this code

Let B be the h-by-v incidence matrix. Any sum of dictionary rectangles is
B^T H B for a matrix H. Suppose its diagonal is one and its entries vanish
when distinct triples meet in one point. For h>=6, fix a triple T and write
z=H t_T. For two points a,b outside T, a third point c outside T exists.
Comparing target triples {i,a,c} and {i,b,c}, for i in T, gives z_a=z_b.
Call this common outside value q. Every target triple meeting T in one point
then gives z_i+2q=0 for i in T. The diagonal gives -6q=1. Hence

    H t_T = (1/2)t_T - (1/6)1

for every triple T. Triple indicators span the full rational coordinate
space, and therefore necessarily

    H = (I-J/9)/2.

For h=14 this has rank 14, so at least 14 nonzero rank-one centres are needed.

## Fourteen centres cannot be a dictionary factorization

With exactly h=14 centres the chosen source rows form a basis. A basis
selected from the P_i and E_i has one of two forms:

1. One choice P_i or E_i for every index i.
2. Both choices at one index, no choice at another, and one choice at all
   remaining indices.

There cannot be fewer than h-1 represented indices, since each additional
index adds at most one dimension beyond the all-ones row.

For case 1, write k for the number of exclusions, I for their index set,
and s_i=-1 on I and +1 elsewhere. The source matrix is singular exactly
when k=3. Otherwise its uniquely forced target rows are

    g_i = (s_i/2) [ e_i + 1_I/(3-k) - 1/(3(3-k)) ].

For h=14 these cannot all be scalar multiples of P_j or E_j. This is
verified with exact rational arithmetic for k=0,...,14; permutations make
these all the cases. In case 2, the uniquely forced target rows at each
ordinary represented index i are scalar multiples of e_i-e_k, where k is
the omitted index. These are outside the dictionary. The exact checker
also checks every exclusion count among the other 12 indices.

Thus a 14-centre factorization is impossible, regardless of read cost.

## Fifteen centres cannot have total cost below 182

Fifteen cheapest centres already cost 15*12=180. A total below 182 permits
at most one centre of cost at least 13; the other 14 must be cheap
P_i P_j or E_i E_j rectangles with i != j.

For every unordered coordinate pair define the linear invariant

    J_ij(M) = M_ij + M_ji - M_ii - M_jj.

For the required H, all 91 invariants equal -1. For a rank-one term fg^T,

    J_ij(fg^T) = -(f_i-f_j)(g_i-g_j).

A cheap rectangle changes the invariant of exactly one unordered pair.
Any dictionary rectangle changes at most h-1=13 pair invariants (the
maximum is attained when both function indices coincide). Consequently
14 cheap rectangles plus one other can cover at most 14+13=27 of the
91 required nonzero pair invariants. This is impossible even with arbitrary
cancellations or complex coefficients.

Sixteen or more centres cost at least 16*12=192. Combining the cases,
every fitting decomposition in the full dictionary costs at least 182.
Since v=binom(14,3)=364, this cannot make the copied two-stage deficit
v(v-2*centre_cost) positive.

The simpler variant allowing only the cheap rectangles requires at least
one centre per unordered pair, costing at least
(h-2)binom(h,2)=3v. It is especially far from the target.

Run `python3 check_rectangles.py --dictionary`; the exact receipt is
`dictionary-results.json`. These checks verify the finite h=14 basis
classification and arithmetic. They are not a general lower bound beyond
the stated dictionary. Written with OpenAI Codex assistance.
