# A temporal bound for the strict source-history model

This is a bound on a specified monotone common-frame model. It is not a bound
on every interchange algorithm, and in particular does not exclude paid
resets, fresh temporary erasure, or more general coupled address operations.
The source-borrowing construction counts data rows among its middle working
roles and therefore is consistent with the bound.

Let N(T) be the set of triples S with |S intersect T|=1.

**Two-neighborhood non-cover lemma.** For h>=5 and triples U,V distinct from T,
N(T) is not contained in N(U) union N(V).

For h5,6,7, `neighborhood_cover.py` exhausts every other target pair with T fixed.
Permutation symmetry covers every T. The respective neighbor counts are3,9,18,
and the maximum numbers covered by two other targets are2,8,15. The program
also checks h8..15, as independent small checks of the argument below.

For h>=8 put R=[h] minus T, so |R|>=5. Write A=U intersect R and B=V intersect R.
Both are nonempty, since U,V differ from T. For c in T, the triples in N(T) that
contain c correspond to every pair from R. The pairs whose triple also belongs
to N(U) are:

- pairs avoiding A, if c belongs to U;
- pairs crossing the cut (A,R minus A), if c does not belong to U.

The analogous statement holds for V. If c belongs to both U and V, choose a
pair meeting both nonempty A and B (using a second point if their selected
points coincide). That pair is in neither allowed family. If c belongs to
neither, the two allowed families are cuts. Assign to each point its two cut
membership bits. Since |R|>=5, two points receive the same bits, and their pair
crosses neither cut. Both cases contradict the proposed cover.

Therefore U intersect T and V intersect T partition T. Neither part can have
size3, as that would give U=T or V=T. One part has size1. Suppose it belongs to
U, with its point c. Then |A|=2. The triangle formed by the two points of A and
any third point of R has no edge avoiding A, so all three edges would have to
cross B. A cut contains no triangle. This is the final contradiction.

**Consequence for monotone source histories.** Consider binary common-frame
CNOTs on roles whose initial source history is empty or one input triple.
At each CNOT, both roles' forced histories acquire the union of their previous
histories, since they must share a frame containing both old frames. Histories
never shrink. Input injections occur at their source lines. A full side output
for T is read in t_T^perp. Because the input streams are independent, its forced
history must contain N(T); because its frame is orthogonal to t_T, that history
cannot contain any triple outside N(T). Thus its history is exactly N(T).

Call a history T-exclusive when it is contained in N(T) and not contained in
N(U) for any other target U. Singleton input histories are not exclusive.
Consider the first gate that creates a T-exclusive history along the relevant
computation. Its two prior histories fit some N(U),N(V), with U,V different
from T. By the non-cover lemma, their union omits an element of N(T). The gate
makes both incident roles T-exclusive, but has not completed the full history.
If only those two T-exclusive roles ever existed, subsequent gates between
them could not add the missing source. A gate bringing in that source through
another role makes that third role T-exclusive too. A fresh line injection
cannot add a distinct source directly after exclusivity, since the recipient
frame already contains at least two distinct source lines.

Hence at least three working roles become T-exclusive before its full output
can be produced. Once a role is T-exclusive, monotonicity prevents it from
being compatible with a different target's output. These groups are disjoint
across T, proving that at least3v working roles are needed in this model.

This bound counts all middle working rows, not just auxiliary rows. Borrowing
v existing X data rows permits an auxiliary count down to2v without violating
it. A producer operating only on v data rows cannot satisfy this strict model,
even though its scalar adjacency matrix is invertible. The argument uses time
order and source histories; it does not assert that every role in an
undirected interaction component inherits the whole component's inputs.
