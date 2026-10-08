"""Finite checks supporting the source-support lemma; the general proof is in the note."""
from itertools import combinations, combinations_with_replacement


def binary_rank(rows):
    basis={}
    for x in rows:
        while x:
            p=x.bit_length()-1
            if p in basis:x^=basis[p]
            else:basis[p]=x;break
    return len(basis)


def support_rank(c):
    return binary_rank(sum(1<<i for i in s) for s in combinations(range(len(c)),3)
                       if sum(c[i] for i in s))


def main():
    count=0
    # Coordinate permutations preserve the claim. Exhaust all multisets of
    # these five coefficient values, so this covers every arrangement as well.
    for h in range(8,17):
        for c in combinations_with_replacement(range(-2,3),h):
            if not any(c):continue
            assert support_rank(c)>=h-1,(h,c)
            count+=1
        for c in [(-2,)+(1,)*(h-1),(-2,-2)+(1,)*(h-2),
                  (1,)+(0,)*(h-1),(1,-1)+(0,)*(h-2)]:
            assert support_rank(c)==h-1
    assert (9-1)**2>84/2
    print('PASS',count,'coefficient multisets and all four hyperplane types at h=8..16')

if __name__=='__main__':main()
