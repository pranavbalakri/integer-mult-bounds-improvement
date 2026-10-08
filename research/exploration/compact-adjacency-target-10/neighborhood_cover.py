"""Exact small checks for the two-neighborhood non-cover lemma.
Permutation symmetry fixes the first target; all other target pairs are checked.
"""

if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')

from itertools import combinations
import json
from pathlib import Path


def check(h):
    triples=list(combinations(range(h),3));target=frozenset(triples[0])
    neighbors=[s for s in triples if len(set(s)&target)==1]
    masks=[sum(1<<i for i,s in enumerate(neighbors) if len(set(s)&set(t))==1)
           for t in triples[1:]]
    maximum=max((a|b).bit_count() for a,b in combinations(masks,2))
    assert maximum<len(neighbors)
    return dict(h=h,target=triples[0],neighbor_count=len(neighbors),
                maximum_covered_by_two_other_targets=maximum,
                all_other_target_pairs_checked=len(masks)*(len(masks)-1)//2)

if __name__=='__main__':
    rows=[check(h) for h in range(5,16)]
    print(json.dumps(rows,indent=2))
    Path(__file__).with_name('neighborhood-cover.json').write_text(json.dumps(rows,indent=2)+'\n')
