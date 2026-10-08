"""Unrealized hypothetical removal of a fraction of centre-return rank.
Numeric diagnostics only: these are NOT certified algorithmic improvements.
"""
import sys
from math import comb
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'independent' / 'bit-improvement'))
import paired_loo,reassociate,rtgm
from collections import Counter
from certificate_round5 import inner_children


def run(h,stages):
 G=paired_loo.build(h);C=rtgm.compile_(G);v=comb(h,3);N=v*v;m=h*h;R=C['roles'];A=v*R;H={};rk=Counter()
 for ss in range(R):
  d=rtgm.chain_dims(G,C,ss);rk.update(b-a for a,b in zip(d,d[1:]))
 def add(w,n):
  if w>0 and n:H[w]=H.get(w,0)+n
 for _ in range(2):
  add(m-2*h,A);add(h,A)
  for r,n in rk.items():
   for w in inner_children(r,h):add(w,v*n)
 for _ in range(stages):
  for w in inner_children(h-1,h):add(w,v*h)
 add(m-4*h+2,2*N);add(h-2,2*N);add(1,(h+1)*2*N)
 for _ in range(2):add(h-2,2*N);add(1,2*N)
 add(1,N);W=2*N+2*A;s=sum(w*n for w,n in H.items());D=W*m-s
 if D<=0:return None
 def moment(a):return sum(n*(w/m)**(1-a) for w,n in H.items())/W
 lo,hi=0.,.1
 for _ in range(70):
  md=(lo+hi)/2
  if moment(md)<1:lo=md
  else:hi=md
 return dict(h=h,R=R,root=lo,deficit=D)
if __name__=='__main__':
 for stages in (2,1,0):
  rows=[run(h,stages) for h in range(5,29) if h!=9]
  rows=[r for r in rows if r]
  print('CENTRE STAGES',stages,'BEST',max(rows,key=lambda r:r['root']),flush=True)
  print(rows[:6],flush=True)
