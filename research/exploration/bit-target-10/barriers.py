"""Optimistic relaxations of the retained-total two-stage bit architecture.
These are not realized constructions or certified lower-bound witnesses.
"""
from math import comb,log


def histogram(h, centres=True, whole_aux=False, ideal_data=False):
 v=comb(h,3);m=h*h;N=v*v;R=3*v+h;A=v*R;H={}
 def add(w,n):
  if w>0 and n:H[w]=H.get(w,0)+n
 for _ in range(2):
  if whole_aux:add(m-h,A)
  else:add(m-2*h,A);add(h,A)
  add(h,A) # optimistic: all remaining side-chain rank in one width-h child
  if centres:add(h-1,v*h) # optimistic: copy loss in one child
 if ideal_data:
  add(m-1,2*N);add(1,N)
 else:
  add(m-4*h+2,2*N);add(h-2,2*N);add(1,(h+1)*2*N)
  for _ in range(2):add(h-2,2*N);add(1,2*N)
  add(1,N)
 W=2*N+2*A;s=sum(w*n for w,n in H.items());D=W*m-s
 def moment(a):return sum(n*(w/m)**(1-a) for w,n in H.items())/W
 lo,hi=0.,1.
 if D<=0:return dict(h=h,root=0.,derivative_bound=0.)
 for _ in range(60):
  mid=(lo+hi)/2
  if moment(mid)<1:lo=mid
  else:hi=mid
 entropy=sum(w*n*log(m/w) for w,n in H.items())
 return dict(h=h,root=lo,derivative_bound=D/entropy,R=R,W=W,D=D)
if __name__=='__main__':
 for centres in (True,False):
  for aux in (False,True):
   for data in (False,True):
    results=[histogram(h,centres,aux,data) for h in range(5,201)]
    best=max(results,key=lambda r:r['root'])
    print('centres',centres,'whole_aux',aux,'ideal_data',data,best,flush=True)
