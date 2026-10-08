"""GF(2) fitting-matrix search. Solutions are scalar candidates, not framed networks."""
import sys,json,time
from itertools import combinations
from fractions import Fraction as Q
import z3


def rank_q(mat):
 a=[[Q(x) for x in row] for row in mat];k=0
 for j in range(len(a[0]) if a else 0):
  pivot=next((i for i in range(k,len(a)) if a[i][j]),None)
  if pivot is None:continue
  a[k],a[pivot]=a[pivot],a[k];vv=a[k][j];a[k]=[x/vv for x in a[k]]
  for i in range(k+1,len(a)):
   vv=a[i][j]
   if vv:a[i]=[x-vv*y for x,y in zip(a[i],a[k])]
  k+=1
 return k


def search(h,r,timeout=60000,symmetric=False):
 T=list(combinations(range(h),3));v=len(T)
 s=z3.Solver();s.set(timeout=timeout)
 U=[[z3.Bool('u_%s_%s'%(i,j)) for j in range(r)] for i in range(v)]
 V=U if symmetric else [[z3.Bool('v_%s_%s'%(i,j)) for j in range(r)] for i in range(v)]
 def dot(a,b):
  terms=[z3.And(x,y) for x,y in zip(a,b)]
  out=terms[0]
  for t in terms[1:]:out=z3.Xor(out,t)
  return out
 for i,a in enumerate(T):
  for j,b in enumerate(T):
   if i==j:s.add(dot(U[i],V[j]))
   elif len(set(a)&set(b))!=1:s.add(z3.Not(dot(U[i],V[j])))
 # Normalize a biorthogonal independent family; complete its dual bases when r is larger.
 I=[i for i,a in enumerate(T) if 0 in a and 1 in a]
 if len(I)<=r:
  for j,i in enumerate(I):
   for k in range(r):s.add(U[i][k]==(j==k),V[i][k]==(j==k))
 start=time.time();res=s.check();print('h',h,'rank',r,'symmetric',symmetric,'result',res,'seconds',time.time()-start,flush=True)
 if res!=z3.sat:return
 model=s.model();uu=[[z3.is_true(model.eval(x)) for x in row] for row in U];vv=[[z3.is_true(model.eval(x)) for x in row] for row in V]
 costs=[];rect=[]
 for k in range(r):
  targets=[i for i in range(v) if uu[i][k]];sources=[j for j in range(v) if vv[j][k]]
  G=[[len(set(T[i])&set(T[j]))-1 for j in sources] for i in targets]
  cost=rank_q(G);costs.append(cost);rect.append(dict(sources=[T[i] for i in sources],targets=[T[i] for i in targets],cross_gram_rank=cost))
 print('CROSS RANK COST',sum(costs),costs,flush=True)
 with open('fitting_%s_%s_%s.json'%(h,r,int(symmetric)),'w') as f:json.dump(dict(h=h,rank=r,rectangles=rect),f,indent=2)
if __name__=='__main__':
 h,r=map(int,sys.argv[1:3]);search(h,r,symmetric='symmetric' in sys.argv)
