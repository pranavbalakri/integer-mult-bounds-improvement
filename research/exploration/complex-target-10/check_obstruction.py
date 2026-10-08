from fractions import Fraction as Q
from math import comb
K=1028
assert Q(4*22*(12*K+36))>Q(3*K+46,3)**2
assert Q(1093,405)>Q(8,3)
assert Q(1,1028)<Q(1,1024)
for h in range(3,15):
 v=comb(h,3)
 assert v*v-2*v*(h*h-h+1)<=0
for h in range(15,10001):
 v=comb(h,3);u=Q(2*(h*h-h+1),v)
 assert u>Q(12,h)
 q=22*h-Q(46,3)+Q(4,h)+Q(8,3)*u
 assert q>K*(1-u)
print('PASS exact arithmetic; general h proof uses AM-GM inequality in note')
