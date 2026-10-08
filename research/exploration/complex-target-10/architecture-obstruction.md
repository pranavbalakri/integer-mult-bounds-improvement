# Complex two-stage architecture obstruction

Scope: the current triple-incidence two-stage complex transfer, copied retained centres, the current endpoint correction and data-frame paths, arbitrary side producer and arbitrary whole-residual batching of its small edges. This is not an upper bound on multiplication algorithms.

Let h>=3, m=h^2, v=C(h,3), N=v^2, W=2N+2vR, L=2v(h^2-h+1), and D=N-L. Every target has a nonzero side correction, so the compiled producer needs at least one side piece per target: R>=v. In fact it needs more, but that is not used.

The data histogram is forced:
- 2N children of width (h-1)^2;
- 4N children of width h-1;
- N endpoint-correction children of width 1.

Each of the 2vR auxiliary roles has one child of width m-h. The remaining auxiliary rank is 2vRh+L, entirely in children of width at most h. These facts follow directly from the two stage lifts and remain true if the internal side producer is changed.

Write Q=sum n_r r ln(m/r). The rank sum is s=Wm-D. A saving a requires

  M(a)=sum n_r r/(Wm) exp(a ln(m/r)) < 1.

Since exp(z)>=1+z, a<D/Q is necessary. If D<=0, no positive a can work. Directly D>0 requires h>=15. For h>=15, ln h>8/3. One elementary certificate is

  ln15 = 4ln2-ln(16/15)
       >4(2/3+2/81)-1/15 =1093/405 >8/3,

using the first two positive terms in the atanh series for ln2 and ln(1+t)<=t.

Also ln(h/(h-1))>1/h. The forced classes therefore give, with u=L/N,

  Q/N > 4(h-1)^2/h + (4h-2)(8/3)
        + 2(R/v)(h-1) + (2(R/v)h+u)(8/3)
      >=22h -46/3 +4/h +(8/3)u.

Here u=12(h^2-h+1)/[h(h-1)(h-2)]>12/h. For K=1028,

  (Q-KD)/N
    >22h -(K+46/3)+(12K+36)/h
    >=2sqrt(22*12372)-3130/3 >0.

The last inequality is exact because

  4*22*12372=1088736 > (3130/3)^2=9796900/9.

Thus every network in this architecture satisfies a_c<1/1028<1/1024.
No changes of side circuit, its role count, or batching alone can reach a_c>2^-10 while preserving the stated data/centre/endpoint architecture. Removing or reducing L, changing endpoint savings, or changing the motif avoids this obstruction and remains open.

# Why one scalar stage fails to give the desired missing factor h

Take m=h and N=v with one forward scalar shear y<-y+x. With input frames (P,I) and output frames (F,E=FP^-1), its data outputs are

  A=FP^-1 x,
  B=FP^-2 x+FP^-1 y.

The direct correction is B<-PB-A, followed by A<-PA. This uses two rank-one recursive phase children per data pair. The data edges saved two ranks per pair, so these children consume the entire saving before any retained-centre loss is charged.

Within the restricted ring of scalar role mixing and P^{+-1} child calls, one phase child cannot suffice: the required correction matrix [[P,0],[-1,P]] has determinant P^2, whereas a single child changes the determinant by P or P^-1. This is a restricted obstruction, not a proof against additional address transforms or a different primitive. The two-stage construction is special because its first output is already -Fy, so only one rank-one correction remains.
