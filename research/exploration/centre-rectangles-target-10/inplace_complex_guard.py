"""Positive exact envelopes for the proposed in-place complex scalar word.

Includes dense early cancellation using |D| <= |J| M and its transpose,
without materializing D. This does not certify the lifted child histogram.
"""
if not __debug__:
    raise RuntimeError("This exact checker requires assertions; do not run Python with -O.")

from fractions import Fraction as Q
from pathlib import Path
import json
import sys

RESEARCH=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(RESEARCH/'independent/complex-twostage'))
from producer import NStar3,compile_roles


def certify(h):
    c=NStar3(h);kk=compile_roles(c);R=kk['size'];v=len(c.triples)
    sources=set(kk['src'].values());tid={t:i for i,t in enumerate(c.triples)}
    M=[]
    for _,ins,outs in kk['gates']:
        M.extend((ins[0],s) for s in ins[1:])
        M.extend((s,ins[0]) for s in outs[1:])
    J=[]
    for i,t in enumerate(c.triples):
        last=h-1 in t
        J.append((i,kk['rout'][('*',)],Q(5-h if last else 2,2)))
        points=[j for j in range(h-1) if j not in t] if last else t
        J.extend((i,kk['rout'][('E',j)],Q(1 if last else -1,2)) for j in points)
    for k,s in kk['pout'].items():
        t,n,coef=c.pieces[k]
        J.append((tid[t],s,Q(coef,2)))
    # M has only additions, so applying it to all ones gives its row sums.
    mrow=[Q(1)]*R
    for a,b in M:mrow[a]+=mrow[b]
    # These include the source columns too, making the bound conservative.
    drow=[Q(0)]*v;dcol=[Q(0)]*R
    for a,b,c in J:drow[a]+=abs(c)*mrow[b];dcol[b]+=abs(c)
    for a,b in reversed(M):dcol[b]+=dcol[a]
    peaks=[]
    for transpose in (False,True):
        if not transpose:
            x=[Q(1)]*R;y=[1+z for z in drow]
        else:
            x=[Q(1) if i in sources else 1+dcol[i] for i in range(R)];y=[Q(1)]*v
        peak=max(x+y)
        for a,b in M:
            if transpose:a,b=b,a
            x[a]+=x[b];peak=max(peak,x[a])
        for a,b,c in J:
            if not transpose:y[a]+=abs(c)*x[b];peak=max(peak,y[a])
            else:x[b]+=abs(c)*y[a];peak=max(peak,x[b])
        for a,b in reversed(M):
            if transpose:a,b=b,a
            x[a]+=x[b];peak=max(peak,x[a])
        # Forming a numerator before division by two costs at most this factor.
        peaks.append(2*peak)
    G=2*peaks[0]*peaks[1] # two disjoint-bank stages and endpoint addition
    B=2 # one dyadic denominator bit per stage
    bits=0
    while G>2**bits:bits+=1
    Jcharge=max(B,bits)+2
    conservative_J=max(64,Jcharge)
    m=h*h;rho=Q(h-1,h)
    C=max(8,2*h+(conservative_J+h-1)//h)
    assert C*(1-rho)>=2+Q(conservative_J,m)
    return dict(h=h,roles_before=R,roles_after=R-v,
                forward_prefix_G=str(peaks[0]),inverse_transpose_prefix_G=str(peaks[1]),
                two_stage_endpoint_G=str(G),two_stage_denominator_bits=B,
                ceil_log2_G=bits,constant_prefix_charge=conservative_J,
                recursive_guard_coefficient=C,
                same_128_layer_guard_numeric_check=C+10<=128,
                scope='Exact scalar envelopes; conditional on the proposed word and inherited whole-residual width interface')


if __name__=='__main__':
    print(json.dumps([certify(h) for h in map(int,sys.argv[1:] or ['8','16','18'])],indent=2))
