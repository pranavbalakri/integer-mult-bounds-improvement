"""Exact stack margins for two choices trading exponent against overhead.

Assumes a separately established C1=1 phase guard and complex h16 witness.
The Gaussian-width constraint at x=2 can use the explicit fixed digit prefactor
b'=128*b^3, p=6*b', b=ceil(log2(n)). See construction.tex.
"""
from fractions import Fraction as Q
from pathlib import Path
import json

AB = Q(36943733, 10**12)
AC = Q(39650268, 10**12)
BASELINE = Q(3666565558019, 10**17)


def evaluate(name, eps, delta, kappa, gap, x=Q(2), beta=Q(1, 16), c1=Q(1)):
    tau, sigma = 1-AB, 1-AC
    chi = tau + (1-beta)*max(sigma-tau, 0)
    leaf = sigma + beta*(1-sigma)
    lam = max(tau, sigma, chi) + gap
    lamp = max(lam, leaf) + gap
    y = eps
    poly = (1+x)*delta
    constraints = {
        'lambda above tau sigma chi': lam-max(tau,sigma,chi),
        'lambda_prime above lambda leaf': lamp-max(lam,leaf),
        'lambda_prime below one': 1-lamp,
        'beta in (0,1)': min(beta,1-beta),
        'guard': 1+x-eps*c1,
        'Gaussian precision': 1+x-2*eps-y,
        'coupling y>=eps': Q(1) if y>=eps else Q(-1),
        'reservations': lamp-(2*eps-1)/eps,
        'record regime': 1-eps,
        'delta in (0,1/8)': min(delta,Q(1,8)-delta),
    }
    margins = {
        'prefix': 1-eps,
        'butterflies': eps*(1-lamp),
        'fine fields Gaussian chirps reserved axes': 1-eps-poly,
        'packed products': eps,
        'CRT reversal': max(eps,1-eps)*(1-tau),
    }
    bad = [k for k,v in constraints.items() if v<=0]
    bad += [k for k,v in margins.items() if v<=kappa]
    if kappa<=BASELINE: bad.append('does not beat comparator')
    return dict(name=name, ok=not bad, bad=bad, exponent_gap=str(gap), a_b=str(AB), a_c=str(AC), eps=str(eps), delta=str(delta),
                x=str(x), beta=str(beta), y=str(y), C1=str(c1), lambda_value=str(lam),
                lambda_prime=str(lamp), kappa=str(kappa), mu=str((1-eps)/(2*(1+x))),
                precision='b_prime=128*b^3; p=768*b^3',
                constraints={k:str(v) for k,v in constraints.items()},
                margins={k:str(v) for k,v in margins.items()},
                gap=str(min(margins.values())-kappa),
                improvement_over_baseline=str(kappa-BASELINE), relative_improvement=str(kappa/BASELINE-1))


def options():
    result = [
        evaluate('higher saving', Q(24999,25000), Q(1,10**6), Q(36942,10**9), Q(1,10**11)),
        evaluate('lower overhead', Q(993,1000), Q(1,500), Q(3668,10**8), Q(1,10**9)),
    ]
    assert all(r['ok'] for r in result), result
    return result


if __name__ == '__main__':
    data = options()
    Path(__file__).with_name('assembly-options.json').write_text(json.dumps(data,indent=2)+'\n')
    for option in data:
        print(option['name'], 'kappa', option['kappa'], 'gap', option['gap'])
