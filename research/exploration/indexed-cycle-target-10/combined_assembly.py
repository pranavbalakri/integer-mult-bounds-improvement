"""Exploratory exact macro-stack arithmetic, not a multiplication theorem.

Combines credited PR84 bit profiles with this repository's h18 complex word
and linear precision guard. The combined global interface is UNREVIEWED.
This uses the enlarged-precision macro stack, not PR84's refined 47-row stack.
Written with OpenAI Codex assistance. Apache-2.0; see README and UPSTREAM-NOTICE.
"""
from argparse import ArgumentParser
from fractions import Fraction as Q
from pathlib import Path
import json
import sys
from independent_moment import profile, moment

HERE = Path(__file__).resolve().parent
AC = Q(59292472, 10**12)
PR84 = Q(408390793141, 7812500000000000)


def bridge(bit_width=132466108):
    axes = {}
    for name, m, child, width in (
        ('bit', 575, 529, bit_width), ('complex', 324, 306, 32011680)
    ):
        degree = 1
        while m**degree <= 2*child**degree:
            degree += 1
        axes[name] = dict(m=m, W=width, maxchild=child,
                          halving_degree=degree, wire_bits=width.bit_length())
    coefficient = sum(x['halving_degree']*x['wire_bits'] for x in axes.values())
    assert coefficient == 568
    row_degree = 1200
    row_gap = row_degree-Q(51,25)*coefficient
    assert row_gap == Q(1032,25) > 0
    return dict(axes=axes, row_coefficient=coefficient, row_degree=row_degree,
                row_degree_gap=row_gap, suffix_slope=4*row_degree,
                guard='128*(d+1)',
                scope='Arithmetic template only; new combined bridge has not been globally integrated')


def evaluate(name, ab, eps, delta, gap, kappa, power=Q(3), *, bit_profile=None):
    beta, x, c1, y = Q(1,16), power-1, Q(1), eps
    tau, sigma = 1-ab, 1-AC
    internal = tau+(1-beta)*max(sigma-tau,Q(0))
    leaf = sigma+beta*(1-sigma)
    lam = max(tau,sigma,internal)+gap
    lamp = max(lam,leaf)+gap
    constraints = {
        'bit saving positive': ab,
        'complex saving above bit': AC-ab,
        'leaf saving above bit': (1-beta)*AC-ab,
        'lambda above tau sigma internal': lam-max(tau,sigma,internal),
        'lambda_prime above lambda leaf': lamp-max(lam,leaf),
        'lambda_prime below one': 1-lamp,
        'beta in (0,1)': min(beta,1-beta),
        'guard exponent': power-eps*c1,
        'Gaussian precision exponent': power-2*eps-y,
        'reservations': lamp-(2*eps-1)/eps,
        'record regime': 1-eps,
        'delta in (0,1/8)': min(delta,Q(1,8)-delta),
    }
    assert y == eps and power > 1
    margins = {
        'prefix': 1-eps,
        'butterflies': eps*(1-lamp),
        'fine fields Gaussian chirps reserved axes': 1-eps-power*delta,
        'packed products': eps,
        'CRT reversal': max(eps,1-eps)*(1-tau),
    }
    assert all(v > 0 for v in constraints.values())
    assert min(margins.values()) > kappa > 0
    _, upper = moment(profile() if bit_profile is None else bit_profile,ab)
    assert upper < 1
    return dict(name=name, bit_saving=ab, complex_saving=AC, epsilon=eps,
                delta=delta, beta=beta, lambda_gap=gap, lambda_value=lam,
                lambda_prime=lamp, x=x, y=y, C1=c1, kappa=kappa,
                precision_power=power, precision='b_prime=128*ceil(b^P); p=768*ceil(b^P)',
                mu=(1-eps)/(2*power), axis_width_power=1-eps,
                Gaussian_uniform_lower_bound_coefficient=Q(15,8),
                Gaussian_lower_bound_power=power-3*eps,
                constraints=constraints, margins=margins,
                minimum_margin=min(margins.values()),
                strict_surplus=min(margins.values())-kappa,
                bit_moment_gap=1-upper, bit_moment_gap_float=float(1-upper))


def verify_complex():
    research = HERE.parents[1]
    sys.path.insert(0,str(research/'independent/complex-twostage'))
    from hist import build, label_dims
    from cert import cert
    from frames import Checker
    data = build(18,copied=True)
    assert (data['m'],data['W'],data['s'],data['maxrank']) == (324,32011680,10371619488,306)
    saving, gap = cert(data['hist'],data['W'],data['m'])
    assert saving == AC
    checker = Checker(data['c'])
    assert {n:len(checker.label(n)) for n in data['c'].active} == label_dims(data['c'])
    frames = checker.run()
    assert frames['bad'] == 0
    sys.path.insert(0,str(research/'independent/guard-improvement'))
    from semantic_guard import certify
    from audit import matrix_identity
    guard = certify(18)
    assert guard['C1'] == 1 and guard['layer_guard'] == '128(d+1)'
    matrix = matrix_identity(18)
    # Keep the guard statement separate from the guard verifier's original
    # cubic-precision example; each exploratory precision is checked above.
    guard = {k:guard[k] for k in ('h','m','roles','G','B','prefix_charge',
             'contraction','recursive_guard_coefficient','layer_guard','C1','assumptions')}
    return dict(complex_moment_gap=gap, frames=frames, guard=guard, all_pairs=matrix)


def combine(check_complex=False):
    result = dict(
        status='Exact macro-stack arithmetic PASS; combined global interface UNREVIEWED',
        source_bit_sha='88ca39571907343a49e97f328971ec7bcd26fbfd',
        bit_credit='Chafik Boukhalfa and the credited PR84 contributors; see UPSTREAM-NOTICE',
        bridge=bridge(), profiles=[
            evaluate('larger saving',PR84,Q(49997,50000),Q(1,10**6),Q(1,10**11),Q(5227,10**8)),
            evaluate('wider analytic slack',PR84,Q(71,100),Q(9,100),Q(1,10**8),Q(37,10**6)),
            evaluate('wider slack and recurrence backoff',Q(522,10**7),Q(71,100),Q(9,100),Q(1,10**8),Q(37,10**6)),
            evaluate('lower precision with wider slack and recurrence backoff',Q(522,10**7),Q(71,100),Q(9,100),Q(1,10**8),Q(37,10**6),Q(43,20)),
        ], caveats=[
            'These are exact inequalities in the enlarged-precision macro stack, not PR84\'s 47-row refined assembly.',
            'Finite words, frames and moments have separate replay evidence; that evidence does not close the global interface.',
            'The global ordered residual compiler, routing, prime, recovery and tape hypotheses remain inherited.',
            'No finite operational crossover or practical runtime claim is made.',
            'All candidate exponents are below the credited PR84 headline; the proposed benefit is lower overhead.',
        ])
    if check_complex:
        result['local_complex_validation'] = verify_complex()
    return result


def serial(x):
    if isinstance(x,Q): return str(x)
    if isinstance(x,dict): return {k:serial(v) for k,v in x.items()}
    if isinstance(x,list): return [serial(v) for v in x]
    return x


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--check-complex',action='store_true',help='Recompute the h18 moment, all label incidences, all ordered scalar pairs and guard')
    args = parser.parse_args()
    print(json.dumps(serial(combine(args.check_complex)),indent=2))
