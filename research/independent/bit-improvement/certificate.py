"""Exact exponential-moment certificate for the reassociated bit producer."""
import importlib.util
import json
from pathlib import Path
from fractions import Fraction as Q
import reassociate

SOURCE = Path(__file__).resolve().parents[1] / 'complex-twostage' / 'cert.py'
SPEC = importlib.util.spec_from_file_location('complex_exact_moment', SOURCE)
EXP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXP)


def certificate(h=23):
    network, graph, compiled = reassociate.histogram(h)
    saving, slack = EXP.cert(network['hist'], network['W'], network['m'], grid=10**12)
    weights = [(Q(w*n, network['W']*network['m']), EXP.ln_up(Q(network['m'], w)))
               for w, n in network['hist'].items()]
    assert saving > Q(36667, 10**9)
    assert slack > 0
    assert EXP.moment(weights, saving + Q(1, 10**12)) >= 1
    result = {key: value for key, value in network.items() if key != 'hist'}
    result.update(hist={str(w): n for w, n in sorted(network['hist'].items())},
                  saving=str(saving), exact_moment_slack=str(slack),
                  reassociated_nodes=graph['reassociated_nodes'],
                  proof_scope='conditional on the inherited bit network and common flag-basis contracts')
    return result


if __name__ == '__main__':
    output = certificate()
    target = Path(__file__).with_name('certificate-h23.json')
    target.write_text(json.dumps(output, indent=2) + '\n')
    print('saving', output['saving'])
    print('W', output['W'], 'm', output['m'], 'roles', output['R'])
    print('exact moment slack', output['exact_moment_slack'])
