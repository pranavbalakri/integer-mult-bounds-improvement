"""Exact assembly inequalities for the borrowed-source PR84 extension.

This does not certify the inherited global interfaces. The finite new word,
its lifted frames and the profile have separate audits. Main notes and the
main multiplication certificate remain unchanged pending wider review.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
import sys

if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'indexed-cycle-target-10'))
import combined_assembly as macro
from borrowed_indexed import profile


def check():
    p = profile()
    profile_for_moment = dict(m=p['m'], W=p['W'],
                             rows={int(r): n for r, n in p['child_multiplicities'].items()})
    saving = Q(p['bit_saving'])
    results = [
        macro.evaluate('larger conditional saving', saving, Q(9999, 10000),
                       Q(1, 10**6), Q(1, 10**11), Q(5487, 10**8),
                       bit_profile=profile_for_moment),
        macro.evaluate('lower precision and wider analytic margins', saving,
                       Q(71, 100), Q(9, 100), Q(1, 10**8), Q(389, 10**7),
                       Q(43, 20), bit_profile=profile_for_moment),
    ]
    assert all(item['kappa'] < Q(1, 1024) for item in results)
    return dict(status='Exact candidate assembly inequalities PASS; complete global theorem remains conditional and unreviewed',
                borrowed_bit_word='borrowed_indexed.py; borrowed-sources.tex',
                original_bit_source_commit='88ca39571907343a49e97f328971ec7bcd26fbfd',
                credit='New wrapper around Chafik Boukhalfa and contributors\' credited PR84 producer; see ../indexed-cycle-target-10/UPSTREAM-NOTICE',
                bridge=macro.bridge(p['W']), profiles=results,
                inherited_complex='Existing h18 two-stage word and linear guard, unchanged',
                target_met=False,
                limitations=[
                    'The ordered residual compiler, common rational flags, compact routing, prime selection, exact recovery and fixed-tape interfaces remain assumptions.',
                    'The new wrapper has separate full-dirty-basis replay and a lifted-frame argument, rather than a complete all-size executable implementation.',
                    'The lower-precision variant deliberately has a smaller kappa; the two parameter choices must not be combined.',
                    'Role counts and precision requirements are not running times or operational crossover thresholds.',
                    'Main certificate and open LaTeX note have not been changed.',
                ])


if __name__ == '__main__':
    result = macro.serial(check())
    (HERE/'borrowed-assembly-results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
