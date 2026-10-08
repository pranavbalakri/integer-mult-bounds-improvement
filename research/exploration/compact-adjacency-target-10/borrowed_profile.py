"""Conditional histogram for the new borrowed-source word; not a global certificate."""

if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')

from collections import Counter
from fractions import Fraction as Q
from importlib.util import spec_from_file_location,module_from_spec
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]
spec=spec_from_file_location('exact_moment',RESEARCH/'independent/complex-twostage/cert.py')
EXP=module_from_spec(spec);spec.loader.exec_module(EXP)

def transform(record):
    h=record['h'];m=record['m'];v=record['v'];N=record['N']
    hist=Counter({int(w):n for w,n in record['hist'].items()})
    def inner(r):return [1]*r if 2*r<=h else [1]*(h-r)+[2*r-h]
    removed=Counter([m-2*h,h,1]+inner(h-1))
    assert sum(w*n for w,n in removed.items())==m
    for w,n in removed.items():hist[w]-=2*N*n
    assert all(n>=0 for n in hist.values())
    hist={w:n for w,n in hist.items() if n}
    W=record['W']-2*N;s=sum(w*n for w,n in hist.items())
    assert s==record['s']-2*N*m
    deficit=W*m-s;assert deficit==record['W']*m-record['s']
    a,slack=EXP.cert(hist,W,m,grid=10**12)
    assert a>Q(record['saving']) and slack>0
    return dict(h=h,m=m,v=v,N=N,W=W,s=s,R=record['R']-v,
                rank_deficit=deficit,hist={str(w):n for w,n in sorted(hist.items())},
                old_saving=record['saving'],candidate_saving=str(a),
                improvement_factor=float(a/Q(record['saving'])),exact_moment_slack=str(slack),
                status='Conditional histogram only: new borrowed-source tensor transfer requires review.')

if __name__=='__main__':
    old=json.loads((RESEARCH/'independent/bit-improvement/certificate-h23.json').read_text())
    record=transform(old)
    (HERE/'borrowed-profile-h23.json').write_text(json.dumps(record,indent=2)+'\n')
    print({k:v for k,v in record.items() if k not in ('hist','exact_moment_slack')})
