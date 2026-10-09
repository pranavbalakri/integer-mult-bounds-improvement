"""Exact local obstruction to cap18 with the old <=18 frames and order fixed."""
if not __debug__:
    raise RuntimeError('Run without -O')

from pathlib import Path
import argparse
import hashlib
import json
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'deferred-truncation-target-10'))
from check_truncation import PINS


def determinant(A):
    """Integer Bareiss elimination, checking every division exactly."""
    A=[row[:] for row in A];n=len(A);previous=1;sign=1
    assert n and all(len(row)==n for row in A)
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if A[i][k]),None)
        if pivot is None:return 0
        if pivot!=k:A[k],A[pivot]=A[pivot],A[k];sign=-sign
        p=A[k][k]
        for i in range(k+1,n):
            q=A[i][k]
            for j in range(k+1,n):
                numerator=p*A[i][j]-q*A[k][j]
                assert numerator%previous==0
                A[i][j]=numerator//previous
            A[i][k]=0
        previous=p
    return sign*A[-1][-1]


def run(source):
    for name,digest in PINS.items():
        assert hashlib.sha256((source/name).read_bytes()).hexdigest()==digest,name
    sys.path.insert(0,str(source/'independent/deferred-readout'))
    import deferred as dr
    from linalg import Echelon,Q31,null_exact,rank_mod
    S=dr.Schedule(*dr.load(23));cov=S.adjoint();position={s:i for i,s in enumerate(S.readout)}
    last={}
    for s in S.readout:
        if S.f[s]<=18:
            for t in cov[s]:last[t]=s
    blocked=[]
    for s in S.readout:
        if S.f[s]>18:
            predecessors=sorted({last[t] for t in cov[s] if t in last})
            rank=rank_mod([r for p in predecessors for r in S.sigma[p]],Q31)
            if rank>18:blocked.append(dict(slot=s,old_dimension=S.f[s],
                targets=sorted(cov[s]),predecessors=predecessors,rank_lower_bound=rank))
    assert len(blocked)==14
    witness=next(row for row in blocked if row['slot']==23806)
    s=witness['slot'];preds=witness['predecessors']
    assert preds==[24415,25074] and S.f[s]==19
    assert all(S.f[p]==18 and position[p]<position[s] for p in preds)
    E=Echelon(Q31);selected=[]
    for p in preds:
        for i,row in enumerate(S.sigma[p]):
            if E.insert(row):selected.append((p,i,row))
    assert len(E)==19
    columns=E.piv
    minor=[[row[j] for j in columns] for p,i,row in selected]
    det=determinant(minor)
    assert det and det%Q31
    # An exact 19-dimensional container proves the matching upper bound.
    H=S.sigma[s];assert len(H)==19 and rank_mod(H,Q31)==19
    null=null_exact(H,S.h)
    assert all(sum(x*y for x,y in zip(row,z))==0
               for p in preds for row in S.sigma[p] for z in null)
    targets=[]
    for t,coefficient in sorted(cov[s].items()):
        p=last[t]
        assert p in preds and coefficient==1 and cov[p][t]==1
        targets.append(dict(target_id=t,target=list(S.trip[t]),predecessor_slot=p,
                            predecessor_dimension=S.f[p],readout_coefficient=coefficient))
    return dict(source_repository='https://github.com/Swapnil-jain/integer-mult-kappa',
      source_commit='741e7aa078392553815df7926ee17ac5e25a8c38',source_sha256=PINS,
      h=S.h,cap=18,blocked_slot=s,old_slot_dimension=19,required_targets=targets,
      predecessor_union_rank_exact=19,predecessor_intersection_dimension=17,
      original_readout_order_retained=True,unchanged_predecessors_at_dimensions_at_most_18=True,
      applies_to_integer_and_F2_readout_supports=True,
      nonzero_integer_minor_dimension=19,minor_columns=columns,
      minor_rows_from_slot_and_basis_row=[[p,i] for p,i,row in selected],
      exact_integer_minor=minor,exact_integer_determinant=det,
      old_slot_annihilator=null,all_local_obstructions=blocked,
      conclusion='No <=18-dimensional replacement for slot23806 can contain both unchanged preceding18-planes. This is an obstruction to this frame-capping ansatz and retained order, not to other algorithms or changed lower frames.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,default=HERE/'cap18-obstruction.json')
    args=parser.parse_args()
    result=run(args.source.resolve())
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in
        ('source_sha256','exact_integer_minor','old_slot_annihilator','all_local_obstructions')},indent=2))
