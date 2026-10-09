"""Zero-reference cleanups at existing producer frames.

A source X visits a nested chain of source-support frames, then its primary
point star.  Every such detour must be charged by a subsequent profile audit.
This checker makes no moment or exponent claim.
"""
if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')

from collections import Counter, defaultdict
from pathlib import Path
import argparse
import hashlib
import json
import reference_plans as rp


def select(M, P, order='small-first'):
    assert not P['refs']
    G=M['G']; C=M['C']; h=M['h']
    groups=defaultdict(list)
    for s in P['anchors']:
        n=C['hold'][s][-1]
        if G['dn'][n] < h-1:
            groups[n].append(s)
    visits=[set() for _ in range(M['v'])]
    chosen={}
    def compatible(i,n):
        a=G['sup'][n]
        return all(a & ~G['sup'][other] == 0 or G['sup'][other] & ~a == 0
                   for other in visits[i])
    def key(n):
        if order=='large-first':return -G['dn'][n],-len(groups[n]),n
        if order=='multiplicity':return -len(groups[n]),G['dn'][n],n
        return G['dn'][n],-len(groups[n]),n
    for n in sorted(groups,key=key):
        # A group may have several retained operand roles; incompatible rows
        # are left at their already valid point-star cleanup.
        for s in groups[n]:
            sources=list(rp.old.bits(M['signal'][s]))
            if all(compatible(i,n) for i in sources):
                chosen[s]=n
                for i in sources:visits[i].add(n)
    Q=dict(P)
    Q['local_nodes']=chosen
    Q['source_chains']=[sorted(chain,key=lambda n:(G['dn'][n],G['sup'][n].bit_count(),n))
                        for chain in visits]
    return Q


def audit(M,P,full=False):
    G=M['G'];C=M['C'];h=M['h']
    assert not P['refs']
    chains=P['source_chains'];chosen=P['local_nodes']
    assert len(chains)==M['v']
    for i,chain in enumerate(chains):
        previous=1<<i
        c=P['primary'][i]
        for n in chain:
            label=G['sup'][n]
            assert not previous&~label
            assert all(c in M['triples'][j] for j in rp.old.bits(label))
            previous=label
    for s,n in chosen.items():
        assert s in P['anchors'] and s not in C['out']
        assert n==C['hold'][s][-1]
        assert G['dn'][n]<h-1
        for i in rp.old.bits(M['signal'][s]):
            assert n in chains[i]
    # The literal cleanup schedule is increasing in dimension, then role ID.
    # A source may have multiple equal-rank nested supports: the represented
    # nondegenerate spaces then coincide, so no frame motion is required.
    clean_order=sorted(chosen,key=lambda s:(G['dn'][chosen[s]],s))
    for i in range(M['v']):
        observed=[chosen[s] for s in clean_order if M['signal'][s]>>i&1]
        ranks=[G['dn'][n] for n in observed]
        assert ranks==sorted(ranks)
    receipt=rp.audit_plan(M,P,full)
    hist=Counter(G['dn'][n] for n in chosen.values())
    sourcehist=Counter()
    source_local_edges=0
    for i,chain in enumerate(chains):
        if i not in P['used_primary']:
            assert not chain
            continue
        d=1
        for n in chain:
            nd=G['dn'][n]
            if nd>d:sourcehist[nd-d]+=1;source_local_edges+=1
            d=nd
        if h-1>d:sourcehist[h-1-d]+=1;source_local_edges+=1
        sourcehist[1]+=1;source_local_edges+=1
    assert sum(w*c for w,c in sourcehist.items())==(h-1)*len(P['used_primary'])
    payload=dict(local_nodes=sorted(chosen.items()),source_chains=chains,
                 base_plan_sha256=receipt['plan_sha256'])
    receipt.update(cleaned_at_existing_last_frame=len(chosen),
                   remaining_point_star_cleanups=len(P['anchors'])-len(chosen),
                   direct_cleanup_by_last_rank=dict(sorted(hist.items())),
                   source_data_local_rank_histogram=dict(sorted(sourcehist.items())),
                   source_data_local_edges=source_local_edges,
                   source_paths_are_nested_existing_support_frames=True,
                   direct_cleanup_order_compatible_with_every_source_path=True,
                   hierarchy_sha256=hashlib.sha256(json.dumps(payload,sort_keys=True,
                     separators=(',',':')).encode()).hexdigest(),
                   status='Exact scalar and nested-support schedule only. All source detours and merged exits need a fixed-basis child profile.')
    return receipt


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('h',nargs='*',type=int,default=[5,6,7,8,23])
    parser.add_argument('--full-max',type=int,default=8)
    parser.add_argument('--order',action='append',choices=['small-first','large-first','multiplicity'])
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('hierarchy-results.json'))
    args=parser.parse_args()
    result=[]
    for h in args.h:
        M=rp.model(h);P=rp.choose(M,[min(t) for t in M['triples']])
        for order in args.order or ('small-first','large-first','multiplicity'):
            Q=select(M,P,order)
            row=audit(M,Q,h<=args.full_max);row['selection_order']=order
            result.append(row)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
