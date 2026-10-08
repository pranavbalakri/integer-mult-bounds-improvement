"""Experimental retained-total word borrowing the source data roles.
All source, side, centre and dirty coefficients are replayed exactly for small h.
The optional large histogram comparison remains conditional on the new word's
transport through the global tensor/recursion interface.
"""

if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')

from collections import Counter
from pathlib import Path
import json, sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'independent/two-stage-bit'))
import rtgm, side_chains


def compile_borrowed(G, replay=True):
    h=G['h'];v=len(G['trip']);args=G['args'];act=G['active']
    users={n:[] for n in act}
    for n in sorted(act):
        if args[n]:
            for pos,k in enumerate(args[n]):users[k].append(('gate',n,pos))
    for target,n in sorted(G['outputs'].items()):users[n].append(('side',target))
    for c,n in sorted(G['retained'].items()):users[n].append(('centre',c))
    edge={};mid=[];history={i:[] for i in range(v)};aux=[];side=[];centres=[]
    def alloc():
        s=2*v+len(aux);aux.append(s);history[s]=[];return s
    def xor(a,b,n):
        assert a!=b;mid.append((a,b,n));history[a].append(n);history[b].append(n)
    for n in sorted(act):
        if args[n]:
            a,b=edge[n,0],edge[n,1];xor(a,b,n);pivot=a
        else:
            assert 1<=n<=v;pivot=n-1;history[pivot].append(n)
        for k,u in enumerate(users[n]):
            s=pivot if k==0 else alloc()
            if k:xor(s,pivot,n)
            if u[0]=='gate':edge[u[1],u[2]]=s
            elif u[0]=='side':side.append((s,u[1]))
            else:centres.append((s,u[1]))
    assert len(aux)==rtgm.compile_(G)['roles']-v
    # Initial and middle frame histories. Existing source-frame theorem applies.
    rk=Counter();rankmass=0
    terminal=dict((s,('side',t)) for s,t in side)
    terminal.update((s,('centre',c)) for s,c in centres)
    assert len(terminal)==len(side)+len(centres)
    for s,hs in history.items():
        previous=(1<<s) if s<v else 0
        dim=1 if s<v else 0
        for n in hs:
            assert not previous&~G['sup'][n]
            nd=G['dn'][n];assert nd>=dim
            if nd>dim:rk[nd-dim]+=1
            previous=G['sup'][n];dim=nd
        if s in terminal:
            typ,label=terminal[s]
            if typ=='side':
                _,t=label
                assert all(len(set(G['trip'][i])&set(t))==1 for i in range(v) if previous>>i&1)
            else:
                assert previous==sum(1<<i for i,t in enumerate(G['trip']) if label in t)
            assert dim<=h-1
            if dim<h-1:rk[h-1-dim]+=1
            dim=h-1
        if dim<h:rk[h-dim]+=1
    assert sum(r*n for r,n in rk.items())==h*len(aux)+(h-1)*v
    # At all centre reads Y is still at D0. Side reads follow, then cleanup at F.
    Jcentre=[(v+i,s) for s,c in centres for i,t in enumerate(G['trip']) if c in t]
    Jside=[(v+G['trip'].index(t),s) for s,(_,t) in side]
    def check_scatter_order(events):
        yranks=[0]*v
        for kind,target in events:
            if kind=='centre':
                assert yranks[target]==0, 'A centre read would return Y to D0'
            else:
                assert yranks[target] in (0,h-1)
                yranks[target]=h-1
        assert yranks==[h-1]*v
    centre_events=[('centre',a-v) for a,b in Jcentre]
    side_events=[('side',a-v) for a,b in Jside]
    check_scatter_order(centre_events+side_events)
    try:check_scatter_order(side_events+centre_events)
    except AssertionError:pass
    else:raise AssertionError('The deliberately reversed scatter order was accepted')
    summary=dict(h=h,v=v,old_aux_roles=len(aux)+v,new_aux_roles=len(aux),borrowed_data_roles=v,
                 middle_cnot_count=len(mid),centre_reads=len(Jcentre),side_reads=len(Jside),
                 all_middle_histories_monotone=True,centre_reads_before_side_reads=True,
                 negative_side_before_centres_rejected=True,
                 middle_rank_histogram=dict(sorted(rk.items())))
    if replay:
        count=2*v+len(aux);initial=[1<<i for i in range(count)];values=initial[:]
        for a,b,n in mid:values[a]^=values[b]
        image=[0]*v
        for a,b in Jcentre+Jside:image[a-v]^=values[b]
        xmask=(1<<v)-1;ymask=xmask<<v
        assert [r&xmask for r in image]==[1<<i for i in range(v)]
        assert all(not r&ymask for r in image)
        # D is the unwanted arbitrary-scratch contribution, cancelled at D0.
        early=[]
        for i,row in enumerate(image):
            row>>=2*v
            while row:
                p=(row&-row).bit_length()-1;early.append((v+i,2*v+p));row&=row-1
        word=early+[(a,b) for a,b,n in mid]+Jcentre+Jside+[(a,b) for a,b,n in reversed(mid)]
        for mode in ('forward','reverse','same_order_transpose','reverse_transpose','reverse_bank_swap'):
            values=initial[:]
            gates=list(reversed(word)) if mode in ('reverse','reverse_transpose','reverse_bank_swap') else word
            def bank_swap(a):
                if a<v:return a+v
                if a<2*v:return a-v
                return a
            for a,b in gates:
                if mode in ('same_order_transpose','reverse_transpose'):a,b=b,a
                if mode=='reverse_bank_swap':a,b=bank_swap(a),bank_swap(b)
                values[a]^=values[b]
            if mode in ('forward','reverse'):
                assert values[:v]==initial[:v]
                assert values[v:2*v]==[initial[v+i]^initial[i] for i in range(v)]
            else:
                assert values[:v]==[initial[i]^initial[v+i] for i in range(v)]
                assert values[v:2*v]==initial[v:2*v]
            assert values[2*v:]==initial[2*v:]
        summary.update(early_cnot_count=len(early),word_length=len(word),
                       full_dirty_basis_forward_and_same_order_transpose_checked=True,
                       reverse_and_reverse_bank_swap_also_checked=True,
                       dirty_basis_dimension=count)
        # A sparse realization of D uses fresh zero temporary source pivots.
        # These are private storage at D0, not arbitrary input roles in W.
        def temporary(a):return count+a if a<v else a
        temp_mid=[(temporary(a),temporary(b)) for a,b,n in mid]
        temp_J=[(a,temporary(b)) for a,b in Jcentre+Jside]
        factored_early=temp_mid+temp_J+list(reversed(temp_mid))
        middle_word=[(a,b) for a,b,n in mid]+Jcentre+Jside+[(a,b) for a,b,n in reversed(mid)]
        for transpose in (False,True):
            values=initial[:]+[0]*v
            for a,b in factored_early:
                if transpose:a,b=b,a
                values[a]^=values[b]
            if not transpose:assert not any(values[count:])
            # The inverse-transpose may leave private T=C^T Y; clear it at D0.
            values[count:]=[0]*v
            for a,b in middle_word:
                if transpose:a,b=b,a
                values[a]^=values[b]
            assert values[2*v:count]==initial[2*v:]
            assert not any(values[count:])
            if transpose:
                assert values[:v]==[initial[i]^initial[v+i] for i in range(v)]
                assert values[v:2*v]==initial[v:2*v]
            else:
                assert values[:v]==initial[:v]
                assert values[v:2*v]==[initial[v+i]^initial[i] for i in range(v)]
        summary.update(fresh_zero_source_temporaries=v,
                       factored_early_gate_count=len(factored_early),
                       factored_word_gate_count=len(factored_early)+len(middle_word),
                       factored_early_public_basis_both_orientations_checked=True)
        # Negative control: omit a nonzero early cancellation gate.
        assert early
        values=initial[:]
        for a,b in word[1:]:values[a]^=values[b]
        assert values[v:2*v]!=[initial[v+i]^initial[i] for i in range(v)]
        summary['negative_missing_dirty_cancellation_rejected']=True
    return summary


if __name__=='__main__':
    optimized='--optimized' in sys.argv
    dims=list(map(int,[x for x in sys.argv[1:] if x!='--optimized'])) or ([23] if optimized else [5,6,7,8])
    if optimized:assert dims==[23], 'The optimized receipt is pinned to h23'
    result=[]
    for h in dims:
        if optimized:
            sys.path.insert(0,str(ROOT/'independent/bit-improvement'))
            import reassociate
            old,G,C=reassociate.histogram(h)
        else:G=rtgm.build(h,retain='search')
        record=compile_borrowed(G)
        if optimized:
            previous=Counter()
            for slot in range(C['roles']):
                chain=rtgm.chain_dims(G,C,slot)
                previous.update(b-a for a,b in zip(chain,chain[1:]))
            current=Counter(record['middle_rank_histogram'])
            assert previous==current+Counter({1:len(G['trip'])})
            from borrowed_profile import transform
            old['saving']='36943733/1000000000000' if h==23 else '0'
            candidate=transform(old)
            v=record['v'];N=v*v;m=h*h;H=Counter()
            def add(w,n):
                if n:H[w]+=n
            def inner(r):return [1]*r if 2*r<=h else [1]*(h-r)+[2*r-h]
            for _ in range(2):
                add(m-2*h,v*record['new_aux_roles']);add(h,v*record['new_aux_roles'])
                for r,n in current.items():
                    for w in inner(r):add(w,v*n)
                for w in inner(h-1):add(w,v*h)
            add(m-4*h+2,2*N);add(h-2,4*N);add(1,2*N*(h+2)+N)
            assert dict(H)=={int(w):n for w,n in candidate['hist'].items()}
            record.update(parent_producer='paired-tree exact-support reassociation',
                          independently_rebuilt_histogram_matches=True,
                          candidate_saving=candidate['candidate_saving'],W=candidate['W'],
                          rank_sum=candidate['s'],rank_deficit=candidate['rank_deficit'])
        result.append(record)
    payload=result[0] if optimized and dims==[23] else result
    print(json.dumps(payload,indent=2))
    name='borrowed-sources-h23-audit.json' if optimized and dims==[23] else 'borrowed-sources-audit.json'
    Path(__file__).with_name(name).write_text(json.dumps(payload,indent=2)+'\n')
