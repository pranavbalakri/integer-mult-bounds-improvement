"""Scalar, framed-word, moment and semantic-guard checks for the new producer."""
from collections import Counter
from fractions import Fraction as Q
from math import comb, exp, log, floor, lcm
from pathlib import Path
from hashlib import sha256
import argparse
import json
from producer import Producer, ROOT, imported

if not __debug__:
    raise RuntimeError("Run without -O: exact assertions are required.")

cert = imported("augmented_exact_moment",ROOT / "exploration/indexed-cycle-target-10/independent_moment.py")


def add_form(left,right,coefficient=1):
    a,b=left; c,d=right
    if coefficient == -1: c,d=d,c
    assert not (a|b)&(c|d)
    return a|c,b|d


def scalar_check(c,k):
    roles=[(0,0)]*k['size']
    for source,slot in k['source'].items(): roles[slot]=(1<<source,0)
    for n,inputs,outputs in k['gates']:
        if c.args[n]:
            a,sign=c.args[n][0]
            assert roles[inputs[0]] == (c.pos[a],c.neg[a])
            value=roles[inputs[0]] if sign == 1 else tuple(reversed(roles[inputs[0]]))
            for slot,(a,sign) in zip(inputs[1:],c.args[n][1:]):
                assert roles[slot] == (c.pos[a],c.neg[a])
                value=add_form(value,roles[slot],sign)
            roles[inputs[0]]=value
        expected=(c.pos[n],c.neg[n])
        assert roles[inputs[0]] == expected
        for slot in outputs[1:]:
            roles[slot]=add_form(roles[slot],expected)
        assert all(roles[slot] == expected for slot in outputs)
    for piece,slot in k['pout'].items():
        n=c.pieces[piece][1]
        assert roles[slot] == (c.pos[n],c.neg[n])
    for i,n in c.retained:
        assert roles[k['rout'][i]] == (c.pos[n],c.neg[n])
        expected=[t-x[i] for x,_,t in c.roots]
        assert c.pos[n] == sum(1<<j for j,x in enumerate(expected) if x==1)
        assert c.neg[n] == sum(1<<j for j,x in enumerate(expected) if x==-1)
    denominator=lcm(2,c.h-3)
    by_target=[[] for _ in range(c.v)]
    for target,n,coefficient in c.pieces: by_target[target].append((n,coefficient))
    for target,(x,_,t) in enumerate(c.roots):
        row=[0]*c.v
        terms=by_target[target]+[(n,-Q(x[i],2)+Q(t,c.h-3)) for i,n in c.retained]
        for n,coefficient in terms:
            integer=denominator*coefficient
            assert integer.denominator == 1
            integer=int(integer)
            if not integer: continue
            for sign,mask in ((1,c.pos[n]),(-1,c.neg[n])):
                while mask:
                    bit=mask & -mask; mask ^= bit
                    row[bit.bit_length()-1] += sign*integer
        assert row == [denominator*int(j==target) for j in range(c.v)], target
    return {'clean_producer_checked_on_all_source_columns':True,
            'complete_shear_matrix_entries_checked':c.v*c.v,
            'scalar_common_denominator':denominator}


def histogram(c,k):
    h,v,R=c.h,c.v,k['size']; m=h*h; n=v*v
    current=[()]*R; local=Counter()
    def grow(slot,frame):
        previous=current[slot]
        assert c.valid_edge(previous,frame),(slot,previous,frame)
        difference=len(frame)-len(previous)
        if difference: local[difference]+=1
        current[slot]=frame
    for source,slot in k['source'].items(): grow(slot,c.label(source+1))
    for node,inputs,outputs in k['gates']:
        for slot in set(inputs+outputs): grow(slot,c.label(node))
    for i,node in c.retained:
        assert current[k['rout'][i]] == c.label(node)
        assert len(c.label(node)) == h-1
    for piece,slot in k['pout'].items(): grow(slot,c.target[c.pieces[piece][0]])
    for slot in range(R): grow(slot,c.full)
    assert sum(rank*count for rank,count in local.items()) == R*h
    for node in c.active:
        assert c.valid_edge((),c.label(node))
        assert c.valid_edge(c.label(node),c.full)
        for child,_ in c.args[node] or ():
            assert c.valid_edge(c.label(child),c.label(node))
    assert all(c.valid_edge((),target) for target in c.target)
    rows=Counter({r:2*v*count for r,count in local.items()})
    rows[m-h]+=2*v*R
    rows[(h-1)**2]+=2*n
    rows[h-1]+=4*n+2*v*h  # data paths and copied centres
    rows[1]+=n
    width=2*n+2*v*R
    loss=2*v*h*(h-1)
    mass=sum(r*count for r,count in rows.items())
    assert mass == width*m-n+loss
    assert max(rows) == m-h
    return dict(m=m,W=width,rows=rows,mass=mass,deficit=width*m-mass,
                N=n,R=R,loss=loss,local=local,
                distinct_frame_edges_checked=len(c.cache))


def elementary_word(c,k):
    v=c.v; offset=2*v
    mixer=[]
    for node,inputs,outputs in k['gates']:
        pivot=offset+inputs[0]
        if c.args[node]:
            coefficient=c.args[node][0][1]
            if coefficient != 1: mixer.append((pivot,None,Q(coefficient)))
            mixer.extend((pivot,offset+slot,Q(coefficient))
                         for slot,(_,coefficient) in zip(inputs[1:],c.args[node][1:]))
        mixer.extend((offset+slot,pivot,Q(1)) for slot in outputs[1:])
    inverse=[(a,b,-q if b is not None else 1/q) for a,b,q in reversed(mixer)]
    read=[]
    for i,_ in c.retained:
        for target,(x,_,t) in enumerate(c.roots):
            coefficient=-Q(x[i],2)+Q(t,c.h-3)
            if coefficient: read.append((v+target,offset+k['rout'][i],coefficient))
    read.extend((v+target,offset+k['pout'][j],coefficient)
                for j,(target,_,coefficient) in enumerate(c.pieces))
    injection=[(offset+slot,source,Q(1)) for source,slot in sorted(k['source'].items())]
    word=(mixer+[(a,b,-q) for a,b,q in read]+inverse+injection+
          mixer+read+inverse+[(a,b,-q) for a,b,q in injection])
    early=(len(mixer),len(mixer)+len(read))
    return word,early


def chronological_trace(c,k,word,profile):
    """Audit every scalar gate at its actual frame, including centre bridges.

    Stage two takes this same chronological trace with inverse-transposed
    scalar gates. It has identical ordinary incidences and reverses only
    each copied-centre bridge, not the original-role frame chronology.
    """
    v=c.v; offset=2*v; full=c.full
    mixer=[]; node_ops=[]
    for node,inputs,outputs in k['gates']:
        pivot=offset+inputs[0]; local=[]
        if c.args[node]:
            first=c.args[node][0][1]
            if first!=1: local.append((pivot,None,Q(first)))
            local.extend((pivot,offset+slot,Q(q)) for slot,(_,q) in zip(inputs[1:],c.args[node][1:]))
        local.extend((offset+slot,pivot,Q(1)) for slot in outputs[1:])
        mixer.extend(local);node_ops.append((node,local))
    inverse=[(a,b,-q if b is not None else 1/q) for a,b,q in reversed(mixer)]
    centres=[]
    for i,node in c.retained:
        rows=[]
        for target,(x,_,t) in enumerate(c.roots):
            q=-Q(x[i],2)+Q(t,c.h-3)
            if q:rows.append((v+target,offset+k['rout'][i],q))
        centres.append((node,rows))
    pieces=[(v+target,offset+k['pout'][j],q) for j,(target,_,q) in enumerate(c.pieces)]
    read=[operation for _,rows in centres for operation in rows]+pieces
    current=[c.label(j+1) for j in range(v)]+[()]*(v+k['size'])
    aux_hist=Counter();data_hist=Counter();bridge_hist=Counter();seen=[];digest=sha256()
    def grow(role,frame):
        old=current[role]
        assert c.valid_edge(old,frame),(role,old,frame)
        rank=len(frame)-len(old)
        if rank:(aux_hist if role>=offset else data_hist)[rank]+=1
        current[role]=frame
    def gate(operation,frame):
        a,b,q=operation
        grow(a,frame)
        if b is not None:grow(b,frame)
        seen.append(operation)
        digest.update(repr((a,b,str(q),frame)).encode())
    for operation in mixer:gate(operation,())
    for a,b,q in read:gate((a,b,-q),())
    for operation in inverse:gate(operation,())
    for source,slot in sorted(k['source'].items()):gate((offset+slot,source,Q(1)),c.label(source+1))
    for node,operations in node_ops:
        for operation in operations:gate(operation,c.label(node))
    for node,rows in centres:
        assert rows
        high=c.label(node)
        assert c.valid_edge((),high)
        assert all(current[a]==() and current[b]==high for a,b,_ in rows)
        bridge_hist[len(high)]+=1
        # Forward: copy high->0, scatter, erase. Transpose: gather at 0,
        # grow the fresh temporary 0->high, add there, erase. Same rank.
        seen.extend(rows)
        digest.update(repr(('centre_bridge',high,rows)).encode())
    for target in range(v):grow(v+target,c.target[target])
    for j,operation in enumerate(pieces):gate(operation,c.target[c.pieces[j][0]])
    for slot in range(k['size']):grow(offset+slot,full)
    for operation in inverse:gate(operation,full)
    for source in range(v):grow(source,full)
    for source,slot in sorted(k['source'].items()):gate((offset+slot,source,Q(-1)),full)
    assert seen==word
    assert aux_hist==profile['local']
    assert data_hist==Counter({c.h-1:2*v})
    assert bridge_hist==Counter({c.h-1:c.h})
    assert current[:v]==[full]*v
    assert current[v:offset]==c.target
    assert current[offset:]==[full]*k['size']
    return {'scalar_gates_with_chronological_frame_trace':len(seen),
            'trace_sha256':digest.hexdigest(),'ordinary_inverse_transpose_incidences_unchanged':True,
            'copied_centre_bridges_per_local_invocation':c.h,
            'local_data_growth_per_role':c.h-1,
            'stage_two_interstage_join_rank':(c.h-1)**2,
            'stage_two_rule':'Same order; G maps to G^-T; D_U=(a_perp tensor F)+(a tensor U); no bank swap'}


def envelope(word,size,transpose=False):
    values=[Q(1)]*size; denominators=[0]*size; peak=Q(1); denominator_bits=0
    for a,b,q in word:
        if transpose and b is not None: a,b=b,a
        if b is None:
            values[a] *= abs(q)
        else:
            values[a] += abs(q)*values[b]
            den=q.denominator
            assert den & (den-1) == 0
            denominators[a]=max(denominators[a],denominators[b]+den.bit_length()-1)
        peak=max(peak,values[a]); denominator_bits=max(denominator_bits,denominators[a])
    return peak,denominator_bits


def guard(c,k,word):
    size=2*c.v+k['size']
    forward,b1=envelope(word,size)
    transposed,b2=envelope(word,size,True)
    # Covers forming a dyadic numerator before its division.
    den=max(q.denominator for _,_,q in word)
    forward*=den; transposed*=den
    combined=2*forward*transposed
    bits=0
    while combined>2**bits: bits+=1
    charge=max(64,b1+b2,bits+2)
    coefficient=2*c.h+(charge+c.h-1)//c.h
    assert Q(coefficient,c.h)>=2+Q(charge,c.h*c.h)
    return dict(forward_scalar_envelope=str(forward),inverse_transpose_scalar_envelope=str(transposed),
                combined_two_stage_endpoint_envelope=str(combined),ceil_log2_combined=bits,
                scalar_denominator_bits_per_stage=[b1,b2],constant_prefix_charge=charge,
                recursive_guard_coefficient=coefficient,
                fits_retained_128_layer_coefficient_budget=coefficient+10<=128,
                scope='Conditional on the inherited whole-residual width and layer interfaces')


def dirty_basis_check(c,k,word,early):
    size=2*c.v+k['size']; peak,b=envelope(word,size)
    peak_t,b_t=envelope(word,size,True)
    scale=1<<max(b,b_t)
    maximum=scale*(max(peak,peak_t)+2)
    stride=maximum.numerator.bit_length()+2
    base=1<<stride
    assert maximum < Q(base,2)
    original=[scale << (stride*j) for j in range(size)]
    def run(transpose=False,omit_early=False):
        values=original.copy()
        for index,(a,b,q) in enumerate(word):
            if omit_early and early[0]<=index<early[1]: continue
            if transpose:
                if b is None: q=1/q
                else: a,b,q=b,a,-q
            if b is None:
                numerator=values[a]*q.numerator
                assert numerator%q.denominator == 0
                values[a]=numerator//q.denominator
            else:
                numerator=values[b]*q.numerator
                assert numerator%q.denominator == 0
                values[a]+=numerator//q.denominator
        return values
    forward=run()
    expected=original.copy()
    for j in range(c.v): expected[c.v+j]+=original[j]
    assert forward == expected
    transposed=run(True)
    expected=original.copy()
    for j in range(c.v): expected[j]-=original[c.v+j]
    assert transposed == expected
    assert run(omit_early=True) != forward
    return dict(basis_dimension=size,balanced_digit_stride_bits=stride,
                numerator_coefficient_bound=str(maximum),all_dirty_columns_checked_exactly=True,
                forward_shear_and_inverse_transpose=True,omitted_cancellation_negative_control=True)


def moment_check(p):
    if p['deficit']<=0:
        return {'positive_rank_deficit':False,'saving_claimed':None}
    def value(a):
        return sum(r*n*exp(a*log(p['m']/r)) for r,n in p['rows'].items())/(p['m']*p['W'])
    lo,hi=0.,.001
    for _ in range(50):
        middle=(lo+hi)/2
        if value(middle)<1: lo=middle
        else: hi=middle
    grid=10**10
    index=floor(lo*grid)
    while cert.moment(p,Q(index,grid))[1]>=1:index-=1
    while cert.moment(p,Q(index+1,grid))[0]<=1:index+=1
    accepted=Q(index,grid); rejected=Q(index+1,grid)
    lower,upper=cert.moment(p,accepted)
    next_lower,_=cert.moment(p,rejected)
    assert upper<1<next_lower
    return {'positive_rank_deficit':True,'saving_lower':str(accepted),
            'saving_upper':str(rejected),'acceptance_gap':str(1-upper),
            'next_grid_rejection_gap':str(next_lower-1)}


def audit(h,dirty=False):
    c=Producer(h); k=c.compile()
    scalar=scalar_check(c,k)
    profile=histogram(c,k)
    word,early=elementary_word(c,k)
    trace=chronological_trace(c,k,word,profile)
    result=dict(h=h,v=c.v,R=k['size'],active_nodes=len(c.active),signed_pieces=len(c.pieces),
                inlined_arguments_to_repair_alternating_residuals=c.inlined_arguments,
                scalar=scalar,elementary_scalar_gates=len(word),
                chronological_trace=trace,
                m=profile['m'],W=profile['W'],total_rank=profile['mass'],rank_deficit=profile['deficit'],
                local_history={str(r):n for r,n in sorted(profile['local'].items())},
                child_histogram={str(r):n for r,n in sorted(profile['rows'].items())},
                distinct_frame_edges_checked=profile['distinct_frame_edges_checked'],
                moment=moment_check(profile),guard=guard(c,k,word))
    if dirty: result['dirty_register_basis_audit']=dirty_basis_check(c,k,word,early)
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('h',type=int,nargs='*',default=[7,19])
    parser.add_argument('--dirty-small',action='store_true')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result={'status':'PASS','profiles':[audit(h,args.dirty_small and h==7) for h in args.h],
            'scope':'Finite complex word and moment only; all-size compiler/global multiplication integration not claimed'}
    value=json.dumps(result,indent=2)+'\n'
    if args.output: args.output.write_text(value)
    print(value,end='')


if __name__=='__main__': main()
