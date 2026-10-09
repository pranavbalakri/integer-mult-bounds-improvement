"""Transfer the linear semantic guard to PR97's actual signed h24 word.

The credited pinned event ledger is replayed in a disposable output directory.
This independent audit then includes every scalar update and copied temporary
in exact prefix envelopes. Whole-residual recursion and global tape/analytic
interfaces remain explicit assumptions; this is not an all-size compiler proof.
"""
if not __debug__:raise RuntimeError('Run with assertions enabled, without -O')
from argparse import ArgumentParser
from array import array
from contextlib import redirect_stdout
from fractions import Fraction as Q
from hashlib import sha256
from io import StringIO
from pathlib import Path
import gzip,importlib.util,json,shutil,sys,tempfile

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec);sys.modules[name]=obj
    spec.loader.exec_module(obj);return obj


def literal(source):
    pkg=source/'research/deferred-signed'
    module=load('pr97_signed_h24_ledger',pkg/'round6_complex_literal_ledger.py')
    with tempfile.TemporaryDirectory(prefix='kappa-signed-guard-') as temp:
        root=Path(temp);lean=root/'swapnil-round7/lean';lean.mkdir(parents=True)
        shutil.copy2(pkg/'swapnil-round7/lean/round6-histograms.json',lean)
        module.HERE=root
        stdout=StringIO()
        with redirect_stdout(stdout):module.main()
        folder=root/'round6-complex-literal-ledger'
        raw=gzip.decompress((folder/'forward-events.i32.gz').read_bytes())
        events=array('i');events.frombytes(raw)
        frames=json.loads((folder/'frames-copies.json').read_text())
        receipt=json.loads((folder/'result.json').read_text())
        expected=json.loads((pkg/'round6-complex-literal-ledger/result.json').read_text())
        assert {k:v for k,v in receipt.items() if k!='elapsed'}=={k:v for k,v in expected.items() if k!='elapsed'}
        return events,frames,receipt,sha256(raw).hexdigest()


def envelope(events,frames,v,R,inverse=False):
    size=2*v+R;bank=lambda s:s+v if s<v else s-v if s<2*v else s
    if inverse:
        current=[None]*size
        for s in range(size):current[bank(s)]=~frames['final'][s]
    else:current=frames['initial'][:]
    val=[Q(1)]*size;den=[0]*size
    peak=Q(1);denpeak=0;updates=0;copies=0;numerator_peak=Q(1);coeffmax=Q(0)
    def add(t,bound,db,n,d):
        nonlocal peak,denpeak,updates,numerator_peak,coeffmax
        coeff=abs(Q(n,d));assert coeff.denominator in (1,2)
        # Charge the literal integer numerator BEFORE dividing, even when
        # n/d was written non-reduced (for example 2/2).
        product=abs(n)*bound;numerator_peak=max(numerator_peak,product)
        peak=max(peak,product)
        val[t]+=coeff*bound
        den[t]=max(den[t],db+coeff.denominator.bit_length()-1)
        peak=max(peak,val[t]);denpeak=max(denpeak,den[t],db)
        updates+=1;coeffmax=max(coeffmax,coeff)
    zero=frames['bases'].index([])
    order=range(len(events)-5,-1,-5) if inverse else range(0,len(events),5)
    for j in order:
        kind,a,b,n,d=events[j:j+5]
        if inverse:a=bank(a)
        if kind==0:
            before,after=(~n,~b) if inverse else (b,n)
            assert current[a]==before
            assert len(frames['bases'][n])-len(frames['bases'][b])==d>0
            current[a]=after
        elif kind==1:
            if inverse:b=bank(b);n=-n
            assert current[a]==current[b]
            add(a,val[b],den[b],n,d)
        else:
            assert kind==2
            cp=frames['copies'][n]
            assert (cp['source'],cp['frame'],cp['target_frame'])==((bank(a) if inverse else a),b,zero)
            assert current[a]==(~b if inverse else b)
            assert d==len(frames['bases'][b])
            # The temporary's scalar row duplicates THIS live source row;
            # its paid frame conversion does not change logical coefficients.
            temporary,temporary_den=val[a],den[a]
            peak=max(peak,temporary);denpeak=max(denpeak,temporary_den)
            for t,coefficient in (reversed(cp['terms']) if inverse else cp['terms']):
                if inverse:t=bank(t);coefficient=-coefficient
                assert current[t]==(~zero if inverse else zero)
                add(t,temporary,temporary_den,coefficient,2)
            copies+=1
    if inverse:
        assert all(current[bank(s)]==~frames['initial'][s] for s in range(size))
    else:assert current==frames['final']
    assert copies==24 and updates==374340
    return dict(inverse=inverse,largest_scalar_prefix=str(peak),
                largest_literal_numerator_product=str(numerator_peak),
                scalar_denominator_bits=denpeak,scalar_updates=updates,
                explicit_copied_temporaries=copies,largest_coefficient=str(coeffmax),
                every_ordinary_gate_has_equal_frames=True,
                every_copied_scatter_uses_its_actual_live_source=True)


def phase_controls():
    # Gaussian rationals, all represented exactly as pairs.
    add=lambda a,b:(a[0]+b[0],a[1]+b[1])
    mul=lambda a,b:(a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
    units=[(Q(1),Q(0)),(Q(0),Q(1)),(Q(-1),Q(0)),(Q(0),Q(-1))]
    z=units[0];zero=(Q(0),Q(0));a=(Q(1,2),Q(1,2));b=(Q(1,2),Q(-1,2))
    C=[[a,b],[b,a]];Ci=[[b,a],[a,b]];I=[[z,zero],[zero,z]]
    X=[[zero,z],[z,zero]];Z=[[z,zero],[zero,units[2]]]
    def matmul(A,B):return [[add(mul(A[i][0],B[0][j]),mul(A[i][1],B[1][j])) for j in range(2)] for i in range(2)]
    assert matmul(C,Ci)==I and matmul(X,C)==Ci
    ZCZ=matmul(matmul(Z,C),Z)
    assert [[mul(units[3],x) for x in row] for row in ZCZ]==Ci
    wrong_inverse=omitted_translation=0
    for u in range(4):
        for w in range(4):
            Xy=units[(w+2)%4];Yx=units[(w-2*u)%4];Yy=units[(w-u)%4]
            assert mul(Yx,units[2*u%4])==units[w]
            assert add(Yy,mul(units[-u%4],Xy))==zero
            wrong_inverse+=add(Yy,mul(units[u],Xy))!=zero
            omitted_translation+=Yx!=units[w]
    assert wrong_inverse==omitted_translation==8
    return dict(exact_one_axis_inverse_identity=True,
                inverse_equals_minus_i_Z_C_Z=True,inverse_equals_address_swap_C=True,
                endpoint_phase_pairs_checked=16,wrong_inverse_correction_rejections=8,
                omitted_pretranslation_rejections=8)


def run(source):
    baseaudit=load('pr97_pin_audit',HERE.parent/'deferred-interface-target-10/audit.py')
    pin_report=baseaudit.pins(source)
    foreign=source/'research/deferred-signed/swapnil-round7/independent/complex-twostage/producer.py'
    own=RESEARCH/'independent/complex-twostage/producer.py'
    assert foreign.read_bytes()==own.read_bytes()
    events,frames,ledger,event_digest=literal(source)
    forward=envelope(events,frames,2024,49208)
    backward=envelope(events,frames,2024,49208,True)
    guard=load('retained_semantic_guard',RESEARCH/'independent/guard-improvement/semantic_guard.py').certify(24)
    gf=Q(forward['largest_scalar_prefix']);gi=Q(backward['largest_scalar_prefix'])
    assert gf==Q(guard['forward_G']) and gi==Q(guard['inverse_G'])
    G=2*gf*gi;B=forward['scalar_denominator_bits']+backward['scalar_denominator_bits']
    assert G==Q(guard['G'])==3876535381862 and B==guard['B']==2
    assert G<2**48 and B<=8
    m=576;rho=Q(552,576);C=51
    assert Q(C)*(1-rho)>=2+Q(64,m) and C>=8 and C+10<=128
    assert max(map(int,ledger['full_histogram']))==552
    precision=Q(43,20);epsilon=Q(71,100)
    assert precision>epsilon and precision-3*epsilon==Q(1,50)>0
    identity=load('pr97_h24_support_identity',source/'research/deferred-signed/check_complex_identity.py')
    text=StringIO()
    with redirect_stdout(text):identity.main()
    assert '4096576 coefficients' in text.getvalue()
    return dict(status='PASS conditional linear-guard transfer to the actual PR97 signed h24 word',
        source_commit=baseaudit.HEAD,pin_report=pin_report,
        producer_sha256=sha256(own.read_bytes()).hexdigest(),producer_bytes_identical=True,
        reproduced_event_count=len(events)//5,reproduced_event_bytes_sha256=event_digest,
        h24_support_identity_coefficients=4096576,forward=forward,inverse_opposite=backward,
        scalar_prefix_product=str(G),scalar_denominator_bits=B,
        copied_endpoint_addition_factor=2,phase_controls=phase_controls(),
        paused_parent_charge='2e+64',largest_child_fraction=str(rho),
        recursive_guard_coefficient=C,layer_guard='128(d+1)',C1=1,
        precision=dict(exponent=str(precision),coefficient=768,digit_coefficient=128,
                       epsilon=str(epsilon),gaussian_power_gap=str(precision-3*epsilon),
                       gaussian_uniform_lower_bound='15/8',numerical_guard_valid_for_b_at_least=1),
        assumptions=['Inherited whole-residual signed complex compiler and exact child return contract',
                     'Depth-first one-active-child schedule and disjoint top-level axis pieces',
                     'Fixed-tape exact routing/address wrappers and copied-temporary implementation',
                     'Inherited normalization, prime selection, Gaussian error and recovery interfaces'],
        all_size_compiler_proved=False,new_multiplication_bound=False)


if __name__=='__main__':
    p=ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,default=HERE/'results.json');args=p.parse_args()
    result=run(args.source);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','scalar_prefix_product','scalar_denominator_bits','recursive_guard_coefficient','layer_guard','precision')},indent=2))
