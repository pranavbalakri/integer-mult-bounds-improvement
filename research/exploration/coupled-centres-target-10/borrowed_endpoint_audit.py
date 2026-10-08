"""Does PR96-style endpoint merging apply to our complete borrowed word?
The audit retains every gate frame, including early cancellation and inverse
cleanup. It reads the pinned source's literal producer but imports no compiler.
"""
if not __debug__:
    raise RuntimeError('Run with assertions enabled, without -O')
from argparse import ArgumentParser
from pathlib import Path
from hashlib import sha256
import gzip,json,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'indexed-cycle-target-10'))
from independent_moment import check_source


def check(root,h):
    raw=gzip.decompress((root/f'certificates/indexed-cycle-word-{h}.json.gz').read_bytes())
    d=json.loads(raw);v,R=d['v'],d['R']
    source_slots=set(d['sources'].values());aux=set(range(R))-source_slots
    state=[1<<i for i in range(R)]
    touched=set()
    for a,b,fr in d['ops']:
        state[a]^=state[b];touched.update((a,b))
    output=[0]*v
    for a,b in d['scatter']:output[a-v]^=state[b-2*v]
    used=0
    for row in output:used|=row
    zero_columns=[s for s in sorted(aux) if not (used>>s&1)]
    assert not zero_columns
    assert aux<=touched
    # Every source pivot's coefficient remains the expected identity column.
    sources={int(i):s for i,s in d['sources'].items()}
    source_mask=sum(1<<s for s in source_slots)
    assert [row&source_mask for row in output]==[1<<sources[i] for i in range(v)]
    dirty_rows=[row&~source_mask for row in output]
    assert all(row for row in dirty_rows)
    # One deliberately omitted early coefficient leaves that dirt in the output.
    victim=next(i for i,row in enumerate(dirty_rows) if row)
    bit=dirty_rows[victim]&-dirty_rows[victim]
    assert (dirty_rows[victim]^bit)^dirty_rows[victim]==bit!=0
    return dict(h=h,old_auxiliary_roles=R,borrowed_auxiliary_roles=len(aux),
                producer_sha256=sha256(raw).hexdigest(),
                early_correction_nonzero_auxiliary_columns=len(aux),
                early_correction_zero_auxiliary_columns=0,
                all_auxiliaries_touched_by_inverse_cleanup=True,
                first_actual_auxiliary_gate_frame='D0 (early correction)',
                last_actual_auxiliary_gate_frame='D1 (inverse cleanup)',
                creation_gains_from_fixed_gate_endpoint_shift=0,
                retirement_gains_from_fixed_gate_endpoint_merging=0,
                omitted_early_coefficient_detected=True,
                scope='Our literal borrowed word with all gate frames fixed; no conclusion about a different circuit or endpoint contract.')

if __name__=='__main__':
    ap=ArgumentParser();ap.add_argument('--source',required=True,type=Path);args=ap.parse_args()
    pinned=check_source(args.source)
    result=dict(source_commit=pinned,axes=[check(args.source,h) for h in (23,25)],
                external_proposal='https://github.com/eumemic/integer-mult-bounds/blob/606d16d6dfc714d2a467190dc91b2f8dcab38d9c/research/merged-exterior/PROOF.md')
    Path(__file__).with_name('borrowed-endpoint-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
