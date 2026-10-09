# Auditing deferred interfaces and trading saving for smaller row stock

**The target kappa > 2^-10 remains unmet.** This package independently checks
the finite numerical inputs of a newer, separately credited construction and
accounts for modifications to its deferred entrance frames. It does not
replace the main multiplication certificate or prove the inherited analytic
and fixed-tape interfaces.

## Source and credit

The source is [Zhihao Chen's PR97](https://github.com/CrocSwap/integer-mult-bounds/pull/97),
frozen at `f5f9c56e637463cac1e300d1589ccf42838f688a`. Its bit and complex
producers are [Swapnil Jain's round-seven/round-six work](https://github.com/Swapnil-jain/integer-mult-kappa/tree/741e7aa078392553815df7926ee17ac5e25a8c38),
frozen at `741e7aa078392553815df7926ee17ac5e25a8c38`. Jain supplies the
deferred readout schedule, lifted frames, late and source copies, common
flag/staircase family, and paired complex producer. Chen supplies the
explicit reflected schedules, fan-order clarification, gauged auxiliary
exits, signed complex correction and balanced assembly used here. Their
preserved source notices credit the earlier contributors and AI assistance.
This package and the adjacent frame modifications were written with OpenAI
Codex assistance. Reusing this package is not a claim to have discovered
the original producers, or to have surpassed their published saving.

## Reproduction finding and independent checks

The pinned PR97 `SOURCE.json` lists 138 files, of which 134 are present and
match exactly. Its four historical validation logs are absent:

```
research/deferred-signed/validation-prior/frames-replay.log
research/deferred-signed/validation-prior/lifted-replay.log
research/deferred-signed/validation-prior/stair-replay.log
research/deferred-signed/validation-prior/stair_control-replay.log
```

All 102 files in its imported-source manifest match. Consequently PR97's
top-level `verify.py --replay-own` stops with `FileNotFoundError` at this
commit. We do **not** report that command as passing, edit its manifests,
or manufacture the missing historical logs.

Our `audit.py` pins both complete manifest bytes, reports the four missing
logs, verifies every other listed file, and independently encloses both
moments using the repository's rational logarithm/exponential checker.
Its optional replay runs the unmodified source scripts in a disposable copy:
both explicit event/frame ledgers, tensor endpoint controls, signed complex
controls, the complete h24 support identity, and balanced assembly. All six
replays pass, and all five regenerated JSON outputs agree with their source
versions after removing only elapsed times. SymPy 1.14.0 is needed by the
symbolic endpoint control; the numerical audit itself uses the standard library.

We also reran the five native source commands: `check_word.py 23`,
`check_frames.py 23`, `check_lifted.py 23`, and `check_stair.py 23` with and
without `--control`. All exited successfully, including the expected
`certified=False` control. Their complete new outputs and hashes are in
`native-replay-results.json`. These are newly generated receipts, not the
absent historical logs. The unchanged geometry audit checks all rational
containments, and the lifted-frame check covers 15,138 distinct high-rank
steps under both common axis bases.

These checks support the published finite inputs. The bit saving
`31987/500000000` and complex saving `36926111/500000000000` both have
strictly positive exact moment gaps. Both moments **reject** `1/1024`.
The source's reported conditional multiplication saving is
`63965813/10^12`, about `0.000063965813`; it is not our new result.

## Smaller entrance frames

`tradeoffs.py` reconstructs the complete bit histogram from the pinned
literal schedule. It changes each entrance dimension and recomputes the
first internal edges, gauged exterior corners, and every target-data chain.
Source-data paths, copied centres, data connectors, scalar word and roles
are unchanged. Every profile has

```
m = 529, W = 108516254, total rank = 57403754177.
```

The identity `total rank = m W - N + 2 v h(h-1)` is checked for every row.
Each displayed saving is an accepted rational lower endpoint, with the root
also proved below that number plus `10^-10`. It is a **finite bit-moment
saving**, not automatically the final multiplication kappa.

| Entrance choice | Bit saving, accepted lower endpoint | Largest child | Least halving depth | Product-row coefficient | Rounded stock degree |
|---|---:|---:|---:|---:|---:|
| Original | 0.0000639830 | 527 | 183 | 5417 | 12000 |
| Cap dimensions at 21 | 0.0000637083 | 525 | 92 | 2960 | 7000 |
| Cap dimensions at 20 | 0.0000631911 | 523 | 61 | 2123 | 5000 |
| Cap dimensions at 19 | 0.0000624868 | 521 | 46 | 1718 | 4000 |
| Common codimension-one cut | 0.0000622931 | 525 | 92 | 2960 | 7000 |
| Common codimension-two cut | 0.0000607136 | 523 | 61 | 2123 | 5000 |
| Common codimension-three cut | 0.0000592269 | 521 | 46 | 1718 | 4000 |
| Common codimension-four cut | 0.0000578258 | 519 | 37 | 1475 | 4000 |

The [frame certificate](../deferred-truncation-target-10/README.md) realizes
the first cap and common cuts 1 through 4 over the rationals, using the same
fixed axis bases throughout. The [coordinated-cap certificate](../deferred-cap-target-10/README.md)
realizes caps 20 and 19, with all required predecessors and shared targets
checked. Cap 19 reuses 1,660 existing frames and introduces one new frame.
The dimension calculation alone does not establish these geometric facts.
Rows for larger common cuts in the machine-readable result are arithmetic
explorations only unless accompanied by a separate frame certificate. The
zero-entrance endpoint is also certified: the explicit common cut
`H=span(e_0)` intersects every original deferred frame trivially. Its bit
saving is `0.0000504973`, halving depth 8, and h24 stock degree 2000. With
that complex side it is weaker on both saving and stock degree than the
earlier borrowed PR84/h18 choice, so it is not presented as an improvement.

The halving depth is the least integer d for which `m^d > 2 maxchild^d`.
It describes a sufficient number of recursive levels to halve every child
dimension. With the unchanged h24 complex side, whose depth is 17, the
retained row bound has coefficient

```
bit_length(W_bit) d_bit + bit_length(W_complex) d_complex.
```

The rounded stock degree is the source's conservative multiple-of-1000
choice exceeding `51/25` times this coefficient. These are sufficient
proof bounds, not measured memory use, running times or crossover sizes.
The role count W does not fall. Nor do these choices dominate every older
construction: our previous borrowed PR84/h18 candidate has lower saving
`0.00005487` but rounded stock degree 1200, below all displayed h24 choices.

## Conditional assembly and retained lower precision

The [actual signed-word guard audit](../deferred-guard-target-10/README.md)
replays all 535,744 h24 events, including copied temporary streams. It
independently obtains scalar envelopes 1,006,577 and 1,925,603, and transfers
the existing `128(d+1)` linear guard to this precise word. Its proof includes
the inverse rank-one correction and address translation. The inherited
whole-residual compiler, tape layout and analytic interfaces remain assumed.

`assembly_options.py` recomputes every selected bit profile and the complex
moment, binds the corresponding finite frame receipts, and checks every
inequality of our existing enlarged-precision macro stack with the actual
h24 complex saving substituted. This is a distinct parameter stack from
PR97's balanced assembly; their numerical margins must not be mixed.

| Entrance choice | Larger-saving candidate kappa | Lower-precision candidate kappa | Rounded stock degree |
|---|---:|---:|---:|
| Original, credited word | 0.00006396 | 0.00004540 | 12000 |
| Cap 21 | 0.00006370 | 0.00004521 | 7000 |
| Cap 20 | 0.00006318 | 0.00004485 | 5000 |
| Cap 19 | 0.00006248 | 0.00004435 | 4000 |

The larger-saving column uses `epsilon=0.9999` and `p=768 b^3`. The
lower-precision column uses `epsilon=0.71` and
`p=768 ceil(b^(43/20))`. These are **alternative parameter choices**; the
larger saving cannot be paired with the smaller precision from the other
column. Here `b=ceil(log2 n)` as in the retained assembly.

All larger-saving candidates exceed `2^-14`, and all are far below `2^-10`.
Cap 19 loses about 2.3% of the original finite bit saving while reducing the
least halving depth from 183 to 46 and the specified rounded stock degree
from 12000 to 4000. These are improvements to sufficient theoretical bounds
relative to this credited deferred construction. The earlier PR84/h18
candidate remains another tradeoff with a smaller stock degree.

## Reproduce

Use a checkout or archive of the exact PR97 commit, preserving its tree.
The command argument below is the repository root, not its imported subfolder.

```
python3 audit.py --source /path/to/pr97 --replay --output audit-results.json
python3 tradeoffs.py --source /path/to/pr97 --output tradeoff-results.json
python3 assembly_options.py --source /path/to/pr97 --check-guard --output assembly-options.json
```

The moment checker is independent of the source's arithmetic implementation.
`tradeoffs.py` deliberately imports the credited schedule/accounting module;
its original histogram must exactly equal the separately pinned literal-ledger
histogram before any changed profile is considered. For realized frames,
also run the adjacent frame checkers; a saved JSON receipt alone is not a
geometric proof. All new entry points fail closed under Python `-O`.
