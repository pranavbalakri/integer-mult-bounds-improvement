# Endpoint merging and the complete borrowed word

This is a scoped compatibility audit of the
[PR96 merged-exterior proposal](https://github.com/eumemic/integer-mult-bounds/blob/606d16d6dfc714d2a467190dc91b2f8dcab38d9c/research/merged-exterior/PROOF.md)
at its pinned source revision. It concerns our literal borrowed word, not the
validity of every possible implementation of that proposal.

The proposal changes auxiliary endpoint frames while retaining all gate frames.
An endpoint shift by minus the exterior projector `E` moves the exterior charge
to the first residual. If the first gate is a nonzero local frame `A_first`, this
can combine the exterior and local entrance into one complementary projector.
An analogous exit merge needs the last actual gate to occur before the local
full-frame cleanup.

Our complete borrowed word is different from its positive middle producer:

```
early D correction; M; retained scatter; side scatter; M^-1.
```

All early correction gates occur at the common base frame `D0`. All inverse
cleanup gates occur at the full stage frame `D1`. The last two phases cannot
be omitted merely because a positive producer role has a smaller first or last
support label.

The exact literal audit finds:

| Axis | Remaining auxiliary roles | Nonzero columns in early correction D |
| --- | ---: | ---: |
| 23 | 25063 | 25063 |
| 25 | 33048 | 33048 |

Therefore every remaining auxiliary must be read by the early correction.
There are no zero columns on which to skip that first incidence. Moving its
exterior to the entrance gives the residual `E` at that first `D0` gate; it
does not combine `E` with a later nonzero local entrance.

Every remaining auxiliary also participates in `M`, hence in the inverse word
performed at `D1`. Its last actual gate is consequently at `D1`. Combining the
previous local growth with the exterior would remove or cross a required
`D1` gate, violating the premise that gate frames are unchanged. Even retiring
each role immediately after its last inverse gate leaves that fact unchanged.

Thus the two fixed-gate endpoint modifications give no new merged child for
this literal word. Applying PR96's selected-edge numbers directly to our
borrowed histogram would be unsupported. A different dirty-cancellation or
restoration circuit may change the conclusion and requires a new audit.

The checker verifies the pinned source revision, replays the scalar middle
word, links every scatter to its role, tests every remaining column of `D`,
and confirms inverse-cleanup incidence on every auxiliary. It also detects a
missing early dirty-cancellation coefficient.

```sh
python3 research/exploration/coupled-centres-target-10/borrowed_endpoint_audit.py --source /path/to/pinned/PR84
```
