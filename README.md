# Exact certificates for the complete six- and seven-node regions

This package proves all layered facets for `(n,k,d)=(6,5,5)` and `(7,6,6)`
without a linear-code assumption. It also includes the smaller-node
certificates needed to make the finite proof chain self-contained.

## Verify

Python 3.10 or later is sufficient. No packages, solver, or network are needed.

```sh
python verify.py
```

Optional exact audit report:

```sh
python verify.py --json verification_results.json
```

## What is proved

For a symmetric exact-repair message polymatroid, let `F=H(E)`, `a=H(O_i)`,
and `b=H(S_ij)`. Define

```text
D_(n,r) = n(n-1)b + nr(r-1)a - r(r+1)F.
```

The child is the `(n-1)`-node system obtained by deleting a node and
conditioning on its entire incoming repair tuple. Each JSON proves

```text
scale * (p D_parent - q D_child) = sum_t weight_t I(X_t;Y_t|R_t),
```

with all weights positive integers. Here all ten files have `q=1`.
The note proves the contraction lemma, the cutset bases, and the induction.
Thus the files are complete finite proof chains, not isolated LP statuses.

| n | r | p | Positive CMI terms |
|---|---|---|---|
| 4 | 2 | 3 | 10 |
| 5 | 2 | 3 | 86 |
| 5 | 3 | 6 | 93 |
| 6 | 2 | 3 | 250 |
| 6 | 3 | 6 | 203 |
| 6 | 4 | 10 | 42 |
| 7 | 2 | 3 | 117 |
| 7 | 3 | 6 | 1144 |
| 7 | 4 | 20 | 501 |
| 7 | 5 | 15 | 213 |

For actual code parameters these imply

```text
n=6:  5 alpha >= B; 15 beta >= B;
      2 alpha + 5 beta >= B;
      6 alpha + 5 beta >= 2B;
      36 alpha + 15 beta >= 10B;
      4 alpha + beta >= B.

n=7:  6 alpha >= B; 21 beta >= B;
      7 alpha + 21 beta >= 3B;
      7 alpha + 7 beta >= 2B;
      42 alpha + 21 beta >= 10B;
      70 alpha + 21 beta >= 15B;
      5 alpha + beta >= B.
```

The endpoint bounds and the final mixed bound in each case follow from
elementary cutset reconstruction sets. Layered codes and space sharing attain
the resulting boundary. There is no general-n theorem in this package.

## Files

- `finite_region_note.tex`: self-contained short note, with references and the
  comparison to Chen and Tian, ISIT 2022.
- `finite_region_note.pdf`: compiled note.
- `certificates/`: ten unchanged certificate JSON files.
- `verify.py`: readable exact checker plus search-family audit and coverage checks.
- `verify_original.py`: unaltered independent checker from the earlier bundle.
- `manifest.json`: exact parameters and SHA-256 hashes of all ten JSON files.
- `verification_results.json`: successful exact audit produced for this release.
- `verification.txt`: output of the main checker and the original checker.

Compile the note with `pdflatex finite_region_note.tex` twice. Its bibliography
is embedded, so no BibTeX file is required.

## Certificate format

Node labels in JSON are zero-based. Bits enumerate the ordered pairs `(i,j)`,
`i != j`, lexicographically. A mask denotes the corresponding message tuple.
`R`, `X`, and `Y` give a CMI. Its four masks are ordered as

```text
R, R union X union Y, R union X, R union Y.
```

For each mask, repair closure adds every outgoing row whose whole incoming
column is present. The saved node permutation must send that closed set
exactly to the saved label. These identifications use only repair and symmetry.

The archived `claimed_sum` and `target_witnesses` use the negative-CMI sign
convention. The readable verifier changes both signs and independently
recomputes `scale*(p D_parent-q D_child)` from `n,r,p,q`. The harmless empty-set
coordinate is retained for literal coefficient matching; `H(empty)=0`.

## Search provenance and scope

The conditioning dictionary allows complete outgoing rows plus at most `n-2`
residual messages after closure and removal of full rows. Residual edges may
be disconnected or cyclic. The CMI arguments are each a message or a full row.
Only CMIs tight on both adjacent complete parity layers are retained. These
properties are rechecked on every term by `verify.py`.

This restriction is on conditioning sets, not all four entropy corners. It is
different from the later all-corner two-piece restrictions. In particular,
these 203- and 1144-term r=3 recursion certificates differ from the later 159-
and 423-term direct r=3 certificates. Neither list is claimed minimal.

The package reproduces exact verification and the support audit. It does not
rerun the original numerical discovery searches or provide their large LP
dictionaries. Failure in a restricted dictionary does not prove impossibility
in the full Shannon cone.

The ten certificates and original checker were copied unchanged from the
September 23, 2026 `recursion_handoff_bundle(1).zip`; its checksum is recorded
in `manifest.json`. The new readable checker, coverage audit, and literature
note were assembled on October 16, 2026. No new LP experiment was performed.
