# Equivariant monad bundles on CICY free quotients — a defect-audited scan

[![DOI](https://zenodo.org/badge/1346059966.svg)](https://doi.org/10.5281/zenodo.22095895)

**Status: closed at tag `v1-monades`.** This repository is complete and is not
under active development. What it establishes is below; what it does *not*
establish is stated just as explicitly, and so is what it never explored.

The full working record — 220 kB, in French — is
[`CY_Landscape_Explorer_PROJET.md`](CY_Landscape_Explorer_PROJET.md). This
README is the English summary.

---

## What this is

A systematic scan of two-term monad bundles

```
0 → V → B = ⊕ᵢ O(bᵢ) → C = ⊕ₐ O(cₐ) → 0
```

over complete intersection Calabi–Yau threefolds admitting a free quotient,
with the aim of finding three-generation heterotic models. Every stage —
cohomology, stability, surjectivity, equivariance under Γ — is computed rather
than assumed, and every stage **declares what it did not determine**.

The scan of record covers **194 CICYs** with `|χ| = 3|Γ|` (Braun's
classification of free quotients), **14 943 candidate identities**, and
**505 601 verdict lines**, computed in a single state of the code.

## Result

| | |
|---|---|
| verdict lines | 505 601 |
| survive all tests | **34 885** — all SO(10), rank 4, `n_gen(X/Γ) = 3` |
| eliminated: `V = ker f` is not a bundle | **28 006** — base locus exhibited, re-substituted, shown to lie **on Y** |
| undetermined | **34 693** — all in one stratum, `rank_C = 1 / rank_V = 4` |
| discarded before evaluation | 36 898 — outside the model's domain, index filter, etc. |
| distinct (B, C) | 2 857, on **91 CICYs**, 691 orbits under Aut(configuration) |

**Two candidates are fully verified**, `#6890` and `#6947`: three generations,
`h⁰(V) = 0` on the Γ-equivariant subspace, complete Hoppe criterion,
surjectivity of `f` certified, anomaly condition checked.

**And a limitation of principle closes the line.** A Wilson line preserves the
rank of the gauge group; SO(10) has rank 5 and the Standard Model rank 4. Both
candidates therefore cap at Pati–Salam or flipped SU(5). The route that could
go further — `|Γ| ≥ 4` — is closed and measured: 574 pairs, 544 killed by
equivariant `h⁰(V)`, 28 with no equivariant `f`, **0 survivors**.

## Position with respect to arXiv:0911.1569

Anderson, Gray, He and Lukas (JHEP 02(2010)054) scanned **7 118 positive**
monads over 4 515 favourable CICYs and ruled the class out on phenomenological
grounds. Their positivity is **strict** — `bᵢʳ > 0` and `cₐʳ > 0` for all
indices, their eq. (2.11).

The two candidates here carry **20 zeros in B and 2 in C** each. They are
semi-positive, and thus lie **outside** those 7 118. This scan imposes only
their eq. (2.6), `cₐ ≥ bᵢ` componentwise.

Three further differences, all checkable:

- their quotient coverage is **five manifolds**; this scan uses the **194**
  free-quotient CICYs of Braun's classification, published four months after
  their paper;
- their footnote 2 sets aside the case-by-case analysis required when (2.6) is
  relaxed. That analysis is done here, on the **constrained** `f` rather than a
  generic one — which matters, since the equivariant `f` is a special point and
  nothing guarantees it avoids the locus where `f` drops rank;
- their three SO(10) models lie on none of the manifolds used here.

The authors themselves note that positivity "is not necessary" for `h¹¹ > 1`,
that semi-positive monads "may well be stable", and that positive monads are
"likely a small sub-set of all stable monad bundles". They stated in 2009 that
a systematic scan of non-positive monads was already underway. **We have not
found such a scan published since, but absence of a search result is not proof
of absence** — this point needs checking by someone with the field's citation
habits, and it is the first thing a reader should challenge.

## What is *not* established

Stated up front rather than in a footnote:

- **All verdicts are computed mod p**, with the reserve stated explicitly
  everywhere. The reserve runs in the favourable direction for eliminations —
  a rank computed mod p can only drop, so a vanishing `h⁰` is conclusive and a
  non-vanishing one is not — and it was measured rather than assumed: 32
  evaluations across four primes for the ℤ₄ verdict, 24/24 for each base-locus
  argument. It is not a characteristic-zero proof.
- **There has been no independent verification.** Everything comes from one
  codebase. The regression suite (47 tests, external references only) is the
  answer offered to that objection, not a substitute for it.
- **34 693 verdict lines remain undetermined**, all in the stratum
  `rank_C = 1 / rank_V = 4`. They are counted and labelled, never silently
  treated as absent.

## What was never explored

Three routes, named — an absence is not a result of absence:

1. **Rank-5 extensions.** Never generated: 1.1·10⁸ ordered tuples at m = 3,
   beyond any enumeration ceiling. This is the one place where something
   unexpected could still live in this domain.
2. **The Čech correction block**, which decides `#6947`: its charge `c₁ + b₄`
   gives `dim(S/I) = 84` against `χ = 76`, and the 8 missing units are exactly
   the Čech classes left unbuilt.
3. **The stratum `rank_C = 1 / rank_V = 4`** — 686 λ-lines, 34 693 in the
   file. Its shape is measured and favourable (a single repeated configuration
   on P¹×P¹×P¹, 91 CICYs) but `f` has five components there, not four: the
   dimension count has to be redone, not copied.

## The part that may be the most useful

Nine implementation defects were found over the course of this work. **None was
found by the code itself** — all by an external reference: Riemann–Roch, Serre
duality, a value known in advance, the definition of a reduction. Two of them
had produced a publishable number that did not exist.

They share one shape, documented eighteen times in §8 of the working record:
**a condition becomes empty, or universally true, and nothing says so.** The
result keeps coming out, with the appearance of a selection. Examples:

- a generator drawing **10 samples from a family of 2 201**, read as "these
  bundles do not exist";
- a results file not carrying the version of the code that wrote it — **4 049
  candidates wrongly discarded**, reinstated by resumption after resumption;
- an empty `groupes_utiles` silently falling back to *all* groups: **95.5 % of
  pairs computed off-target**, some marked as survivors;
- a word — `indetermine` — covering "not computed", "computed without
  concluding" and "the object is not a bundle", three situations requiring
  three opposite actions.

The discipline that catches them is in §8: count both sides, never substitute
zero for what was not computed, persist discarded lines with their reason —
and **break the filter in both directions**, so that a module accepting
everything and a module rejecting everything each fail a different assertion.
Nine of the 47 tests freeze two opposite verdicts for exactly that purpose.

## Repository contents

```
cy_landscape/core/     the mathematical core: cohomology, sections, monads,
                       covariant rings, Hoppe, exact pruning
equivariance_f.py      the full chain at polynomial level, and the sweep
lieu_de_base_rv3.py    exact base locus, stratum rank_C = 1 / rank_V = 3
lieu_de_base_rc2.py    exact base locus at rank_C = 2 (2×2 minors)
rencontre_F_Y.py       is the witness actually ON Y? (mandatory guard)
empreinte_code.py      the code version, written into every result line
tests_regression.py    47 tests — run before anything else
ancres_port.py         does the sweep reproduce the established verdicts?
compter_strates.py     cost of a stage, measured before it is paid
comparer_scans.py      two sweeps compared, in BOTH directions
```

Start with `python tests_regression.py` (≈ 5 min). Usual commands are in §9 of
the working record.

## Data provenance

- `cicylist.txt` — the CICY list, from the [CICY page at Oxford](http://www-thphys.physics.ox.ac.uk/projects/CalabiYau/cicylist/).
- `cicyquotients.m` — free quotients, from V. Braun's classification
  ([arXiv:1003.3235](https://arxiv.org/abs/1003.3235)).

Both are third-party data, redistributed here for reproducibility only, and
credited to their authors. If you are one of them and would rather they were
not mirrored, open an issue and they will be replaced by a download script.

## The two sweep archives

The two reference sweeps are too large for git and are archived separately,
with their checksums recorded here and in §0 of the working record — an archive
on a disk is a file that *claims* to be the right one:

```
scan_wilson5.zip  sha256  4C1D3C563CC49EF57C416D566CE48752FA925802F53011026591F9F91DB70C9D
scan_wilson6.zip  sha256  250A52A9A507CA2E6FE6E1F5A08A5533D9221197EF15C45E82EA4800F97B88E9
```

`scan_wilson6` is the file of record; `scan_wilson5` is kept as the comparison
base that established that **152 survivors were gained and none lost**.
Deposited at [https://doi.org/10.5281/zenodo.22099905](https://doi.org/10.5281/zenodo.22099905).

## Successor

The Standard Model is out of reach in this domain, for the rank reason above.
The route that reaches it uses a **rank-5 sum of line bundles** with structure
group `S(U(1)⁵) ⊂ SU(5)` — SU(5) has rank 4, so a Wilson line preserves the
rank all the way down. That is a different mathematical core, hence a different
repository.

## Licence and citation

Code: MIT ([`LICENSE`](LICENSE)). Documentation and produced data: CC BY 4.0
([`LICENSE-docs`](LICENSE-docs)). Third-party data keep their own terms.

Please cite via [`CITATION.cff`](CITATION.cff) or the DOI badge above.

## A note on how this was made

The implementation was carried out with substantial AI assistance. The
verification protocol, the external references used, and every defect found —
including how it was found and what it had cost — are documented in the working
record. The results stand or fall on that record, not on who or what typed the
code.
