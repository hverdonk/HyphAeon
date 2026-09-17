# epistaeon

Benchmark data and structure validation for the mutational-order / epistasis
white paper. No analysis engine is implemented here yet; this directory holds
the empirical ground truth, the coordinate set, and the scripts that validate
both. See [`data/README.md`](data/README.md) for the full audit.

---

> # ⚠️ READ BEFORE RUNNING ANY STRUCTURAL ANALYSIS
>
> ## PDB author numbering is **not** alignment numbering
>
> Three of the four benchmark scaffolds are mature proteins whose initiator
> methionine was cleaved. Their PDB author numbering is therefore shifted by
> **one residue** from the UniProt-canonical numbering that alignment columns
> and site tensors use.
>
> ```
> uniprot_position = pdb_author_position + offset
> ```
>
> | Structure | Chain | Protein | UniProt | **Offset** |
> | --- | --- | --- | --- | --- |
> | `2HHB` | A, C | α-globin | P69905 | **+1** |
> | `2HHB` | B, D | β-globin | P68871 | **+1** |
> | `1MBO` | A | myoglobin | P02185 | **+1** |
> | `1U19` | A, B | rhodopsin | P02699 | **0** |
>
> Skipping this shifts **every** hemoglobin and myoglobin contact call by one
> residue. It fails silently — no error, no gap, no mismatch — because the
> off-by-one lands on a real neighbouring residue. It is the single most
> likely way to get plausible-looking but wrong contact classifications.
>
> The steroid receptors are a separate scheme entirely: `2Q1H`, `2Q1V`,
> `2Q3Y`, `3RY9` and `3GN8` use **local LBD numbering** (`-2` … `247`), not
> UniProt. To reach human GR (NR3C1) numbering, add **+531**. Positions
> `-2`, `-1` and `0` are expression-tag residues, not biology. `1U19` carries
> an `ACE` acetyl cap at author position `0` that must be skipped when
> iterating residues.
>
> **Do not hardcode these.** Load
> [`data/numbering_offsets.json`](data/numbering_offsets.json), the
> machine-readable source of truth, and regenerate it with
> `python scripts/numbering.py` if the coordinate set ever changes.

---

## Layout

| Path | Contents |
| --- | --- |
| `data/coordinates/` | Eight mmCIF files: four case scaffolds + the Thornton ancestral receptor series |
| `data/experimental/` | Published source papers and supplements for the four cases |
| `data/uniprot/` | Canonical reference sequences, so the numbering check runs offline |
| `data/numbering_offsets.json` | **Verified numbering crosswalk — load this, don't hardcode** |
| `data/README.md` | Full structure audit and contact-mapping fitness report |
| `scripts/` | Reproduces every number in `data/README.md` |

## Scripts

Each runs standalone with no arguments (except `validate_cif.py`, which takes
mmCIF paths). They need `biopython`, `numpy` and `scipy`.

| Script | Purpose |
| --- | --- |
| `numbering.py` | Verifies the offsets in the box above against UniProt |
| `validate_cif.py` | Entity, chain, organism, ligand and gap inventory |
| `probe.py` | Residue identity at the positions the white paper names |
| `align_ladder.py` | AncCR → AncGR1 → AncGR2 substitution ladder; derives the +531 crosswalk |
| `contacts.py` | Reproduces the §4 Table 2 distance and contact-count claims |
| `storz.py` | Maps the real Storz 2009 deer-mouse sites onto 2HHB |
| `mapfitness.py` | Distance-class composition, copy stability, interface ambiguity, cofactor bridging |
| `bands.py` | Band coverage gaps and statistical power for the §5 benchmark |
| `interface_states.py` | Blocker 2 evidence: interface disjointness and T/R state dependence |
| `_fetch.py` | Resolves a structure locally, else caches a download from RCSB |

## Open decisions

1. **Distance classes** — a third intermediate class and a per-protein upper
   bound will be implemented at analysis time (the ≤5 Å / 10–25 Å scheme
   leaves 29–32% of globin pairs and 54% of rhodopsin pairs unclassified).
2. **Hemoglobin interface convention** — open. The α1β1 and α1β2 contact sets
   are *disjoint*, and α1β2 turns over half its contacts between T and R.
   Recommendation and literature conventions in `data/README.md`.
3. **Cofactor-bridged coupling** — open. 41–49% of cofactor-lining pairs sit
   >10 Å apart and would be mis-called allosteric. Recommendation and
   literature conventions in `data/README.md`.
