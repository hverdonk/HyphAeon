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
> Two of the four benchmark structures, 2HHB and 1MBO, are mature proteins
> whose initiator methionine was cleaved. Their PDB author numbering is therefore shifted by
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
> | `7PRX` | A | glucocorticoid receptor | P04150 | **0** |
>
> Skipping this shifts **every** hemoglobin and myoglobin contact call by one
> residue. It fails silently — no error, no gap, no mismatch — because the
> off-by-one lands on a real neighbouring residue. It is the single most
> likely way to get plausible-looking but wrong contact classifications.
>
> **Case 3 has a second trap.** The white paper names receptor residues in
> *ancestral* (AncCR) numbering — Tyr27, Ser106, Leu111. 7PRX uses human GR
> numbering, so add **+531** (Ser106 → GR 637, Leu111 → GR 642). The backup
> ancestral structures in `data/ancestral_coordinates/` use local LBD
> numbering (`-2` … `247`), and in 3GN8 the +531 rule breaks after about
> position 212.
>
> `1U19` carries an `ACE` acetyl cap at author position `0` that must be
> skipped when iterating residues.
>
> **Alignment side, too.** The TOGA2 model inputs in
> `data/alignments/toga2/` report *alignment columns*, and the NR3C1 reference
> is the GRγ isoform (+1 Arg at 452), so human positions after 451 are one
> higher than 7PRX. See [`data/alignments/toga2/README.md`](data/alignments/toga2/README.md).
>
> **Do not hardcode these.** Load
> [`data/numbering_offsets.json`](data/numbering_offsets.json), the
> machine-readable source of truth, and regenerate it with
> `python scripts/numbering.py` if the coordinate set ever changes.

---

## Layout

| Path | Contents |
| --- | --- |
| `data/coordinates/` | The active set: `2HHB`, `1MBO`, `7PRX`, `1U19` |
| `data/alternate_coordinates/` | Backup: `4LSJ` (human GR; domain-swapped, superseded by 7PRX) |
| `data/ancestral_coordinates/` | Backup: Thornton ancestral receptor series (`2Q1H`, `2Q1V`, `2Q3Y`, `3RY9`, `3GN8`) |
| `data/experimental/` | Published source papers and supplements for the four cases |
| `data/uniprot/` | Canonical reference sequences, so the numbering check runs offline |
| `data/alignments/toga2/` | TOGA2 codon alignments + species tree — the model inputs |
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
| `interface_states.py` | Hemoglobin interface disjointness and T/R state dependence |
| `_fetch.py` | Resolves a structure locally, else caches a download from RCSB |

## Contact classification

Binary: a detected pair is in direct contact if any heavy atoms are within
5.0 Å, otherwise not. Apply the numbering offset first, exclude pairs with
|i − j| < 5, and count hemoglobin α–β pairs as contacts at either interface.
Baselines and details in [`data/README.md`](data/README.md).
