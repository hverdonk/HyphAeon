# Deer mouse globin sequences (added 2026-09-22)

*Peromyscus maniculatus* globin CDS downloaded from GenBank and merged into the
TOGA2 alignments, which contain no deer mouse globin sequence of their own.

## Downloads

| File | Source | Records |
| --- | --- | --- |
| `beta_globin_storz2009_cds.fasta` | GenBank `GQ139365`–`GQ139470`, deposited by Storz et al. 2009 | 106 CDS (53 `HBB-T1`, 53 `HBB-T2`) |
| `alpha_globin_cds.fasta` | GenBank search, *P. maniculatus* α-globin gene records | 455 CDS (448 `HBA`, plus `HBAT1-3`, `HBQ`, `HBZ` and two non-globin hits) |

Fetched as `rettype=fasta_cds_na`, so coding sequence only; the GenBank records
are genomic and contain introns.

## Isoform validation

Translating and reading the diagnostic sites **in mature numbering** (paper
site *n* = translated position *n*+1, initiator Met removed) recovers the
published isoforms exactly:

| Isoform | Sites | State | Count |
| --- | --- | --- | --- |
| αI | 50, 57, 60, 64, 71 | `PGAGS` | 226 |
| αII | 50, 57, 60, 64, 71 | `HGAGS` | 83 |
| αIII | 50, 57, 60, 64, 71 | `HAGDG` | 72 |
| βI (high-altitude) | 62, 72, 128, 135 | `GGAA` | 29 per paralog |
| βII (low-altitude) | 62, 72, 128, 135 | `ASSS` | 20 per paralog |

This matches Storz 2009's definitions and confirms both the download and the
reading frame.

## What was merged, and the caveats

One representative per gene, both the **high-altitude** allele class:

| Alignment | Accession | Allele | Named in alignment |
| --- | --- | --- | --- |
| `HBA1`, `HBA2` | `KJ726380.1` | αI (`PGAGS`) | `HLperManSon3` |
| `HBB` | `GQ139406.1` | βI `HBB-T1` (`GGAA`) | `HLperManSon3` |

Method: the hg38 row of each alignment was un-gapped to give an in-frame
reference, `cawlign -t codon -s BLOSUM62 -f refmap` mapped the deer mouse CDS
onto it (refmap discards insertions relative to the reference, so no new
columns are created), and the result was padded back with `-` at every column
where hg38 has a gap. Verified afterwards: every pre-existing sequence is
byte-identical, all rows retain equal length, and the merged row still reads
`PGAGS` / `GGAA` at the diagnostic sites in reference coordinates.

**Caveats, each a deliberate choice:**

1. **Named after an existing tree tip.** `aeon_core`'s loader drops any
   sequence absent from `speciesTree.nh`. The tree already contains two sister
   *P. maniculatus* assemblies, so the row is named `HLperManSon3`
   (*P. m. sonoriensis*) rather than editing the tree. The species assignment
   is correct; the specific assembly is a stand-in for branch placement, and
   the true accession is recorded above.
2. **α paralog assignment is arbitrary.** Human `HBA1` and `HBA2` encode
   identical proteins, so the deer mouse sequence is 85.2% identical to both
   and orthology cannot be resolved by similarity. The same sequence was added
   to both α alignments; treat them as alternatives, not independent evidence.
3. **One allele, not the polymorphism.** The studies' epistasis is between
   *within-species* allele classes (αI/αII/αIII, βI/βII). Only the
   high-altitude allele was merged, so the alignment carries the species, not
   the polymorphism. The other haplotypes are in the FASTA files above if a
   different design is wanted.
4. **`HBA2` trailing codon.** That alignment's human reference ends in a masked
   codon, so the merged row translates a stop at its final position.
