# Empirical data and coordinate inventory

These are published source files collected for reviewing
`white_paper_mutational_order_timing_epistasis.pdf`. Files are grouped by
empirical case, not by a claim that they constitute complete mutational-path
ground truth. No `axomeme_site_tensors.h5` was found under the parent directory
when this inventory was prepared (2026-09-15).

## Experimental sources

| Case | Local files | Published source and scope |
| --- | --- | --- |
| Hemoglobin | `experimental/hemoglobin_2009_storz.pdf` | Storz et al., *PNAS* 2009, [doi:10.1073/pnas.0905224106](https://doi.org/10.1073/pnas.0905224106). High-altitude hemoglobin functional data; not a complete ancestral genotype panel. |
| Hemoglobin | `experimental/hemoglobin_2015_andean_waterfowl.pdf`, `hemoglobin_2015_andean_waterfowl_oxygen_data.docx`, `hemoglobin_2015_andean_waterfowl_isoform_data.docx`, `hemoglobin_2015_andean_waterfowl_alignment.pdf` | Natarajan et al., *PLOS Genetics* 2015, [doi:10.1371/journal.pgen.1005681](https://doi.org/10.1371/journal.pgen.1005681), [S2 oxygen-affinity table](https://doi.org/10.1371/journal.pgen.1005681.s013), [S1 isoform table](https://doi.org/10.1371/journal.pgen.1005681.s012), [S2 alignment figure](https://doi.org/10.1371/journal.pgen.1005681.s002). The oxygen-affinity data use KCl and IHP conditions; they are not 2,3-DPG assays. |
| Myoglobin | `experimental/myoglobin_2013_article.pdf`, `myoglobin_2013_supplement.pdf` | Mirceta et al., *Science* 2013, [doi:10.1126/science.1234192](https://doi.org/10.1126/science.1234192). Supplement includes species traits, myoglobin concentration/charge data, and ancestral inference; it does not supply a complete measured recombinant ancestral-mutation panel. |
| Steroid receptor | `experimental/steroid_receptor_2007_ortlund.pdf`, `steroid_receptor_2009_bridgham.pdf` | Ortlund et al., *Science* 2007, [doi:10.1126/science.1142819](https://doi.org/10.1126/science.1142819); Bridgham et al., *Nature* 2009, [doi:10.1038/nature08249](https://doi.org/10.1038/nature08249). Reconstructed receptor mutagenesis, ligand specificity, and epistasis. |
| Steroid receptor | `experimental/steroid_receptor_2011_plos.pdf`, `steroid_receptor_2011_ec50_data.doc`, `steroid_receptor_2011_sequence_accessions.doc`, `steroid_receptor_2011_ancestral_sequences.doc`, `steroid_receptor_2011_stability_data.doc` | Carroll et al., *PLOS Genetics* 2011, [doi:10.1371/journal.pgen.1002117](https://doi.org/10.1371/journal.pgen.1002117), [S1 EC50 data](https://doi.org/10.1371/journal.pgen.1002117.s003), [S2 accession list](https://doi.org/10.1371/journal.pgen.1002117.s004), [S3 ancestral sequences](https://doi.org/10.1371/journal.pgen.1002117.s005), [S5 stability calculations](https://doi.org/10.1371/journal.pgen.1002117.s007). The stability table is calculated, not measured thermal unfolding. |
| Rhodopsin | `experimental/rhodopsin_2017_cetacean_supplement.pdf` | Dungan & Chang, *Proceedings of the Royal Society B* 2017, [doi:10.1098/rspb.2016.2743](https://doi.org/10.1098/rspb.2016.2743), [published supplement](https://doi.org/10.6084/m9.figshare.4653931.v1). Includes mutant spectral and retinal-release results; retinal release after photoactivation is not a direct dark-noise rate. The [Dryad raw archive](https://doi.org/10.5061/dryad.5k0s6) was located but returned HTTP 401 on download. |

The *Science* 2013 hemoglobin combinatorial experiment of Natarajan et al.,
[doi:10.1126/science.1236862](https://doi.org/10.1126/science.1236862),
assayed eight combinations of three allele blocks and additional single/double
variants. Its publication is available online, but its supporting genotype
data were not obtainable as a local file from the public endpoints tested.
Likewise, the 2008 ancestral-rhodopsin PNAS study of Yokoyama et al.,
[doi:10.1073/pnas.0802426105](https://doi.org/10.1073/pnas.0802426105),
was located online, but its publisher supplement could not be downloaded.

## Coordinate sources

All coordinate files are RCSB/wwPDB mmCIF downloads from
`https://files.rcsb.org/download/<PDB-ID>.cif`. Resolution values were checked
against each file's `_refine.ls_d_res_high` field.

| File | Structure | Resolution | Role |
| --- | --- | --- | --- |
| `coordinates/2HHB.cif` | Human deoxyhemoglobin | 1.74 Å | Section 4 accession |
| `coordinates/1MBO.cif` | Sperm-whale oxymyoglobin | 1.60 Å | Section 4 accession |
| `coordinates/1U19.cif` | Bovine rhodopsin | 2.20 Å | Section 4 accession |
| `coordinates/4LS6.cif` | *Bacillus subtilis* FabF I108F | 1.56 Å | Section 4 accession, but **not** a steroid receptor |
| `coordinates/2Q1H.cif` | Ancestral corticoid receptor–aldosterone | 1.90 Å | Relevant receptor replacement |
| `coordinates/3RY9.cif` | AncGR1–DOC | 1.95 Å | Relevant receptor replacement |
| `coordinates/3GN8.cif` | AncGR2–dexamethasone | 2.50 Å | Relevant receptor replacement |
| `coordinates/4P6X.cif` | Human GR–cortisol | 2.50 Å | Relevant receptor replacement |

These structures do not by themselves establish a tensor-to-coordinate
crosswalk. It requires the tensor's reference sequence and site-index
metadata, a sequence alignment to a chosen PDB polymer chain, and explicit
handling of alignment gaps and unresolved residues.
