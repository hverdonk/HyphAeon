# ChronAeon

**The Velocity of Time** — Ultra-fast molecular clock dating, phylodynamics, and genomic surveillance.

ChronAeon operates along the **TAXON / LINEAGE** axis of the Aeon foundation model, providing sub-second execution for tree-free continuous manifold dating, attention-derived covariance, and planetary-scale genomic screening.

## Target Audience

Public health agencies (CDC, WHO, UKHSA), outbreak epidemiologists, hospital infection control teams.

## Installation

```bash
pip install chronaeon
```

## CLI Subcommands

| Command | Aliases | Description |
| :--- | :--- | :--- |
| `chronaeon date` | `dating`, `clock`, `mrca` | Molecular clock calibration, tMRCA dating |
| `chronaeon autoclock` | `deconvolve`, `multiclock` | Hierarchical multi-clock deconvolution |
| `chronaeon triage` | `radar`, `sieve`, `qc` | Streaming genomic QC triage / outbreak radar |
| `chronaeon phylogeo` | `geo`, `spatial`, `dispersal` | Discrete phylogeography |
| `chronaeon dynamics` | `r0`, `rt`, `growth` | Phylodynamic R₀/Rₜ estimation |
| `chronaeon sketch` | `cluster`, `bin`, `centrifuge` | MinHash sketching & binning |
| `chronaeon align` | `thread`, `codon-align` | Reference-guided codon alignment |

## Key Capabilities

- **tMRCA Dating** — Continuous sequence manifold dating using attention-derived covariance (`compute_neural_covariance_kernel`)
- **AutoClock** — Hierarchical multi-clock deconvolution for complex evolutionary scenarios
- **Triage/Radar** — Stream 100k genomes in minutes, detect emerging clades, flag anomalous spillover branches
- **Phylogeography** — Continuous spatial dispersal reconstruction
- **Phylodynamics** — R₀/Rₜ growth rate estimation from heterochronous sequences

## The Radar & Microscope Flywheel

ChronAeon is the **Radar**: rapidly screens genomes, detects emerging clades, and infers origin dates. HyphAeon is the **Microscope**: dissects *why* flagged clades emerged — identifying positive selection bursts and epistatic rewiring. Together they form a collaborative flywheel for genomic surveillance and deep evolutionary analysis.

## Documentation

- [Dating Guide](DATING_GUIDE.md)
- [AutoClock Guide](AUTOCLOCK_GUIDE.md)
- [Triage Production Design](CHRONAEON_TRIAGE_PRODUCTION_DESIGN.md)
- [MRCA Dating Report](MRCA_DATING_REPORT.md)
