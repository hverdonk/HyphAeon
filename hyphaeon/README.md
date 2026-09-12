# HyphAeon

**The Depth of Evolution** — Neural selection inference using a phylogenetic foundation model.

HyphAeon operates along the **SITE** axis of the Aeon foundation model, providing episodic positive selection detection, 3D macromolecular epistasis, in silico DMS sweeps, and phenotype attribution.

## Target Audience

Comparative evolutionary biologists, structural virologists, protein biochemists.

## Installation

```bash
pip install hyphaeon
```

## CLI Subcommands

| Command | Aliases | Description |
| :--- | :--- | :--- |
| `hyphaeon meme` | `site-lrt` | Episodic positive selection (MEME/BUSTED) |
| `hyphaeon epistasis` | `cesi` | 3D macromolecular epistasis |
| `hyphaeon dms` | `essm`, `digital-dms` | In silico DMS sweeps |
| `hyphaeon busted` | `omnibus`, `gene-selection` | Gene-level selection |
| `hyphaeon disease` | `pathogenicity`, `variant`, `clinvar` | Disease pathogenicity prediction |
| `hyphaeon filter` | `mask`, `qc`, `clean` | Alignment QC / masking |
| `hyphaeon splits` | | Spectral bisection |
| `hyphaeon temporal` | | Temporal selection dynamics |
| `hyphaeon evaluate` | | HyPhy concordance evaluation |
| `hyphaeon list-models` | | List available model variants |

## Key Capabilities

- **Selection Inference** — MEME/BUSTED-style episodic positive selection via neural LRTs
- **Epistasis** — 3D macromolecular epistatic sector detection (CESI)
- **DMS** — In silico deep mutational scanning sweeps (ESSM)
- **Phenotype Attribution** — PhyloWAS/PARS phenotype-genotype associations
- **Spectral Splits** — Phylogenetic spectral bisection for clade analysis

## The Radar & Microscope Flywheel

HyphAeon is the **Microscope**: takes clades flagged by ChronAeon (the **Radar**) and dissects *why* they emerged — identifying positive selection bursts, functional epistatic rewiring, or host-jump phenotypic signatures. Together they form a collaborative flywheel for genomic surveillance and deep evolutionary analysis.

## Documentation

- [Architecture](ARCHITECTURE.md)
- [Migration Guide](MIGRATION_GUIDE.md)
- [Spectral Splits Benchmark](SPECTRAL_SPLITS_BENCHMARK.md)
- [Temporal Analysis Guide](TEMPORAL_ANALYSIS_GUIDE.md)
