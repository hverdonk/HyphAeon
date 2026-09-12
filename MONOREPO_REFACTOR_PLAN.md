# Monorepo Refactor Plan: Shared Core + Independent Packages

**Status:** Ready for implementation
**Base branch:** `feature/dating-mrca-module` on `veg/HyphAeon` (upstream remote)
**Target:** Refactor monolithic `hyphaeon` package into three independent PyPI packages within the same git repo: `aeon-core`, `hyphaeon`, and `chronaeon`.

---

## 1. Background

The `veg/HyphAeon` repository currently contains a single Python package (`hyphaeon/`) that serves two scientifically distinct purposes:

1. **Neural selection inference** (HyphAeon): Uses a phylogenetic transformer model (`PhyloAxialTransformer`) for detecting episodic positive selection, phenotype-genotype associations, and epistatic sectors. Targets deep-time multi-species alignments (200 Myr scale).

2. **Molecular clock & phylodynamics** (ChronAeon): Uses the same transformer model for attention-PGLS dating and latent manifold dating, plus classical-statistics methods for viral surveillance — heterochronous molecular clock dating, hierarchical multi-clock deconvolution, streaming genomic QC, phylogeography, R0/Rt estimation, MinHash sketching, and reference-guided codon alignment. Targets viral/pathogen surveillance (days to decades).

Both packages share the same neural model architecture, pretrained weights, and model utilities (model loading, alignment tensor preparation, device management, cross-taxa attention extraction). The model is a **hard dependency** for both — it is always used with pretrained weights.

The goal is to extract the shared model + infrastructure into `aeon-core`, and give each tool its own package, `pyproject.toml`, CLI entry point, version, and PyPI release cadence — all within the same git repo. Future `*aeon` tools can depend on `aeon-core` for the model and shared utilities.

### 1.1 Strategic Positioning: One Foundation Backbone, Two Orthogonal Axes

Both tools query different geometric projections of the same underlying 384-dimensional foundation model (`aeon-core`):

**HyphAeon — "The Depth of Evolution"**
- **Primary Axis:** Operates along the **SITE** axis.
- **Audience:** Comparative evolutionary biologists, structural virologists, protein biochemists.
- **Domain:** Episodic positive selection (MEME/BUSTED), 3D macromolecular epistasis (CESI), in silico DMS sweeps (ESSM), and phenotype attribution (PhyloWAS/PARS).
- **Persona:** Foundation model replacing multi-month Slurm cluster runs and numerical CTMC likelihood optimization.
- **Answers:** "Which codons are mutating, under selection, or structurally coupled?"

**ChronAeon — "The Velocity of Time"**
- **Primary Axis:** Operates along the **TAXON / LINEAGE** axis.
- **Audience:** Public health agencies (CDC, WHO, UKHSA), outbreak epidemiologists, hospital infection control teams.
- **Domain:** Real-time pathogen surveillance, tMRCA dating, phylodynamics (R₀/Rₜ), continuous spatial dispersal, and planetary-scale mega-screens (100k genomes in minutes).
- **Persona:** Sub-second execution, tree-free continuous manifold geometry, and zero-install client-side deployment.
- **Answers:** "When did this emerge, how fast is it spreading, and where is it going?"

**The Collaborative Flywheel ("Radar & Microscope")**

ChronAeon is the **Radar**: Rapidly screens 100k genomes in minutes, detects emerging clades, infers origin dates, and flags anomalous spillover branches. HyphAeon is the **Microscope**: Takes the specific clades flagged by ChronAeon and dissects *why* they emerged — identifying positive selection bursts, functional epistatic rewiring, or host-jump phenotypic signatures.

ChronAeon is not a classical distance/TempEst script; its core breakthrough is continuous sequence manifold dating and attention-derived covariance (`compute_neural_covariance_kernel`), which directly require the foundation model. This is why the model is a hard dependency for both packages.

---

## 2. Current State (on `feature/dating-mrca-module`)

### 2.1 File Inventory

**Top-level files:**
```
hyphaeon/                    # single package directory (18 .py files on main, 25 on dating branch)
tests/                       # 24 test files (dating branch)
model_eval/                  # concordance/calibration/invariance test suite for hyphaeon
scripts/                     # 12 scripts (some hyphaeon-specific, some chronaeon-specific)
examples/                    # shared example data
benchmarks/                  # hyphaeon benchmark data (e.g. benchmark_splits_vs_ml.csv)
docs/                        # GitHub Pages site
train.py                     # hyphaeon training script (top-level)
pyproject.toml               # single package config
model.safetensors            # pretrained weights
model_config.json            # model architecture config
```

### 2.2 Module Classification

**HyphAeon-specific modules** (stay in `hyphaeon` package):
| File | Purpose | Internal imports |
| :--- | :--- | :--- |
| `phenotype.py` | PhyloWAS phenotype-genotype association | `aeon_core.dataset`, `aeon_core.model`, `aeon_core.stats`, `aeon_core.inference`, `aeon_core.weights`, `.epistasis` |
| `epistasis.py` | Epistatic sector mining | `aeon_core.dataset`, `aeon_core.model`, `aeon_core.weights`, `aeon_core.stats`, `aeon_core.inference`, `aeon_core._progress` |
| `disease.py` | Disease pathogenicity prediction | `aeon_core.dataset`, `aeon_core.model`, `aeon_core.inference`, `aeon_core._progress` |
| `filter.py` | Alignment QC / masking | `aeon_core.dataset`, `aeon_core.model`, `aeon_core.weights`, `aeon_core.inference`, `aeon_core.stats`, `aeon_core.inference` (for `compute_adaptive_safe_batch_size`, was `.epistasis`) |
| `evaluation.py` | HyPhy MEME concordance evaluation | `aeon_core.io` |
| `attribution.py` | Branch attribution | `aeon_core.dataset`, `aeon_core._progress` |
| `training_data.py` | Training data loading | `aeon_core.dataset` (currently absolute `from hyphaeon.dataset import ...`) |
| `temporal.py` | Temporal selection dynamics (4 functions) | `aeon_core.inference`, `aeon_core.stats`, `aeon_core.io`, `aeon_core.temporal`, `.epistasis` |
| `inference.py` | Site-level LRT prediction (`predict_site_lrts` only) | `aeon_core.inference`, `aeon_core.model`, `aeon_core._progress` |
| `splits.py` | Spectral bisection CLI logic (`spectral_bisection`, `run_spectral_splits`, etc.) | `aeon_core.splits`, `aeon_core.inference`, `aeon_core.dataset` |

**ChronAeon modules** (move to `chronaeon` package):
| File | Lines | Purpose | Internal imports (after refactor) |
| :--- | ---: | :--- | :--- |
| `dating.py` | 3,194 | Molecular clock calibration, MRCA dating | `aeon_core.dataset`, `aeon_core.inference`, `aeon_core.splits`, `aeon_core.temporal`, `aeon_core.io` |
| `autoclock.py` | 1,878 | Hierarchical multi-clock deconvolution | `aeon_core.dataset`, `.dating`; conditional: `aeon_core.inference`, `aeon_core.splits`, `.dating` (for `compute_neural_covariance_kernel`) |
| `triage.py` (was `sieve.py`) | 803 | Streaming genomic QC triage / outbreak radar | `aeon_core.dataset`, `aeon_core.temporal`, `.dating` |
| `geo.py` | 921 | Discrete phylogeography | `aeon_core.dataset`; conditional: `aeon_core.inference`, `aeon_core.splits` |
| `r0.py` | 935 | Phylodynamic R0/Rt estimation | none |
| `sketch.py` | 302 | MinHash sketching & binning | none |
| `alignment.py` | 200 | Reference-guided codon threading | none |

**Shared infrastructure** (extract to `aeon-core` package):
| File | Functions to extract | Used by |
| :--- | :--- | :--- |
| `dataset.py` | **All functions** — `parse_alignment_sequences`, `compute_tn93_distance_matrix`, `compute_tn93_cross_distance_matrix`, `parse_beast_xml`, `load_alignment_and_tree`, `GENETIC_CODE`, `CODON_TO_AA`, `AA_MAP`, `get_codon_token`, `get_aa_token`, `extract_tree_from_string_or_file`, `has_nonzero_branch_lengths`, `estimate_tree_branch_lengths_hyphy`, `enforce_nonzero_branch_lengths`, `compute_fast_dist_matrix`, `compute_mds_coordinates`, `downsample_taxa_faith_pd`, `prune_identical_sequences` | hyphaeon (phenotype, epistasis, filter, training_data), chronaeon (dating, autoclock, triage, geo) |
| `model.py` | **All** — `PhyloAxialTransformer`, `BustedMultiTaskHead`, `decode_soft_ordinal_lrt` | hyphaeon (inference, phenotype, epistasis, disease, filter, attribution), chronaeon (dating, autoclock, geo) |
| `weights.py` | **All** — `resolve_weights_path`, `load_arch_config`, `load_weights`, `load_model_config`, `list_available_variants`, `get_variant_filename`, `HF_REPO_ID`, `DEFAULT_VARIANT`, `CACHE_DIR` | hyphaeon (inference, cli), chronaeon (dating, autoclock, geo) |
| `inference.py` | `get_device`, `get_device_memory_budget`, `compute_adaptive_safe_batch_size`, `load_model`, `prepare_alignment` | hyphaeon (cli, epistasis, phenotype, disease, filter), chronaeon (dating, autoclock, geo) |
| `splits.py` | `extract_cross_taxa_attentions_and_embeddings`, `compute_fused_affinity_matrix` | hyphaeon (cli for `splits` subcommand), chronaeon (dating) |
| `temporal.py` | `parse_date_to_decimal`, `parse_dates_from_auspice_json`, `extract_date_from_string` | chronaeon (dating, triage), hyphaeon (temporal dynamics) |
| `io.py` | `ensure_parent_directory`, `write_json`, `write_csv`, `format_pq` | hyphaeon (cli), chronaeon (dating) |
| `stats.py` | `pvals_from_lrt_meme`, `pvals_from_lrt_self_liang`, `benjamini_hochberg`, `cauchy_combination_p` | hyphaeon (cli, epistasis, phenotype) |
| `_progress.py` | `ChunkProgress` | hyphaeon (cli), aeon-core (inference) |

**Note on `inference.py` split:** `predict_site_lrts` is the only function that stays in `hyphaeon/inference.py` — it runs the site-level LRT prediction loop specific to selection inference. All other functions (`get_device`, `load_model`, `prepare_alignment`, `compute_adaptive_safe_batch_size`, `get_device_memory_budget`) are generic model utilities that move to `aeon_core/inference.py`.

**Note on `splits.py` split:** `extract_cross_taxa_attentions_and_embeddings` and `compute_fused_affinity_matrix` are model-attention utilities used by both packages — they move to `aeon_core/splits.py`. The spectral bisection functions (`spectral_bisection`, `tree_dict_to_newick`, `get_all_clade_taxa`, `run_spectral_splits`) are HyphAeon-specific CLI logic — they stay in `hyphaeon/splits.py`.

**Note on `temporal.py` split:** The 3 date-parsing functions move to `aeon_core/temporal.py`. The 4 selection-dynamics functions (`parse_temporal_metadata`, `infer_root_sequence`, `run_temporal_surveillance`, `render_temporal_4panel_figure`) stay in `hyphaeon/temporal.py`.

**Note on `dataset.py`:** Moves entirely to `aeon_core/dataset.py`. The `load_alignment_and_tree` function prepares PyTorch tensors for the transformer model, but since the model itself is in `aeon-core`, this tensor preparation belongs there too.

### 2.3 CLI Subcommands (current `hyphaeon/cli.py`)

**HyphAeon subcommands** (stay):
- `meme` (aliases: `predict`, `site-selection`)
- `phenotype` (aliases: `phylowas`, `trait`)
- `epistasis` (aliases: `coselection`, `sector`, `network`)
- `dms` (aliases: `essm`, `digital-dms`)
- `busted` (aliases: `omnibus`, `gene-selection`)
- `disease` (aliases: `pathogenicity`, `variant`, `clinvar`)
- `filter` (aliases: `mask`, `qc`, `clean`)
- `splits`
- `temporal`
- `evaluate`
- `list-models`

**ChronAeon subcommands** (move to `chronaeon/cli.py`):
- `date` (aliases: `dating`, `clock`, `mrca`)
- `autoclock` (aliases: `deconvolve`, `multiclock`)
- `triage` (aliases: `radar`, `sieve`, `qc`)
- `phylogeo` (aliases: `geo`, `spatial`, `dispersal`)
- `dynamics` (aliases: `r0`, `rt`, `growth`)
- `sketch` (aliases: `cluster`, `bin`, `centrifuge`)
- `align` (aliases: `thread`, `codon-align`)

### 2.4 Test Classification

**HyphAeon tests** (stay in `hyphaeon/tests/`):
- `test_busted.py`, `test_disease_filter_cli.py` (note: currently empty), `test_epistasis.py`, `test_evaluation.py`, `test_integration.py`, `test_phenotype.py`, `test_splits.py` (spectral bisection tests only), `test_temporal.py` (selection dynamics tests only), `test_temporal_experimental.py`, `test_training_data.py`, `test_training.py`
- `conftest.py` (imports from `aeon_core.model`, creates `dummy_weights` fixture)
- `test_imports.py` (new — pre-migration import smoke test)

**ChronAeon tests** (move to `chronaeon/tests/`):
- `test_dating.py`, `test_geo.py`, `test_r0.py`, `test_triage.py` (was `test_sieve.py`), `test_sketch_and_alignment.py`

**Shared tests** (move to `aeon-core/tests/`):
- `test_weights.py` → `aeon-core/tests/test_weights.py`
- `test_batch_size.py` → `aeon-core/tests/test_inference.py`
- `test_gpu_mem_guard.py` → `aeon-core/tests/test_inference.py`
- `test_tokenization.py` → `aeon-core/tests/test_model.py`
- `test_distance_mds.py` → `aeon-core/tests/test_dataset.py`
- `test_stats.py` → `aeon-core/tests/test_stats.py`
- `test_load_alignment.py` → `aeon-core/tests/test_dataset.py` (moves entirely — `dataset.py` including `load_alignment_and_tree` is now in `aeon-core`)
- `test_parsing.py` → `aeon-core/tests/test_dataset.py` (moves entirely — all parsing functions are in `aeon-core`)
- `test_splits.py` — split: tests for `extract_cross_taxa_attentions_and_embeddings`, `compute_fused_affinity_matrix` → `aeon-core/tests/test_splits.py`; tests for `spectral_bisection`, `run_spectral_splits` → stay in `hyphaeon/tests/test_splits.py`
- `test_temporal.py` — split: tests for the 3 date-parsing functions → `aeon-core/tests/test_temporal.py`; tests for selection dynamics → stay in `hyphaeon/tests/test_temporal.py`
- Create `aeon-core/tests/test_io.py` with tests for I/O helpers
- Create `aeon-core/tests/test_imports.py` (new — pre-migration import smoke test)
- Create `aeon-core/tests/conftest.py` (imports from `aeon_core.model` — needs torch)

### 2.5 Script Classification

**HyphAeon scripts** (stay):
- `build_training_npz.py`, `calc_masked_fractions.py`, `monitor_and_review_prs.py`, `plot_avian_artifact_context_multipanel.py`, `plot_yokoyama_rhodopsin_benchmark.py`, `recreate_masked_alignments.py`, `run_avian_genome_filter.py`, `run_orthomam_masked_selection.py`

**ChronAeon scripts** (move to `chronaeon/scripts/`):
- `audit_empirical_sus_genomes.py`, `harvest_sars2_bvbrc.py`, `run_sars2_streaming_triage.py` (was `run_sars2_streaming_sieve.py`), `stream_10k_sars2_bvbrc.py`

### 2.6 Documentation Classification

**HyphAeon docs** (stay at repo root or move to `hyphaeon/`):
- `ARCHITECTURE.md`, `MIGRATION_GUIDE.md`, `SPECTRAL_SPLITS_BENCHMARK.md`, `TEMPORAL_ANALYSIS_GUIDE.md`, `PREPRINT_REVIEW.md`

**ChronAeon docs** (move to `chronaeon/` or keep at root with `chronaeon-` prefix):
- `DATING_GUIDE.md`, `AUTOCLOCK_GUIDE.md`, `CHRONAEON_TRIAGE_PRODUCTION_DESIGN.md` (was `CHRONAEON_SIEVE_PRODUCTION_DESIGN.md`), `MRCA_DATING_REPORT.md`

**Shared** (stay at root):
- `README.md` (rewrite as monorepo overview), `LICENSE`, `CITATION.cff` (split into independent citation entries per package)

**Independent package READMEs** (new, per package):
- `aeon-core/README.md` — foundation model description, shared utilities
- `hyphaeon/README.md` — selection inference focus, MEME/BUSTED/epistasis tutorials
- `chronaeon/README.md` — phylodynamics/outbreak surveillance focus, dating/R₀/triage tutorials

---

## 3. Target Directory Layout

```
HyphAeon/                              # monorepo root (git repo stays here)
├── aeon-core/                         # shared infrastructure + model package
│   ├── pyproject.toml                 # name = "aeon-core", deps include torch, huggingface-hub, safetensors
│   ├── README.md
│   ├── src/
│   │   └── aeon_core/
│   │       ├── __init__.py
│   │       ├── model.py               # PhyloAxialTransformer, BustedMultiTaskHead
│   │       ├── weights.py             # HuggingFace weight resolution, resolve_weights_path, load_model
│   │       ├── inference.py           # get_device, load_model, prepare_alignment, batch sizing
│   │       ├── splits.py              # extract_cross_taxa_attentions_and_embeddings, compute_fused_affinity_matrix
│   │       ├── dataset.py             # alignment parsing, TN93, genetic code, BEAST XML, load_alignment_and_tree
│   │       ├── temporal.py            # date parsing (3 functions)
│   │       ├── io.py                  # file I/O helpers
│   │       ├── stats.py               # p-value utilities
│   │       └── _progress.py           # ChunkProgress
│   └── tests/
│       ├── conftest.py
│       ├── test_dataset.py
│       ├── test_imports.py
│       ├── test_inference.py
│       ├── test_io.py
│       ├── test_model.py
│       ├── test_splits.py
│       ├── test_stats.py
│       ├── test_temporal.py
│       └── test_weights.py
│
├── hyphaeon/                           # neural selection inference package
│   ├── pyproject.toml                  # name = "hyphaeon", deps = ["aeon-core>=0.1"]
│   ├── ARCHITECTURE.md
│   ├── MIGRATION_GUIDE.md
│   ├── SPECTRAL_SPLITS_BENCHMARK.md
│   ├── TEMPORAL_ANALYSIS_GUIDE.md
│   ├── PREPRINT_REVIEW.md
│   ├── src/
│   │   └── hyphaeon/
│   │       ├── __init__.py
│   │       ├── cli.py                  # hyphaeon CLI (meme, phenotype, epistasis, etc.)
│   │       ├── inference.py            # predict_site_lrts only (site-level LRT loop)
│   │       ├── splits.py               # spectral_bisection, run_spectral_splits (CLI logic)
│   │       ├── temporal.py             # temporal selection dynamics (4 functions)
│   │       ├── phenotype.py
│   │       ├── epistasis.py
│   │       ├── disease.py
│   │       ├── filter.py
│   │       ├── evaluation.py
│   │       ├── attribution.py
│   │       ├── training_data.py
│   │       └── model_eval/
│   │           ├── __init__.py
│   │           ├── _harness.py
│   │           ├── _sim.py
│   │           ├── conftest.py
│   │           ├── calibration/
│   │           ├── concordance/
│   │           ├── invariance/
│   │           ├── reports/
│   │           └── stability/
│   └── tests/
│       ├── conftest.py
│       ├── test_busted.py
│       ├── test_disease_filter_cli.py
│       ├── test_epistasis.py
│       ├── test_evaluation.py
│       ├── test_imports.py
│       ├── test_integration.py
│       ├── test_phenotype.py
│       ├── test_splits.py             # spectral bisection tests only
│       ├── test_temporal.py            # selection dynamics tests only
│       ├── test_temporal_experimental.py
│       ├── test_training_data.py
│       └── test_training.py
│   └── benchmarks/                     # hyphaeon benchmark data
│       └── benchmark_splits_vs_ml.csv
│
├── chronaeon/                          # molecular clock & phylodynamics package
│   ├── pyproject.toml                  # name = "chronaeon", deps = ["aeon-core>=0.1"] (no optional extras needed)
│   ├── README.md
│   ├── src/
│   │   └── chronaeon/
│   │       ├── __init__.py
│   │       ├── cli.py                  # chronaeon CLI (date, autoclock, triage, phylogeo, etc.)
│   │       ├── dating.py
│   │       ├── autoclock.py
│   │       ├── triage.py               # was sieve.py — streaming QC triage / outbreak radar
│   │       ├── geo.py
│   │       ├── r0.py
│   │       ├── sketch.py
│   │       └── alignment.py
│   ├── examples/                       # chronaeon-specific example data
│   │   ├── H1N1_2009_pandemic.fasta
│   │   ├── H1N1_2009_pandemic.nwk
│   │   ├── H5N1_HA_geo.fasta
│   │   ├── H5N1_HA.nwk
│   │   ├── H5N1_HA_metadata.csv
│   │   └── korber_env_gp160.fasta
│   ├── scripts/
│   │   ├── audit_empirical_sus_genomes.py
│   │   ├── harvest_sars2_bvbrc.py
│   │   ├── run_sars2_streaming_triage.py
│   │   └── stream_10k_sars2_bvbrc.py
│   ├── DATING_GUIDE.md
│   ├── AUTOCLOCK_GUIDE.md
│   ├── CHRONAEON_TRIAGE_PRODUCTION_DESIGN.md
│   ├── MRCA_DATING_REPORT.md
│   ├── mrca_dating_benchmark.png
│   ├── mrca_dating_diagnostic.pdf
│   ├── mrca_dating_diagnostic.png
│   ├── ebola_r0_example.png
│   ├── h1n1_r0_example.png
│   ├── h5n1_geo_example.png
│   └── tests/
│       ├── conftest.py
│       ├── test_autoclock_imports.py
│       ├── test_dating.py
│       ├── test_geo.py
│       ├── test_imports.py
│       ├── test_r0.py
│       ├── test_triage.py
│       └── test_sketch_and_alignment.py
│
├── examples/                           # shared example data (stays at root)
├── scripts/                            # hyphaeon-specific scripts stay here
├── docs/                               # shared documentation site
├── .github/workflows/                  # CI
├── train.py                            # training script (stays at root, imports from aeon_core)
├── model.safetensors                   # pretrained weights (stays at root, used by aeon_core.weights)
├── model_config.json                   # model config (stays at root, used by aeon_core.weights)
├── README.md                           # monorepo overview
├── LICENSE
├── CITATION.cff
└── .gitignore
```

---

## 4. Detailed Refactor Steps

### Step 1: Create `aeon-core/` package

**1.1. Create directory structure:**
```
aeon-core/
├── pyproject.toml
├── src/aeon_core/
│   ├── __init__.py
│   ├── model.py
│   ├── weights.py
│   ├── inference.py
│   ├── splits.py
│   ├── dataset.py
│   ├── temporal.py
│   ├── io.py
│   ├── stats.py
│   └── _progress.py
└── tests/
    ├── conftest.py
    ├── test_dataset.py
    ├── test_imports.py
    ├── test_inference.py
    ├── test_io.py
    ├── test_model.py
    ├── test_progress.py
    ├── test_splits.py
    ├── test_stats.py
    ├── test_temporal.py
    └── test_weights.py
```

**1.2. Create `aeon-core/pyproject.toml`:**
```toml
[build-system]
requires = ["setuptools>=61.0"]
build-backend = "setuptools.build_meta"

[project]
name = "aeon-core"
version = "0.1.0"
authors = [
  { name="Sergei L. Kosakovsky Pond", email="spond@temple.edu" },
]
description = "Aeon-Core: Shared bioinformatics infrastructure for Aeon-family packages"
readme = "README.md"
requires-python = ">=3.8"
license = { text = "MIT" }
classifiers = [
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]
dependencies = [
    "torch>=2.0.0",
    "biopython>=1.80",
    "numpy>=1.22.0",
    "scipy>=1.9.0",
    "pandas>=1.5.0",
    "huggingface-hub>=0.20.0",
    "safetensors>=0.4.0",
]

[project.optional-dependencies]
tn93 = ["tn93>=1.2.0"]
dev = ["pytest>=7.0"]

[tool.setuptools]
packages = ["aeon_core"]

[tool.setuptools.package-dir]
aeon_core = "src/aeon_core"
```

**1.3. Move `model.py` to `aeon_core/model.py`:**

Move the entire file. It contains `PhyloAxialTransformer`, `BustedMultiTaskHead`, `decode_soft_ordinal_lrt`, and related classes. No internal imports — only `torch`, `numpy`, `math`.

**1.4. Move `weights.py` to `aeon_core/weights.py`:**

Move the entire file. It contains `resolve_weights_path`, `load_arch_config`, `load_weights`, `load_model_config`, `list_available_variants`, `get_variant_filename`, `HF_REPO_ID`, `DEFAULT_VARIANT`, `DEFAULT_CONFIG_FILENAME`, `CACHE_DIR`, `_DEFAULT_ARCH`. Imports: `os`, `json`, `pathlib`, `huggingface_hub`, `safetensors`.

**1.5. Split `inference.py`:**

Move these functions to `aeon_core/inference.py`:
- `get_device(cpu: bool = False) -> torch.device`
- `get_device_memory_budget(device: torch.device) -> float`
- `compute_adaptive_safe_batch_size(...)`
- `load_model(weights=None, variant=None, device=None, strict=False)`
- `prepare_alignment(alignment_path, tree_path=None, model=None, device=None, ...)`

Keep in `hyphaeon/inference.py`:
- `predict_site_lrts(model, c, a, d, z, inv, tree_cache=None, ...)` — this is the site-level LRT prediction loop specific to selection inference.

`aeon_core/inference.py` imports: `aeon_core.model`, `aeon_core.weights`, `aeon_core.dataset`, `aeon_core._progress`.
`hyphaeon/inference.py` imports: `aeon_core.inference` (for `get_device`, etc.), `aeon_core.model`.

**1.6. Split `splits.py`:**

Move these functions to `aeon_core/splits.py`:
- `extract_cross_taxa_attentions_and_embeddings(model, msa_codons, msa_aas, tree_cache, device=None, ...)`
- `compute_fused_affinity_matrix(cross_attn, taxa_repr, ...)`

Keep in `hyphaeon/splits.py`:
- `spectral_bisection(...)`
- `tree_dict_to_newick(tree_dict)`
- `get_all_clade_taxa(node)`
- `run_spectral_splits(...)`

`aeon_core/splits.py` imports: `aeon_core.inference` (for `load_model`, `get_device`), `aeon_core.dataset` (for `load_alignment_and_tree`).
`hyphaeon/splits.py` imports: `aeon_core.splits` (for `extract_cross_taxa_attentions_and_embeddings`, `compute_fused_affinity_matrix`), `aeon_core.inference`.

**1.7. Move `dataset.py` to `aeon_core/dataset.py`:**

Move the entire file including `load_alignment_and_tree`. It contains:
- `GENETIC_CODE`, `CODON_TO_AA`, `AA_MAP` (constant dicts)
- `get_codon_token(codon: str) -> int`
- `get_aa_token(codon: str) -> int`
- `parse_alignment_sequences(filepath: str) -> Dict[str, str]`
- `compute_tn93_distance_matrix(...)`
- `compute_tn93_cross_distance_matrix(...)` (dating-branch addition)
- `parse_beast_xml(filepath: Union[str, Path]) -> Dict[str, Any]` (dating-branch addition)
- `extract_tree_from_string_or_file(source: str) -> Optional[Phylo.BaseTree.Tree]`
- `has_nonzero_branch_lengths(tree, ...) -> bool`
- `estimate_tree_branch_lengths_hyphy(seq_dict, tree_obj) -> Optional[Phylo.BaseTree.Tree]`
- `enforce_nonzero_branch_lengths(tree_obj, min_len, default_missing) -> Phylo.BaseTree.Tree`
- `compute_fast_dist_matrix(tree, taxa) -> np.ndarray`
- `compute_mds_coordinates(dist_matrix, n_components) -> np.ndarray`
- `downsample_taxa_faith_pd(dist_mat, taxa, max_species) -> Tuple[np.ndarray, List[str]]`
- `prune_identical_sequences(seq_dict, taxa) -> Tuple[List[str], Dict[str, List[str]], int]`
- `load_alignment_and_tree(...)` — prepares PyTorch tensors for the transformer model
- `_warn_internal_stops(seq_dict, L, taxa) -> None`
- `_parse_numeric_or_calendar_date(v_str) -> Optional[float]`

**1.8. Split `temporal.py`:**

Move these 3 functions to `aeon_core/temporal.py`:
- `parse_date_to_decimal(val: Any, time_units: str = "years") -> float`
- `parse_dates_from_auspice_json(json_path: Union[str, Path]) -> Dict[str, float]`
- `extract_date_from_string(name: str, time_units: str = "years") -> float`

Keep in `hyphaeon/temporal.py`:
- `parse_temporal_metadata(...)`
- `infer_root_sequence(...)`
- `run_temporal_surveillance(...)`
- `render_temporal_4panel_figure(...)`

The 3 date-parsing functions have no internal imports (only stdlib + numpy). The 4 selection-dynamics functions import from `.epistasis` (`REV_AA_MAP`, `AA_MAP`), `.inference` (`load_model`, `get_device`, `prepare_alignment`, `compute_adaptive_safe_batch_size`), `.stats` (`pvals_from_lrt_self_liang`, `benjamini_hochberg`), `.io` (`ensure_parent_directory`, `write_json`, `write_csv`), and optionally `matplotlib`. After the split, `hyphaeon/temporal.py` should change these to:
```python
from aeon_core.inference import load_model, get_device, prepare_alignment, compute_adaptive_safe_batch_size
from aeon_core.stats import pvals_from_lrt_self_liang, benjamini_hochberg
from aeon_core.io import ensure_parent_directory, write_json, write_csv
from .epistasis import REV_AA_MAP, AA_MAP  # stays in hyphaeon
```

**1.9. Move `io.py` to `aeon_core/io.py`:**

Move all 4 functions:
- `ensure_parent_directory(path)`
- `write_json(path, data, label="JSON results")`
- `write_csv(path, df_or_records, label="CSV results")`
- `format_pq(value)`

**1.10. Move `stats.py` to `aeon_core/stats.py`:**

Move all 4 functions:
- `pvals_from_lrt_meme(lrts: np.ndarray) -> np.ndarray`
- `pvals_from_lrt_self_liang(lrts: np.ndarray) -> np.ndarray`
- `benjamini_hochberg(pvals: np.ndarray) -> np.ndarray`
- `cauchy_combination_p(pvals: np.ndarray) -> float`

**1.11. Move `_progress.py` to `aeon_core/_progress.py`:**

Move the entire file. It contains `ChunkProgress` class. No internal imports — only `sys` and `time`.

**1.12. Create `aeon_core/__init__.py`:**

```python
"""
Aeon-Core: Shared model and bioinformatics infrastructure for Aeon-family packages.
"""

from .model import PhyloAxialTransformer, BustedMultiTaskHead, decode_soft_ordinal_lrt
from .weights import (
    resolve_weights_path, load_arch_config, load_weights,
    load_model_config, list_available_variants, get_variant_filename,
    HF_REPO_ID, DEFAULT_VARIANT, CACHE_DIR,
)
from .inference import (
    get_device, get_device_memory_budget,
    compute_adaptive_safe_batch_size,
    load_model, prepare_alignment,
)
from .splits import (
    extract_cross_taxa_attentions_and_embeddings,
    compute_fused_affinity_matrix,
)
from .dataset import (
    GENETIC_CODE, CODON_TO_AA, AA_MAP,
    get_codon_token, get_aa_token,
    parse_alignment_sequences,
    compute_tn93_distance_matrix,
    compute_tn93_cross_distance_matrix,
    parse_beast_xml,
    extract_tree_from_string_or_file,
    has_nonzero_branch_lengths,
    estimate_tree_branch_lengths_hyphy,
    enforce_nonzero_branch_lengths,
    compute_fast_dist_matrix,
    compute_mds_coordinates,
    downsample_taxa_faith_pd,
    prune_identical_sequences,
    load_alignment_and_tree,
)
from .temporal import (
    parse_date_to_decimal,
    parse_dates_from_auspice_json,
    extract_date_from_string,
)
from .io import ensure_parent_directory, write_json, write_csv, format_pq
from .stats import (
    pvals_from_lrt_meme,
    pvals_from_lrt_self_liang,
    benjamini_hochberg,
    cauchy_combination_p,
)

__version__ = "0.1.0"
```

**1.13. Move shared tests to `aeon-core/tests/`:**

- Move `test_stats.py` entirely to `aeon-core/tests/test_stats.py`
- Move `test_weights.py` to `aeon-core/tests/test_weights.py`
- Move `test_batch_size.py` to `aeon-core/tests/test_inference.py` (tests `compute_adaptive_safe_batch_size`)
- Move `test_gpu_mem_guard.py` to `aeon-core/tests/test_inference.py` (tests `get_device_memory_budget`)
- Move `test_tokenization.py` to `aeon-core/tests/test_model.py` (tests `PhyloAxialTransformer` tokenization)
- Move `test_distance_mds.py` to `aeon-core/tests/test_dataset.py` (tests `compute_mds_coordinates`, `compute_fast_dist_matrix`)
- Split `test_load_alignment.py`: tests for `parse_alignment_sequences`, `compute_tn93_distance_matrix`, `load_alignment_and_tree` → `aeon-core/tests/test_dataset.py`
- Split `test_parsing.py`: tests for generic parsing → `aeon-core/tests/test_dataset.py`
- Split `test_splits.py`: tests for `extract_cross_taxa_attentions_and_embeddings`, `compute_fused_affinity_matrix` → `aeon-core/tests/test_splits.py`; tests for `spectral_bisection`, `run_spectral_splits` → stay in `hyphaeon/tests/test_splits.py`
- Create `aeon-core/tests/test_temporal.py` with tests for the 3 date-parsing functions
- Create `aeon-core/tests/test_io.py` with tests for I/O helpers
- Create `aeon-core/tests/conftest.py` — imports from `aeon_core.model` (needs torch for model instantiation tests)

---

### Step 1b: Pre-Migration Test Gap Analysis

Before starting the refactor, assess test coverage to ensure regressions will be caught. The migration is primarily an import-rewrite exercise, so **import breakage** is the primary risk — but some modules have thin functional coverage that could hide subtle regressions.

**Coverage audit (on `feature/dating-mrca-module`):**

| Module | Test file(s) | # tests | Coverage assessment |
| :--- | :--- | :--- | :--- |
| `dataset.py` | `test_load_alignment.py`, `test_parsing.py`, `test_distance_mds.py` | 25+ | **Good** — parsing, TN93, MDS, tensor prep |
| `model.py` | `test_busted.py`, `test_tokenization.py` | 15+ | **Good** — BustedMultiTaskHead, tokenization |
| `weights.py` | `test_weights.py` | 10 | **Good** — variant resolution, caching, safetensors |
| `inference.py` | `test_batch_size.py`, `test_gpu_mem_guard.py` | 15+ | **Good** — adaptive batch, memory budget |
| `splits.py` | `test_splits.py` | 10 | **Good** — symmetry, clustering, bisection |
| `stats.py` | `test_stats.py` | 10 | **Good** — LRT p-values, BH, Cauchy |
| `temporal.py` | `test_temporal.py`, `test_temporal_experimental.py` | 10+ | **Good** — date parsing + dynamics |
| `io.py` | *(none)* | 0 | **Gap** — no dedicated tests |
| `_progress.py` | *(none)* | 0 | **Gap** — no dedicated tests (low risk, simple) |
| `dating.py` | `test_dating.py` | 12 | **Good** — OLS, spline, MRCA, BEAST XML |
| `autoclock.py` | `test_dating.py` (2 tests) | 2 | **Thin** — 1,878 lines, 2 tests |
| `triage.py` (was `sieve.py`) | `test_triage.py` (was `test_sieve.py`) | 1 | **Thin** — 803 lines, 1 test |
| `geo.py` | `test_geo.py` | 6 | **Adequate** |
| `r0.py` | `test_r0.py` | 7 | **Good** |
| `sketch.py` | `test_sketch_and_alignment.py` | 3 | **Adequate** |
| `alignment.py` | `test_sketch_and_alignment.py` | 1 | **Thin** — 200 lines, 1 test |
| `phenotype.py` | `test_phenotype.py` | 10+ | **Good** |
| `epistasis.py` | `test_epistasis.py` | 10+ | **Good** |
| `disease.py` | `test_disease_filter_cli.py` | 0 | **Gap** — test file exists but empty |
| `filter.py` | *(none)* | 0 | **Gap** — covered indirectly by integration tests |
| `evaluation.py` | `test_evaluation.py` | 5 | **Adequate** |
| `attribution.py` | *(none)* | 0 | **Gap** — covered indirectly by epistasis tests |
| `training_data.py` | `test_training_data.py` | 10+ | **Good** |
| `training.py` | `test_training.py` | 7 | **Good** |
| Integration | `test_integration.py` | 6 | **Good** — full CLI pipelines with dummy weights |
| model_eval | 10 files across calibration/concordance/invariance/stability | 20+ | **Good** — real-weight evaluation suite |

**Pre-migration test additions (required):**

1. **`aeon-core/tests/test_io.py`** — Add 4 smoke tests for `write_json`, `write_csv`, `ensure_parent_directory`, `format_pq`. Currently zero coverage; if import paths break, only CLI integration tests catch it indirectly.

2. **Import smoke test for each package** — Add a `test_imports.py` to each package's test suite that imports every module from its new location:
   ```python
   # aeon-core/tests/test_imports.py
   def test_aeon_core_imports():
       from aeon_core.model import PhyloAxialTransformer
       from aeon_core.weights import resolve_weights_path
       from aeon_core.inference import load_model, get_device
       from aeon_core.splits import extract_cross_taxa_attentions_and_embeddings
       from aeon_core.dataset import load_alignment_and_tree
       from aeon_core.temporal import parse_date_to_decimal
       from aeon_core.io import write_json, write_csv
       from aeon_core.stats import pvals_from_lrt_meme
       from aeon_core._progress import ChunkProgress

   # chronaeon/tests/test_imports.py
   def test_chronaeon_imports():
       from chronaeon.dating import run_mrca_dating, run_pgls_dating
       from chronaeon.autoclock import run_autoclock_deconvolution
       from chronaeon.triage import ChronAeonSieve
       from chronaeon.geo import run_phylogeography_analysis
       from chronaeon.r0 import run_r0_analysis
       from chronaeon.sketch import CanonicalMinHashSketcher
       from chronaeon.alignment import ReferenceCodonAligner

   # hyphaeon/tests/test_imports.py
   def test_hyphaeon_imports():
       from hyphaeon.inference import predict_site_lrts
       from hyphaeon.splits import run_spectral_splits
       from hyphaeon.temporal import run_temporal_surveillance
       from hyphaeon.phenotype import run_phenotype_association
       from hyphaeon.epistasis import run_epistasis_analysis
       from hyphaeon.disease import predict_disease_pathogenicity
       from hyphaeon.filter import run_alignment_filter
       from hyphaeon.evaluation import command as evaluate_command
       from hyphaeon.attribution import attribute_selection
       from hyphaeon.training_data import build_training_directory
   ```

**Pre-migration test additions (recommended but not blocking):**

3. **`aeon-core/tests/test_progress.py`** — 2 smoke tests for `ChunkProgress` (instantiation + update/finish cycle). Low risk but trivial to add.

4. **`chronaeon/tests/test_autoclock_imports.py`** — 1 test that explicitly exercises the conditional import block in `autoclock.py` (lines 376-378) where `from .inference` / `from .splits` are inside a `try` block. The existing 2 tests do trigger this path, but a dedicated import-only test makes the migration boundary explicit.

**Verdict:** The existing tests are **sufficient for catching import-path regressions** — every module has at least one test that imports it, and `test_integration.py` runs full CLI pipelines. The thin functional coverage on `autoclock.py`, `triage.py`, and `alignment.py` is a pre-existing issue, not a migration-specific risk. The `io.py` gap is the only real blind spot.

---

### Step 2: Create `chronaeon/` package

**2.1. Create directory structure:**
```
chronaeon/
├── pyproject.toml
├── src/chronaeon/
│   ├── __init__.py
│   ├── cli.py
│   ├── dating.py
│   ├── autoclock.py
│   ├── triage.py
│   ├── geo.py
│   ├── r0.py
│   ├── sketch.py
│   └── alignment.py
└── tests/
    ├── conftest.py
    ├── test_autoclock_imports.py
    ├── test_dating.py
    ├── test_geo.py
    ├── test_imports.py
    ├── test_r0.py
    ├── test_triage.py
    └── test_sketch_and_alignment.py
```

**2.2. Create `chronaeon/pyproject.toml`:**
```toml
[build-system]
requires = ["setuptools>=61.0"]
build-backend = "setuptools.build_meta"

[project]
name = "chronaeon"
version = "0.1.0"
authors = [
  { name="Sergei L. Kosakovsky Pond", email="spond@temple.edu" },
]
description = "ChronAeon: Ultra-Fast Molecular Clock Dating, Phylodynamics, and Genomic Surveillance"
readme = "README.md"
requires-python = ">=3.8"
license = { text = "MIT" }
classifiers = [
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]
dependencies = [
    "aeon-core>=0.1.0",
    "torch>=2.0.0",
    "biopython>=1.80",
    "numpy>=1.22.0",
    "scipy>=1.9.0",
    "pandas>=1.5.0",
    "scikit-learn>=1.0",
    "matplotlib>=3.5",
]

[project.optional-dependencies]
tn93 = ["tn93>=1.2.0"]
dev = ["pytest>=7.0"]

[project.scripts]
chronaeon = "chronaeon.cli:main"

[tool.setuptools]
packages = ["chronaeon"]

[tool.setuptools.package-dir]
chronaeon = "src/chronaeon"
```

**2.3. Move ChronAeon modules:**

Move these files from `hyphaeon/` to `chronaeon/src/chronaeon/`:
- `dating.py`
- `autoclock.py`
- `triage.py` (was `sieve.py`)
- `geo.py`
- `r0.py`
- `sketch.py`
- `alignment.py`

**2.4. Update imports in moved files:**

In `dating.py`, change:
```python
# OLD
from .dataset import (
    parse_alignment_sequences,
    compute_tn93_distance_matrix,
    compute_tn93_cross_distance_matrix,
    load_alignment_and_tree,
    parse_beast_xml,
    GENETIC_CODE,
    CODON_TO_AA,
)
from .inference import (
    load_model,
    get_device,
    prepare_alignment,
)
from .splits import (
    extract_cross_taxa_attentions_and_embeddings,
    compute_fused_affinity_matrix,
)
from .temporal import parse_date_to_decimal, parse_dates_from_auspice_json, extract_date_from_string
from .io import ensure_parent_directory, write_json, write_csv

# NEW
from aeon_core.dataset import (
    parse_alignment_sequences,
    compute_tn93_distance_matrix,
    compute_tn93_cross_distance_matrix,
    load_alignment_and_tree,
    parse_beast_xml,
    GENETIC_CODE,
    CODON_TO_AA,
)
from aeon_core.inference import (
    load_model,
    get_device,
    prepare_alignment,
)
from aeon_core.splits import (
    extract_cross_taxa_attentions_and_embeddings,
    compute_fused_affinity_matrix,
)
from aeon_core.temporal import parse_date_to_decimal, parse_dates_from_auspice_json, extract_date_from_string
from aeon_core.io import ensure_parent_directory, write_json, write_csv
```

No `try/except` guards needed — `aeon-core` is a hard dependency that includes the model, weights, and inference utilities.

In `autoclock.py`, change:
```python
# OLD
from .dataset import (
    compute_tn93_distance_matrix,
    compute_tn93_cross_distance_matrix,
    parse_alignment_sequences,
)
from .dating import (
    run_mrca_dating,
    parse_sample_dates,
    verify_coding_alignment,
)
# Inside a conditional block:
from .inference import load_model, prepare_alignment, get_device
from .splits import extract_cross_taxa_attentions_and_embeddings
from .dating import compute_neural_covariance_kernel

# NEW
from aeon_core.dataset import (
    compute_tn93_distance_matrix,
    compute_tn93_cross_distance_matrix,
    parse_alignment_sequences,
)
from aeon_core.inference import load_model, prepare_alignment, get_device
from aeon_core.splits import extract_cross_taxa_attentions_and_embeddings
from .dating import (
    run_mrca_dating,
    parse_sample_dates,
    verify_coding_alignment,
    compute_neural_covariance_kernel,
)
```

In `triage.py` (was `sieve.py`), change:
```python
# OLD
from .dataset import parse_alignment_sequences
from .temporal import parse_date_to_decimal, extract_date_from_string
from .dating import run_ols_dating, compute_fieller_mrca_interval

# NEW
from aeon_core.dataset import parse_alignment_sequences
from aeon_core.temporal import parse_date_to_decimal, extract_date_from_string
from .dating import run_ols_dating, compute_fieller_mrca_interval
```

In `geo.py`, change:
```python
# OLD
from .dataset import compute_tn93_distance_matrix, parse_alignment_sequences
# Inside a conditional block:
from .inference import get_device, load_model, prepare_alignment
from .splits import extract_cross_taxa_attentions_and_embeddings

# NEW
from aeon_core.dataset import compute_tn93_distance_matrix, parse_alignment_sequences
from aeon_core.inference import get_device, load_model, prepare_alignment
from aeon_core.splits import extract_cross_taxa_attentions_and_embeddings
```

`r0.py`, `sketch.py`, `alignment.py` — no internal imports to change (they import nothing from the package).

**2.5. Create `chronaeon/src/chronaeon/__init__.py`:**

```python
"""
ChronAeon: Ultra-Fast Molecular Clock Dating, Phylodynamics, and Genomic Surveillance.
"""

from .dating import (
    run_mrca_dating,
    run_ols_dating,
    run_pgls_dating,
    run_restricted_spline_clock_dating,
    run_powerlaw_clock_dating,
    verify_coding_alignment,
    parse_sample_dates,
    generate_consensus_sequence,
    generate_time_decay_consensus_sequence,
)
from .geo import (
    run_phylogeography_analysis,
    estimate_spatial_pgls_epicenter,
    parse_geo_metadata,
)
from .r0 import (
    run_r0_analysis,
    compute_reproduction_numbers,
    plot_r0_diagnostics,
    PATHOGEN_PRESETS,
)
from .autoclock import (
    AutoClockDeconvolution,
    run_autoclock_deconvolution,
    HierarchicalAutoClock,
    run_hierarchical_autoclock,
)
from .triage import ChronAeonSieve
from .sketch import (
    CanonicalMinHashSketcher,
    AlignmentFreeBinner,
    AlignmentFreeCentrifuge,
)
from .alignment import (
    ReferenceCodonAligner,
    ReferenceGuidedCodonThreader,
)

__version__ = "0.1.0"
__all__ = [
    "run_mrca_dating", "run_ols_dating", "run_pgls_dating",
    "run_restricted_spline_clock_dating", "run_powerlaw_clock_dating",
    "verify_coding_alignment", "parse_sample_dates",
    "generate_consensus_sequence", "generate_time_decay_consensus_sequence",
    "run_phylogeography_analysis", "estimate_spatial_pgls_epicenter",
    "parse_geo_metadata",
    "run_r0_analysis", "compute_reproduction_numbers",
    "plot_r0_diagnostics", "PATHOGEN_PRESETS",
    "AutoClockDeconvolution", "run_autoclock_deconvolution",
    "HierarchicalAutoClock", "run_hierarchical_autoclock",
    "ChronAeonSieve",  # class name kept for backward compat; module renamed to triage
    "CanonicalMinHashSketcher", "AlignmentFreeBinner", "AlignmentFreeCentrifuge",
    "ReferenceCodonAligner", "ReferenceGuidedCodonThreader",
]
```

**2.6. Create `chronaeon/src/chronaeon/cli.py`:**

Extract the ChronAeon subcommand handlers (`cmd_date`, `cmd_autoclock`, `cmd_triage`, `cmd_phylogeo`, `cmd_dynamics`, and the `sketch`/`align` handlers) from the current `hyphaeon/cli.py` into a new `chronaeon/cli.py`. Create a new `main()` function with its own `argparse` parser registering only ChronAeon subcommands.

The CLI structure should mirror the current one — copy the argument parser definitions for `date`, `autoclock`, `triage`, `phylogeo`, `dynamics`, `sketch`, `align` subcommands verbatim (renaming from the old `dating`, `sieve`, `geo`, `r0` names), updating import paths:
- `from hyphaeon.dating import ...` → `from chronaeon.dating import ...`
- `from hyphaeon.autoclock import ...` → `from chronaeon.autoclock import ...`
- `from hyphaeon.sieve import ...` → `from chronaeon.triage import ...`
- `from hyphaeon.geo import ...` → `from chronaeon.geo import ...`
- `from hyphaeon.r0 import ...` → `from chronaeon.r0 import ...`
- `from hyphaeon.sketch import ...` → `from chronaeon.sketch import ...`
- `from hyphaeon.alignment import ...` → `from chronaeon.alignment import ...`
- Any `from hyphaeon.inference/model/weights/dataset/stats/io` → `from aeon_core.<module> import ...`

Register `sieve`, `dating`, `geo`, `r0` as backward-compatible aliases on the new primary subcommands (`triage`, `date`, `phylogeo`, `dynamics` respectively).

**2.7. Move ChronAeon tests:**

Move these test files from `tests/` to `chronaeon/tests/`:
- `test_dating.py`
- `test_geo.py`
- `test_r0.py`
- `test_triage.py` (was `test_sieve.py`)
- `test_sketch_and_alignment.py`

Update imports in these test files from `hyphaeon.dating` → `chronaeon.dating`, etc. Also update any `from hyphaeon.dataset/inference/splits/temporal/io/stats` → `from aeon_core.<module>`.

Create `chronaeon/tests/conftest.py` — likely minimal (no torch needed for most tests). Check if any ChronAeon tests import from `hyphaeon` and update accordingly.

Also add new test files (see Step 1b):
- `chronaeon/tests/test_imports.py` — import smoke test
- `chronaeon/tests/test_autoclock_imports.py` — exercises conditional import block in `autoclock.py`

**2.8. Move ChronAeon scripts:**

Move these from `scripts/` to `chronaeon/scripts/`:
- `audit_empirical_sus_genomes.py`
- `harvest_sars2_bvbrc.py`
- `run_sars2_streaming_triage.py` (was `run_sars2_streaming_sieve.py`)
- `stream_10k_sars2_bvbrc.py`

Update imports in these scripts:
- `from hyphaeon.sieve import ...` → `from chronaeon.triage import ...`
- `from hyphaeon.dating import ...` → `from chronaeon.dating import ...`
- `from hyphaeon.temporal import ...` → `from aeon_core.temporal import ...`
- `from hyphaeon.dataset/inference/model/weights/stats/io` → `from aeon_core.<module> import ...`

**2.9. Move ChronAeon docs:**

Move these to `chronaeon/`:
- `DATING_GUIDE.md` → `chronaeon/DATING_GUIDE.md`
- `AUTOCLOCK_GUIDE.md` → `chronaeon/AUTOCLOCK_GUIDE.md`
- `CHRONAEON_TRIAGE_PRODUCTION_DESIGN.md` → `chronaeon/CHRONAEON_TRIAGE_PRODUCTION_DESIGN.md`
- `MRCA_DATING_REPORT.md` → `chronaeon/MRCA_DATING_REPORT.md`

---

### Step 3: Update `hyphaeon/` package

**3.1. Create new `hyphaeon/pyproject.toml`:**

```toml
[build-system]
requires = ["setuptools>=61.0"]
build-backend = "setuptools.build_meta"

[project]
name = "hyphaeon"
version = "0.2.0"
authors = [
  { name="Sergei L. Kosakovsky Pond", email="spond@temple.edu" },
]
description = "HyphAeon: A Deep-Time Phylogenetic Foundation Model for Multi-Scale Evolutionary, Structural, and Clinical Genomics"
readme = "README.md"
requires-python = ">=3.8"
license = { text = "MIT" }
classifiers = [
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]
dependencies = [
    "aeon-core>=0.1.0",
    "torch>=2.0.0",
    "biopython>=1.80",
    "numpy>=1.22.0",
    "scipy>=1.9.0",
    "pandas>=1.5.0",
    "networkx>=3.0",
    "matplotlib>=3.5",
]

[project.optional-dependencies]
tn93 = ["tn93>=1.2.0"]
dev = ["pytest>=7.0"]
model_eval = ["pytest>=7.0", "scikit-learn>=1.0"]
all = ["tn93>=1.2.0", "pytest>=7.0", "scikit-learn>=1.0"]

[project.scripts]
hyphaeon = "hyphaeon.cli:main"

[tool.setuptools]
packages = ["hyphaeon", "hyphaeon.model_eval", "hyphaeon.model_eval.calibration", "hyphaeon.model_eval.concordance", "hyphaeon.model_eval.invariance", "hyphaeon.model_eval.stability"]

[tool.setuptools.package-dir]
hyphaeon = "src/hyphaeon"
"hyphaeon.model_eval" = "src/hyphaeon/model_eval"

[tool.pytest.ini_options]
testpaths = ["tests"]

[project.urls]
"Homepage" = "https://github.com/veg/HyphAeon"
"Bug Tracker" = "https://github.com/veg/HyphAeon/issues"
```

**3.2. Move hyphaeon source to `hyphaeon/src/hyphaeon/`:**

Move these files from the old `hyphaeon/` directory to `hyphaeon/src/hyphaeon/`:
- `cli.py`, `phenotype.py`, `epistasis.py`, `disease.py`, `filter.py`, `evaluation.py`, `attribution.py`, `training_data.py`

Create a slimmed `inference.py` containing only `predict_site_lrts`. It imports from `aeon_core`:
```python
from aeon_core.inference import get_device, load_model, prepare_alignment, compute_adaptive_safe_batch_size
from aeon_core.model import PhyloAxialTransformer
```

Create a slimmed `splits.py` containing only `spectral_bisection`, `tree_dict_to_newick`, `get_all_clade_taxa`, `run_spectral_splits`. It imports from `aeon_core`:
```python
from aeon_core.splits import extract_cross_taxa_attentions_and_embeddings, compute_fused_affinity_matrix
from aeon_core.inference import load_model, get_device
from aeon_core.dataset import load_alignment_and_tree
```

Create a slimmed `temporal.py` with only the 4 selection-dynamics functions, importing date utilities from `aeon_core.temporal`.

Move `model_eval/` to `hyphaeon/src/hyphaeon/model_eval/`.

**3.3. Update imports in remaining hyphaeon files:**

In every file that previously did `from .dataset import X`, `from .model import X`, `from .inference import X`, `from .weights import X`, `from .splits import X`, `from .stats import X`, or `from .io import X` where those are now in `aeon_core`, change to `from aeon_core.<module> import X`.

Specifically:
- `inference.py` (slimmed): `from aeon_core.inference import get_device, load_model, prepare_alignment, compute_adaptive_safe_batch_size`; `from aeon_core.model import PhyloAxialTransformer`; `from aeon_core._progress import ChunkProgress`
- `splits.py` (slimmed): `from aeon_core.splits import extract_cross_taxa_attentions_and_embeddings, compute_fused_affinity_matrix`; `from aeon_core.inference import load_model, get_device`; `from aeon_core.dataset import load_alignment_and_tree`
- `phenotype.py`: `from aeon_core.dataset import ...`; `from aeon_core.model import PhyloAxialTransformer`; `from aeon_core.stats import cauchy_combination_p, benjamini_hochberg`; `from aeon_core.inference import get_device, load_model`; `from aeon_core.weights import ...`; `from .epistasis import compute_transformer_attributions`
- `epistasis.py`: `from aeon_core.dataset import ...`; `from aeon_core.model import PhyloAxialTransformer`; `from aeon_core.weights import load_weights, load_arch_config`; `from aeon_core.stats import pvals_from_lrt_self_liang, benjamini_hochberg`; `from aeon_core.inference import get_device, load_model, get_device_memory_budget, compute_adaptive_safe_batch_size`; `from aeon_core._progress import ChunkProgress`
- `disease.py`: `from aeon_core.dataset import ...`; `from aeon_core.model import PhyloAxialTransformer`; `from aeon_core.inference import get_device, load_model, compute_adaptive_safe_batch_size`; `from aeon_core._progress import ChunkProgress`
- `filter.py`: `from aeon_core.dataset import ...`; `from aeon_core.model import PhyloAxialTransformer`; `from aeon_core.weights import resolve_weights_path, DEFAULT_VARIANT`; `from aeon_core.inference import get_device, load_model, compute_adaptive_safe_batch_size`; `from aeon_core.stats import pvals_from_lrt_meme, benjamini_hochberg, cauchy_combination_p` (note: `compute_adaptive_safe_batch_size` was previously imported from `.epistasis` — now import from `aeon_core.inference`)
- `evaluation.py`: `from aeon_core.io import ensure_parent_directory`
- `attribution.py`: `from aeon_core.dataset import CODON_TO_AA, AA_MAP, GENETIC_CODE`; `from aeon_core._progress import ChunkProgress`
- `training_data.py`: `from aeon_core.dataset import load_alignment_and_tree` (currently `from hyphaeon.dataset import load_alignment_and_tree` — change to `aeon_core`)
- `cli.py`: `from aeon_core.inference import get_device, load_model, prepare_alignment, compute_adaptive_safe_batch_size`; `from aeon_core.model import PhyloAxialTransformer, BustedMultiTaskHead`; `from aeon_core.weights import ...`; `from aeon_core.dataset import load_alignment_and_tree`; `from aeon_core.stats import ...`; `from aeon_core.io import ...`; `from aeon_core._progress import ChunkProgress`; `from .phenotype import ...`; `from .epistasis import ...`; `from .temporal import run_temporal_surveillance`; `from .disease import ...`; `from .filter import ...`; `from .splits import run_spectral_splits`; `from .evaluation import ...`; `from .attribution import ...`

**3.4. Update `hyphaeon/__init__.py`:**

Remove all ChronAeon imports and re-export model/weights/dataset from `aeon_core` for backward compatibility:
```python
# Remove these imports entirely (modules moved to chronaeon):
from .dating import ...
from .geo import ...
from .r0 import ...
from .autoclock import ...
from .sketch import ...
from .alignment import ...

# Change these from .model/.dataset to aeon_core re-exports:
# OLD: from .model import PhyloAxialTransformer
# NEW:
from aeon_core.model import PhyloAxialTransformer
# OLD: from .dataset import load_alignment_and_tree, parse_alignment_sequences, compute_tn93_distance_matrix, parse_beast_xml
# NEW:
from aeon_core.dataset import (
    load_alignment_and_tree, parse_alignment_sequences,
    compute_tn93_distance_matrix, parse_beast_xml,
)

# These imports stay (modules remain in hyphaeon):
from .phenotype import run_phenotype_association, resolve_phenotype_vector, PRESETS
from .epistasis import (
    run_epistasis_analysis, run_digital_dms_analysis, run_dms_analysis,
    compute_phylogenetic_branch_attributions, compute_branch_coselection_network,
    compute_selection_dms_essm, extract_epistatic_sectors,
)
from .disease import predict_disease_pathogenicity

# Update __all__ to remove all ChronAeon names:
# Remove: run_mrca_dating, run_ols_dating, run_pgls_dating,
#   run_restricted_spline_clock_dating, run_powerlaw_clock_dating,
#   verify_coding_alignment, parse_sample_dates,
#   generate_consensus_sequence, generate_time_decay_consensus_sequence,
#   run_phylogeography_analysis, estimate_spatial_pgls_epicenter,
#   parse_geo_metadata, run_r0_analysis, compute_reproduction_numbers,
#   plot_r0_diagnostics, PATHOGEN_PRESETS, AutoClockDeconvolution,
#   run_autoclock_deconvolution, HierarchicalAutoClock,
#   run_hierarchical_autoclock, CanonicalMinHashSketcher,
#   AlignmentFreeBinner, AlignmentFreeCentrifuge,
#   ReferenceCodonAligner, ReferenceGuidedCodonThreader
```

**3.5. Remove ChronAeon subcommands from `hyphaeon/cli.py`:**

Remove `cmd_date`, `cmd_autoclock`, `cmd_triage`, `cmd_phylogeo`, `cmd_dynamics`, and their `add_parser` calls (including any alias registrations for `dating`, `sieve`, `geo`, `r0`). Remove `sketch` and `align` subcommands.

**3.6. Move hyphaeon tests to `hyphaeon/tests/`:**

Move all non-ChronAeon test files from `tests/` to `hyphaeon/tests/` (see Section 2.4 for the exact classification). Update `conftest.py` to import from `aeon_core.model` instead of `hyphaeon.model`. Also add `test_imports.py` (new pre-migration import smoke test — see Step 1b).

**3.7. Update `model_eval/` imports:**

In `model_eval/_harness.py`, change:
```python
# OLD
from hyphaeon import dataset as ds
from hyphaeon.stats import pvals_from_lrt_meme as pvals_from_lrt
from hyphaeon.attribution import attribute_selection
from hyphaeon.inference import predict_site_lrts

# NEW
from aeon_core import dataset as ds
from aeon_core.stats import pvals_from_lrt_meme as pvals_from_lrt
from hyphaeon.attribution import attribute_selection  # stays in hyphaeon
from hyphaeon.inference import predict_site_lrts  # stays in hyphaeon
```

In `model_eval/conftest.py`, change:
```python
# OLD
from hyphaeon.model import PhyloAxialTransformer
from hyphaeon import dataset as ds
from hyphaeon.inference import load_model
from hyphaeon.weights import resolve_weights_path

# NEW
from aeon_core.model import PhyloAxialTransformer
from aeon_core import dataset as ds
from aeon_core.inference import load_model
from aeon_core.weights import resolve_weights_path
```

Also update `sys.path` manipulation in `model_eval/conftest.py` — `REPO_ROOT` calculation stays the same (still the repo root), but the package is now at `hyphaeon/src/hyphaeon/` not `hyphaeon/`. With `pip install -e hyphaeon/` this `sys.path` insert is unnecessary and can be removed.

---

### Step 4: Update root-level files

**4.1. Update root `README.md`:**

Rewrite as a monorepo overview describing all three packages with install instructions for each. Structure: brief description of the `aeon-core` foundation model, then separate sections for HyphAeon ("The Depth of Evolution" — site-axis selection inference) and ChronAeon ("The Velocity of Time" — taxon/lineage-axis phylodynamics and outbreak surveillance). Include the "Radar & Microscope" collaborative framing (see Section 1.1).

**4.2. Update `train.py`:**

Change `from hyphaeon.model import PhyloAxialTransformer, decode_soft_ordinal_lrt` to `from aeon_core.model import PhyloAxialTransformer, decode_soft_ordinal_lrt`. The `from hyphaeon.training_data import GeneTensorsDataset` import stays — `training_data.py` remains in the `hyphaeon` package.

**4.3. Move ChronAeon docs:**

Move `DATING_GUIDE.md`, `AUTOCLOCK_GUIDE.md`, `CHRONAEON_TRIAGE_PRODUCTION_DESIGN.md`, `MRCA_DATING_REPORT.md` into `chronaeon/`.

**4.3a. Create independent package READMEs:**

- `chronaeon/README.md` — Tailored specifically to molecular clock dating, phylodynamics, and outbreak surveillance. Target audience: public health agencies, outbreak epidemiologists, hospital infection control teams. Include: ChronAeon positioning ("The Velocity of Time"), CLI subcommand overview (`date`, `autoclock`, `triage`, `phylogeo`, `dynamics`, `sketch`, `align`), install instructions (`pip install chronaeon`), and example workflows for tMRCA dating and R₀ estimation. Do **not** include selection/MEME tutorials — those belong in `hyphaeon/` docs.
- `hyphaeon/README.md` — Tailored to selection inference, epistasis, and phenotype attribution. Target audience: comparative evolutionary biologists, structural virologists. Include: HyphAeon positioning ("The Depth of Evolution"), CLI subcommand overview, install instructions (`pip install hyphaeon`), and example workflows for MEME/BUSTED and epistasis analysis.
- `aeon-core/README.md` — Brief technical README describing the shared foundation model (`PhyloAxialTransformer`), weight resolution from HuggingFace, and the shared utilities (dataset, inference, splits, temporal, io, stats). Install: `pip install aeon-core`.

**4.4. Move ChronAeon example data:**

Move ChronAeon-specific examples to `chronaeon/examples/`:
- `examples/H1N1_2009_pandemic.fasta`, `examples/H1N1_2009_pandemic.nwk`
- `examples/H5N1_HA_geo.fasta`, `examples/H5N1_HA.nwk`, `examples/H5N1_HA_metadata.csv`
- `examples/korber_env_gp160.fasta`

Keep shared examples (HIV1_RT, bat_oas1, camelid, Smc6, RHO) at root `examples/`.

**4.5. Move benchmark data:**

Move `benchmarks/` directory to `hyphaeon/benchmarks/` (contains HyphAeon benchmark data like `benchmark_splits_vs_ml.csv`).

**4.6. Move ChronAeon result images:**

Move these root-level images to `chronaeon/`:
- `mrca_dating_benchmark.png`, `mrca_dating_diagnostic.pdf`, `mrca_dating_diagnostic.png`
- `ebola_r0_example.png`, `h1n1_r0_example.png`, `h5n1_geo_example.png`

**4.7. Update script imports:**

HyphAeon scripts (stay in `scripts/`):
- `build_training_npz.py`: `from hyphaeon.training_data import build_training_directory` — stays (training_data is in hyphaeon)
- `run_avian_genome_filter.py`: `from hyphaeon.weights import ...` → `from aeon_core.weights import ...`; `from hyphaeon.model import ...` → `from aeon_core.model import ...`; `from hyphaeon.dataset import ...` → `from aeon_core.dataset import ...`
- `run_orthomam_masked_selection.py`: `from hyphaeon.model import ...` → `from aeon_core.model import ...`; `from hyphaeon.weights import ...` → `from aeon_core.weights import ...`; `from hyphaeon.dataset import ...` → `from aeon_core.dataset import ...`; `from hyphaeon.stats import ...` → `from aeon_core.stats import ...`

ChronAeon scripts (move to `chronaeon/scripts/`):
- `run_sars2_streaming_triage.py` (was `run_sars2_streaming_sieve.py`): `from hyphaeon.sieve import ChronAeonSieve` → `from chronaeon.triage import ChronAeonSieve`
- `stream_10k_sars2_bvbrc.py`: `from hyphaeon.sieve import ChronAeonSieve` → `from chronaeon.triage import ChronAeonSieve`; `from hyphaeon.temporal import ...` → `from aeon_core.temporal import ...`

---

### Step 5: Update CI

**5.1. Replace `.github/workflows/tests.yml`:**

Create a matrix that tests all three packages:

```yaml
name: tests

on:
  push:
    branches: [main]
  pull_request:
  workflow_dispatch:

jobs:
  test-aeon-core:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ["3.8", "3.9", "3.10", "3.11", "3.12"]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}
      - run: pip install torch --index-url https://download.pytorch.org/whl/cpu
      - run: pip install -e aeon-core/[dev,tn93]
      - run: python -m pytest aeon-core/tests/ -v

  test-hyphaeon:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ["3.8", "3.9", "3.10", "3.11", "3.12"]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}
      - run: pip install torch --index-url https://download.pytorch.org/whl/cpu
      - run: pip install -e aeon-core/ -e hyphaeon/[dev,tn93]
      - run: python -m pytest hyphaeon/tests/ -v

  test-chronaeon:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ["3.8", "3.9", "3.10", "3.11", "3.12"]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}
      - run: pip install torch --index-url https://download.pytorch.org/whl/cpu
      - run: pip install -e aeon-core/ -e chronaeon/[dev,tn93]
      - run: python -m pytest chronaeon/tests/ -v
```

**5.2. Add publish workflow `.github/workflows/publish.yml`:**

```yaml
name: publish

on:
  push:
    tags:
      - "aeon-core-v*"
      - "hyphaeon-v*"
      - "chronaeon-v*"

jobs:
  publish:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.11"
      - name: Determine package from tag
        id: pkg
        run: |
          TAG=${GITHUB_REF#refs/tags/}
          if [[ "$TAG" == aeon-core-v* ]]; then
            echo "dir=aeon-core" >> $GITHUB_OUTPUT
          elif [[ "$TAG" == hyphaeon-v* ]]; then
            echo "dir=hyphaeon" >> $GITHUB_OUTPUT
          elif [[ "$TAG" == chronaeon-v* ]]; then
            echo "dir=chronaeon" >> $GITHUB_OUTPUT
          fi
      - name: Build and publish
        run: |
          pip install build twine
          cd ${{ steps.pkg.outputs.dir }}
          python -m build
          twine upload dist/* -u __token__ -p ${{ secrets.PYPI_API_TOKEN }}
```

---

### Step 6: Verify

**6.1. Local smoke test:**
```bash
# Install all three in editable mode
pip install -e aeon-core/ -e hyphaeon/ -e chronaeon/

# Test aeon-core (needs torch)
pip install torch --index-url https://download.pytorch.org/whl/cpu
python -m pytest aeon-core/tests/ -v

# Test hyphaeon
python -m pytest hyphaeon/tests/ -v

# Test chronaeon
python -m pytest chronaeon/tests/ -v

# Test CLIs
hyphaeon --help
chronaeon --help

# Test cross-package import
python -c "from aeon_core.dataset import parse_alignment_sequences; print('ok')"
python -c "from aeon_core.model import PhyloAxialTransformer; print('model ok')"
python -c "from aeon_core.inference import load_model; print('inference ok')"
python -c "from chronaeon.dating import run_mrca_dating; print('dating ok')"
python -c "from chronaeon.dating import run_pgls_dating; print('pgls ok')"  # uses aeon_core model
```

**6.2. Check for stale references:**
```bash
# Search for any remaining `from hyphaeon.dating` or `from hyphaeon.geo` etc.
grep -r "from hyphaeon\.\(dating\|autoclock\|sieve\|triage\|geo\|r0\|sketch\|alignment\)" --include="*.py" .
# Should return nothing

# Search for any remaining `from .dataset import` in chronaeon or hyphaeon
grep -r "from \.\(dataset\|model\|weights\|stats\|io\|_progress\) import" chronaeon/ hyphaeon/ --include="*.py"
# Should return nothing (should all be from aeon_core)

# Search for remaining `from .inference import` or `from .splits import` in chronaeon
# (these are valid in hyphaeon/ since slimmed versions exist there)
grep -r "from \.\(inference\|splits\|temporal\) import" chronaeon/ --include="*.py"
# Should return nothing (should all be from aeon_core)

# Search for remaining `from hyphaeon.model` or `from hyphaeon.weights` etc. in scripts
grep -r "from hyphaeon\.\(model\|weights\|dataset\|stats\|io\)" scripts/ --include="*.py"
# Should return nothing (should all be from aeon_core)
```

---

## 5. Key Decisions & Rationale

### 5.1 Why the model, weights, and inference live in `aeon-core`

The `PhyloAxialTransformer` model, HuggingFace weights management (`weights.py`), and model utilities (`load_model`, `prepare_alignment`, `get_device`) are used by **both** packages:

- **HyphAeon** uses them for site-level selection inference (MEME, BUSTED, phenotype, epistasis, disease, filter).
- **ChronAeon** uses them for attention-PGLS dating, latent manifold dating, and autoclock deconvolution.

Since the model is intended to be used with pretrained weights in both packages, it is a **hard dependency** for both — not an optional extra. Placing it in `aeon-core` means `pip install chronaeon` automatically gets the model + weights, and `pip install hyphaeon` does too. No `[neural]` extra is needed.

The trade-off is that `aeon-core` (and therefore both packages) requires PyTorch. This is acceptable because the model is central to both packages' intended use.

### 5.2 Why `inference.py` and `splits.py` are split

`inference.py` contains two categories of functions:
- **Generic model utilities** (`get_device`, `load_model`, `prepare_alignment`, `compute_adaptive_safe_batch_size`, `get_device_memory_budget`) — used by both packages → `aeon_core/inference.py`
- **Selection-specific logic** (`predict_site_lrts`) — runs the site-level LRT prediction loop, only used by HyphAeon's CLI → `hyphaeon/inference.py`

`splits.py` similarly contains:
- **Model-attention utilities** (`extract_cross_taxa_attentions_and_embeddings`, `compute_fused_affinity_matrix`) — used by both packages → `aeon_core/splits.py`
- **Spectral bisection CLI logic** (`spectral_bisection`, `run_spectral_splits`, `tree_dict_to_newick`, `get_all_clade_taxa`) — HyphAeon-specific → `hyphaeon/splits.py`

### 5.3 Why `dataset.py` moves entirely to `aeon-core`

`load_alignment_and_tree()` prepares PyTorch tensors for the transformer model. Since the model itself is in `aeon-core`, the tensor preparation belongs there too. There is no reason to split it out — both packages use `load_alignment_and_tree` when they need model-ready tensors.

### 5.4 Why `stats.py` goes to `aeon-core` despite selection-specific names

The functions `pvals_from_lrt_meme` and `pvals_from_lrt_self_liang` have selection-specific names but are generic statistical utilities (chi-square survival function, Liang self-adjusted p-values). Any future `*aeon` tool doing LRT-based inference would use them. The names can be kept as-is — they describe the statistical method, not the biological application.

### 5.5 Why `temporal.py` is split

The 3 date-parsing functions are generic (used by ChronAeon for sampling dates, by HyphAeon for temporal surveillance). The 4 selection-dynamics functions (`parse_temporal_metadata`, `infer_root_sequence`, `run_temporal_surveillance`, `render_temporal_4panel_figure`) are HyphAeon-specific. Clean split.

### 5.6 Why `src/` layout per package

Using `src/hyphaeon/`, `src/chronaeon/`, `src/aeon_core/` within each package directory prevents accidental imports from the wrong location during development. With `pip install -e`, the package is importable by its installed name, not by its directory path.

---

## 6. File-by-File Import Change Reference

### `dating.py`
| Old import | New import |
| :--- | :--- |
| `from .dataset import parse_alignment_sequences` | `from aeon_core.dataset import parse_alignment_sequences` |
| `from .dataset import compute_tn93_distance_matrix` | `from aeon_core.dataset import compute_tn93_distance_matrix` |
| `from .dataset import compute_tn93_cross_distance_matrix` | `from aeon_core.dataset import compute_tn93_cross_distance_matrix` |
| `from .dataset import parse_beast_xml` | `from aeon_core.dataset import parse_beast_xml` |
| `from .dataset import load_alignment_and_tree` | `from aeon_core.dataset import load_alignment_and_tree` |
| `from .dataset import GENETIC_CODE, CODON_TO_AA` | `from aeon_core.dataset import GENETIC_CODE, CODON_TO_AA` |
| `from .inference import load_model, get_device, prepare_alignment` | `from aeon_core.inference import load_model, get_device, prepare_alignment` |
| `from .splits import extract_cross_taxa_attentions_and_embeddings, compute_fused_affinity_matrix` | `from aeon_core.splits import extract_cross_taxa_attentions_and_embeddings, compute_fused_affinity_matrix` |
| `from .temporal import parse_date_to_decimal, ...` | `from aeon_core.temporal import parse_date_to_decimal, ...` |
| `from .io import ensure_parent_directory, write_json, write_csv` | `from aeon_core.io import ensure_parent_directory, write_json, write_csv` |

### `autoclock.py`
| Old import | New import |
| :--- | :--- |
| `from .dataset import compute_tn93_distance_matrix, compute_tn93_cross_distance_matrix, parse_alignment_sequences` | `from aeon_core.dataset import ...` |
| `from .inference import load_model, prepare_alignment, get_device` (guarded) | `from aeon_core.inference import load_model, prepare_alignment, get_device` (direct, no guard) |
| `from .splits import extract_cross_taxa_attentions_and_embeddings` (guarded) | `from aeon_core.splits import extract_cross_taxa_attentions_and_embeddings` (direct, no guard) |
| `from .dating import run_mrca_dating, parse_sample_dates, verify_coding_alignment` | `from .dating import ...` (stays — `dating` is in `chronaeon`) |
| `from .dating import compute_neural_covariance_kernel` (guarded) | `from .dating import compute_neural_covariance_kernel` (stays — `dating` is in `chronaeon`) |

### `sieve.py`
| Old import | New import |
| :--- | :--- |
| `from .dataset import parse_alignment_sequences` | `from aeon_core.dataset import parse_alignment_sequences` |
| `from .temporal import parse_date_to_decimal, extract_date_from_string` | `from aeon_core.temporal import parse_date_to_decimal, extract_date_from_string` |
| `from .dating import run_ols_dating, compute_fieller_mrca_interval` | `from .dating import ...` (stays) |

### `geo.py`
| Old import | New import |
| :--- | :--- |
| `from .dataset import compute_tn93_distance_matrix, parse_alignment_sequences` | `from aeon_core.dataset import ...` |
| `from .inference import get_device, load_model, prepare_alignment` (guarded) | `from aeon_core.inference import get_device, load_model, prepare_alignment` (direct, no guard) |
| `from .splits import extract_cross_taxa_attentions_and_embeddings` (guarded) | `from aeon_core.splits import extract_cross_taxa_attentions_and_embeddings` (direct, no guard) |

### `r0.py`, `sketch.py`, `alignment.py`
No internal imports — no changes needed.

### `hyphaeon/cli.py`
| Old import | New import |
| :--- | :--- |
| `from .model import PhyloAxialTransformer, BustedMultiTaskHead` | `from aeon_core.model import PhyloAxialTransformer, BustedMultiTaskHead` |
| `from .weights import ...` | `from aeon_core.weights import ...` |
| `from .inference import get_device, load_model, prepare_alignment, predict_site_lrts, compute_adaptive_safe_batch_size` | `from aeon_core.inference import get_device, load_model, prepare_alignment, compute_adaptive_safe_batch_size`; `from .inference import predict_site_lrts` |
| `from .dataset import load_alignment_and_tree` | `from aeon_core.dataset import load_alignment_and_tree` |
| `from .stats import pvals_from_lrt_meme, ...` | `from aeon_core.stats import ...` |
| `from .io import ensure_parent_directory, ...` | `from aeon_core.io import ...` |
| `from ._progress import ChunkProgress` | `from aeon_core._progress import ChunkProgress` |
| `from .temporal import run_temporal_surveillance` | `from .temporal import run_temporal_surveillance` (stays — temporal dynamics stay in hyphaeon) |
| `from .phenotype import run_phenotype_association, PRESETS` | `from .phenotype import ...` (stays) |
| `from .epistasis import run_epistasis_analysis, ...` | `from .epistasis import ...` (stays) |
| `from .disease import predict_disease_pathogenicity` | `from .disease import ...` (stays) |
| `from .filter import run_alignment_filter` | `from .filter import ...` (stays) |
| `from .splits import run_spectral_splits` | `from .splits import ...` (stays) |
| `from .evaluation import command as evaluate_command` | `from .evaluation import ...` (stays) |
| `from .attribution import attribute_selection` | `from .attribution import ...` (stays) |
| `from .training_data import build_training_directory` | `from .training_data import ...` (stays) |

### `hyphaeon/inference.py` (slimmed)
| Old import | New import |
| :--- | :--- |
| `from .model import PhyloAxialTransformer` | `from aeon_core.model import PhyloAxialTransformer` |
| `from .weights import resolve_weights_path, load_arch_config, load_weights` | `from aeon_core.weights import resolve_weights_path, load_arch_config, load_weights` |
| `from .dataset import load_alignment_and_tree` | `from aeon_core.dataset import load_alignment_and_tree` |
| `from ._progress import ChunkProgress` | `from aeon_core._progress import ChunkProgress` |

Note: `hyphaeon/inference.py` only contains `predict_site_lrts`. All other functions moved to `aeon_core/inference.py`.

### `hyphaeon/temporal.py` (slimmed)
| Old import | New import |
| :--- | :--- |
| `from .inference import load_model, get_device, prepare_alignment, compute_adaptive_safe_batch_size` | `from aeon_core.inference import load_model, get_device, prepare_alignment, compute_adaptive_safe_batch_size` |
| `from .stats import pvals_from_lrt_self_liang, benjamini_hochberg` | `from aeon_core.stats import pvals_from_lrt_self_liang, benjamini_hochberg` |
| `from .io import ensure_parent_directory, write_json, write_csv` | `from aeon_core.io import ensure_parent_directory, write_json, write_csv` |
| `from .epistasis import REV_AA_MAP, AA_MAP` | `from .epistasis import REV_AA_MAP, AA_MAP` (stays) |

Note: `hyphaeon/temporal.py` only contains the 4 selection-dynamics functions. The 3 date-parsing functions moved to `aeon_core/temporal.py`.

### `hyphaeon/model_eval/_harness.py`
| Old import | New import |
| :--- | :--- |
| `from hyphaeon import dataset as ds` | `from aeon_core import dataset as ds` |
| `from hyphaeon.stats import pvals_from_lrt_meme` | `from aeon_core.stats import pvals_from_lrt_meme` |
| `from hyphaeon.attribution import attribute_selection` | `from hyphaeon.attribution import attribute_selection` (stays) |
| `from hyphaeon.inference import predict_site_lrts` | `from hyphaeon.inference import predict_site_lrts` (stays) |

### `hyphaeon/model_eval/conftest.py`
| Old import | New import |
| :--- | :--- |
| `from hyphaeon.model import PhyloAxialTransformer` | `from aeon_core.model import PhyloAxialTransformer` |
| `from hyphaeon import dataset as ds` | `from aeon_core import dataset as ds` |
| `from hyphaeon.inference import load_model` | `from aeon_core.inference import load_model` |
| `from hyphaeon.weights import resolve_weights_path` | `from aeon_core.weights import resolve_weights_path` |

### `train.py` (root)
| Old import | New import |
| :--- | :--- |
| `from hyphaeon.model import PhyloAxialTransformer, decode_soft_ordinal_lrt` | `from aeon_core.model import PhyloAxialTransformer, decode_soft_ordinal_lrt` |
| `from hyphaeon.training_data import GeneTensorsDataset` | `from hyphaeon.training_data import GeneTensorsDataset` (stays) |

### `scripts/build_training_npz.py`
| Old import | New import |
| :--- | :--- |
| `from hyphaeon.training_data import build_training_directory` | `from hyphaeon.training_data import build_training_directory` (stays) |

### `scripts/run_avian_genome_filter.py`
| Old import | New import |
| :--- | :--- |
| `from hyphaeon.weights import load_arch_config, load_weights` | `from aeon_core.weights import load_arch_config, load_weights` |
| `from hyphaeon.model import PhyloAxialTransformer` | `from aeon_core.model import PhyloAxialTransformer` |
| `from hyphaeon.dataset import load_alignment_and_tree, CODON_TO_AA, parse_alignment_sequences` | `from aeon_core.dataset import load_alignment_and_tree, CODON_TO_AA, parse_alignment_sequences` |

### `scripts/run_orthomam_masked_selection.py`
| Old import | New import |
| :--- | :--- |
| `from hyphaeon.model import PhyloAxialTransformer` | `from aeon_core.model import PhyloAxialTransformer` |
| `from hyphaeon.weights import load_arch_config, load_weights` | `from aeon_core.weights import load_arch_config, load_weights` |
| `from hyphaeon.dataset import load_alignment_and_tree` | `from aeon_core.dataset import load_alignment_and_tree` |
| `from hyphaeon.stats import pvals_from_lrt_meme, benjamini_hochberg, cauchy_combination_p` | `from aeon_core.stats import pvals_from_lrt_meme, benjamini_hochberg, cauchy_combination_p` |

### `scripts/run_sars2_streaming_triage.py` (was `run_sars2_streaming_sieve.py`, moves to `chronaeon/scripts/`)
| Old import | New import |
| :--- | :--- |
| `from hyphaeon.sieve import ChronAeonSieve` | `from chronaeon.triage import ChronAeonSieve` |

### `scripts/stream_10k_sars2_bvbrc.py` (moves to `chronaeon/scripts/`)
| Old import | New import |
| :--- | :--- |
| `from hyphaeon.sieve import ChronAeonSieve` | `from chronaeon.triage import ChronAeonSieve` |
| `from hyphaeon.temporal import parse_date_to_decimal, extract_date_from_string` | `from aeon_core.temporal import parse_date_to_decimal, extract_date_from_string` |

---

## 7. Release Order

Dependencies must be released in order:

1. **`aeon-core v0.1.0`** — no dependencies on other packages (includes model, weights, inference, dataset, splits, temporal, io, stats)
2. **`hyphaeon v0.2.0`** — depends on `aeon-core>=0.1.0`
3. **`chronaeon v0.1.0`** — depends on `aeon-core>=0.1.0`

Tag format: `aeon-core-v0.1.0`, `hyphaeon-v0.2.0`, `chronaeon-v0.1.0`

---

## 8. Post-Refactor Cleanup

- Delete old `hyphaeon/` directory at repo root (replaced by `hyphaeon/src/hyphaeon/`)
- Delete old `tests/` directory at repo root (split into per-package `tests/`)
- Delete root `pyproject.toml` (replaced by per-package configs)
- Update `.gitignore` for new directory structure
- Update `CITATION.cff` to define independent citation metadata for all three packages. The ChronAeon citation should reference the ChronAeon phylodynamics/outbreak surveillance paper, and the HyphAeon citation should reference the HyphAeon foundation model/selection inference paper. `aeon-core` should be cited as the shared foundation model. Use `preferred-citation` with multiple `references` entries in the CFF format.
- Update `docs/` site to cover all three packages with distinct landing pages for ChronAeon (phylodynamics/outbreak surveillance) and HyphAeon (selection inference), linked by a shared `aeon-core` foundation page
- Add `aeon-core/README.md`, `chronaeon/README.md` (hyphaeon README stays at root or moves to `hyphaeon/`)
- Update `scripts/` — move ChronAeon scripts to `chronaeon/scripts/`, update imports in remaining HyphAeon scripts
