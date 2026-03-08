# Melendez-2021-SciRep-mDA-Multiomics

Reproducible analysis code and workflow notes for the study:

Melendez-Ramirez et al. 2021, *Dynamic landscape of chromatin accessibility and transcriptomic changes during differentiation of human embryonic stem cells into dopaminergic neurons*  
Scientific Reports  
DOI: [10.1038/s41598-021-96263-1](https://doi.org/10.1038/s41598-021-96263-1)

This repository is organized around practical reuse:

- reproduce Figure 2 from RNA-seq matrices
- reproduce Figure 3 from RNA-seq matrices
- build count / FPKM matrices from HTSeq and Cufflinks outputs
- run the RNA-seq processing workflow from FASTQ with Nextflow
- run the RNA-seq workflow through Docker without installing Nextflow on the host
- keep detailed paper and workflow notes in a separate methods document

## Overview

This repository currently focuses on the RNA-seq side of the paper and on figure-level reproduction from processed matrices.

What is already implemented:

- RNA-seq matrix generation from `counts/` and `fpkm_cufflinks/`
- Figure 2 reproduction with DESeq2-based DEG calling
- Figure 3 reproduction from the same matrices
- RNA-seq Nextflow workflow
- Docker image recipe for the RNA-seq workflow

What is not yet a full standalone pipeline here:

- a complete ATAC-seq automation pipeline comparable to the RNA-seq Nextflow workflow

ATAC workflow notes and caveats are still documented, but the repository is currently strongest for RNA-seq and figure reproduction.

## What This Repository Contains

| Task | Main entry point | Typical use case |
|---|---|---|
| Reproduce Figure 2 | `scripts/reproduce_fig2.R` | You already have count and FPKM matrices |
| Reproduce Figure 3 | `scripts/reproduce_fig3.R` | You already have count and FPKM matrices |
| Build matrices | `scripts/generate_matrices.py` | You have per-sample HTSeq and Cufflinks outputs |
| Run RNA-seq from FASTQ | `nextflow/rnaseq/` | You want a full processing workflow |
| Run RNA-seq with Docker only | `docker/rnaseq/` | You do not want to install Nextflow locally |
| Read detailed workflow notes | `Melendez-2021-SciRep-mDA-Multiomics.md` | You want the paper-oriented methods walkthrough |

## Repository Layout

```text
.
├── README.md
├── Melendez-2021-SciRep-mDA-Multiomics.md
├── scripts/
│   ├── generate_matrices.py
│   ├── reproduce_fig2.R
│   └── reproduce_fig3.R
├── nextflow/
│   └── rnaseq/
│       ├── README.md
│       ├── main.nf
│       ├── nextflow.config
│       ├── assets/
│       └── conf/
├── docker/
│   └── rnaseq/
│       ├── Dockerfile
│       ├── entrypoint.sh
│       └── envs/
└── results/
```

Directory roles:

- `scripts/`: figure reproduction and matrix construction
- `nextflow/rnaseq/`: production RNA-seq workflow
- `docker/rnaseq/`: containerized RNA-seq runtime
- `results/`: generated outputs
- `Melendez-2021-SciRep-mDA-Multiomics.md`: detailed paper and workflow notes

## Quick Start

Choose the path that matches your current input.

### 1. You already have `GSE153005_count_mat.txt` and `GSE153005_FPKM_mat.txt`

This is the fastest path if your goal is to reproduce the paper figures.

Create a `data/` directory and place the two matrices there:

```text
data/
├── GSE153005_count_mat.txt
└── GSE153005_FPKM_mat.txt
```

Then run either script from the repository root.

For Figure 2:

```bash
Rscript scripts/reproduce_fig2.R
```

For Figure 3:

```bash
Rscript scripts/reproduce_fig3.R
```

These scripts are written to be RStudio-friendly. If you prefer interactive use, open the script, edit the parameter block at the top, and click `Run All`.

If your matrices are not stored under `data/`, edit the top variables in the script before running:

- `count_file`
- `fpkm_file`
- `outdir`

### 2. You have per-sample HTSeq counts and Cufflinks FPKM outputs

If you do not yet have the two merged matrices, build them first.

Expected input layout:

```text
ref/
└── gencode.v29.annotation.gtf

counts/
├── SRR12070831_counts.txt
├── SRR12070832_counts.txt
└── ...

fpkm_cufflinks/
├── SRR12070831/
│   └── genes.fpkm_tracking
├── SRR12070832/
│   └── genes.fpkm_tracking
└── ...
```

Run:

```bash
python3 scripts/generate_matrices.py
```

This writes:

```text
GSE153005_count_mat.txt
GSE153005_FPKM_mat.txt
```

Then move or copy those files into `data/` and run the figure scripts.

### 3. You want to start from raw RNA-seq FASTQ files

Use the Nextflow workflow in `nextflow/rnaseq/`.

Start here:

- [nextflow/rnaseq/README.md](./nextflow/rnaseq/README.md)

Typical run:

```bash
nextflow run nextflow/rnaseq \
  -profile standard \
  --outdir results/rnaseq_nf \
  -resume
```

The workflow can auto-download the reference FASTA, GTF, and Trimmomatic adapter FASTA if they are missing.

### 4. You want to run the RNA-seq workflow with Docker only

If you do not want to install Nextflow on the host, build the provided all-in-one image and run the pipeline from inside the container.

Build:

```bash
docker build -t melendez-rnaseq:latest docker/rnaseq
```

Run:

```bash
docker run --rm \
  -v "$PWD":/workspace \
  -w /workspace \
  melendez-rnaseq:latest run \
  --samplesheet nextflow/rnaseq/assets/gse153005_rnaseq_samplesheet.tsv \
  --outdir results/rnaseq_nf \
  -resume
```

Useful helper commands:

```bash
docker run --rm -v "$PWD":/workspace -w /workspace melendez-rnaseq:latest help
docker run --rm -v "$PWD":/workspace -w /workspace melendez-rnaseq:latest config
docker run --rm melendez-rnaseq:latest nextflow -version
```

On Apple Silicon hosts, use `linux/amd64` because `cufflinks` is not available from Bioconda for `linux-aarch64`:

```bash
docker build --platform linux/amd64 -t melendez-rnaseq:latest docker/rnaseq

docker run --rm --platform linux/amd64 \
  -v "$PWD":/workspace \
  -w /workspace \
  melendez-rnaseq:latest run \
  --samplesheet nextflow/rnaseq/assets/gse153005_rnaseq_samplesheet.tsv \
  --outdir results/rnaseq_nf \
  -resume
```

## Data Sources

Primary public sources:

- GEO series: `GSE153005`
- BioProject: `PRJNA641151`
- SRA study: `SRP268344`

RNA-seq samples released in the study:

- D0: 3 replicates
- D14: 4 replicates
- D28: 3 replicates

ATAC-seq samples released in the study:

- D0: 1 library
- D14: 1 library
- D28: 1 library

Important implication:

- RNA-seq can be reproduced much more directly from public data
- ATAC-seq public release is more limited in replicate structure, so some paper-level ATAC statistics cannot be reproduced under exactly the same evidence structure as the original internal analysis

## Reproducing Figure 2

Script:

- `scripts/reproduce_fig2.R`

R packages used by the script:

- `DESeq2`
- `dplyr`
- `readr`
- `tidyr`
- `tibble`
- `stringr`
- `ggplot2`
- `VennDiagram`
- `pheatmap`
- `grid`

Default inputs:

- `data/GSE153005_count_mat.txt`
- `data/GSE153005_FPKM_mat.txt`

Default output directory:

- `results/fig2/`

Main outputs:

- `results/fig2/Fig2A_venn.png`
- `results/fig2/Fig2B_up_down_barplot.png`
- `results/fig2/Fig2C_heatmap.png`
- `results/fig2/Fig2C_right_panel_own_enrichment.png`
- `results/fig2/Fig2_summary_counts.tsv`
- `results/fig2/run_parameters.tsv`

Method assumptions in the current script:

- DESeq2 pairwise comparisons:
  - `D0 vs D14`
  - `D14 vs D28`
  - `D0 vs D28`
- DEG thresholds:
  - `FDR < 0.05`
  - `FC > 4`
  - at least one sample with `FPKM >= 1`
- DESeq2 settings:
  - `fitType = "local"`
  - `sfType = "poscounts"`
  - `independentFiltering = TRUE`
  - `cooksCutoff = TRUE`

The script is a single entry point for:

- Fig2A Venn diagram
- Fig2B up/down protein-coding vs lncRNA barplot
- Fig2C heatmap
- Fig2C right-side enrichment panel

## Reproducing Figure 3

Script:

- `scripts/reproduce_fig3.R`

R packages used by the script:

- `DESeq2`
- `dplyr`
- `readr`
- `tidyr`
- `tibble`
- `stringr`
- `ggplot2`
- `forcats`

Default inputs:

- `data/GSE153005_count_mat.txt`
- `data/GSE153005_FPKM_mat.txt`

Default output directory:

- `results/fig3/`

Main outputs:

- `results/fig3/Fig3A_top_degs.png`
- `results/fig3/Fig3B_lnc_classes.png`
- `results/fig3/Fig3C_scatter_combined.png`
- `results/fig3/Fig3A_top_counts.tsv`
- `results/fig3/Fig3B_lnc_class_counts.tsv`
- `results/fig3/Fig3C_pairs_all.tsv`

Current Figure 3 implementation uses the same DEG backbone as Figure 2 and applies a stricter paper-style interpretation for Fig3A:

- top-ranked DEGs are ranked by absolute `log2FoldChange`
- the hypergeometric enrichment background for Fig3A defaults to comparison-specific DEGs

## Building Matrices from Sample-Level Outputs

Script:

- `scripts/generate_matrices.py`

What it does:

- parses `gencode.v29.annotation.gtf`
- keeps protein-coding genes and paper-relevant lncRNA classes
- merges HTSeq count tables
- merges Cufflinks `genes.fpkm_tracking`
- collapses pseudoautosomal duplicate gene IDs
- renames `SRR` columns to paper-style sample names such as `d0_r1`, `d14_r4`, `d28_r3`

Default required paths:

- `ref/gencode.v29.annotation.gtf`
- `counts/*_counts.txt`
- `fpkm_cufflinks/*/genes.fpkm_tracking`

Run:

```bash
python3 scripts/generate_matrices.py
```

Outputs:

- `GSE153005_count_mat.txt`
- `GSE153005_FPKM_mat.txt`

## RNA-seq Workflow

The full FASTQ-to-matrix workflow lives in:

- [nextflow/rnaseq/README.md](./nextflow/rnaseq/README.md)

In brief, the workflow performs:

1. FASTQ input from ENA URLs or local files
2. raw FastQC and MultiQC
3. Trimmomatic trimming
4. trimmed FastQC and MultiQC
5. STAR genome indexing or reuse
6. STAR alignment
7. HTSeq gene counting
8. Cufflinks FPKM estimation
9. matrix aggregation with `scripts/generate_matrices.py`

Use the Nextflow documentation if:

- you want reproducible batch execution
- you want automatic reference download
- you want one command from FASTQ to matrices
- you want Dockerized execution

## Workflow Documentation

Detailed paper-aligned notes remain in:

- [Melendez-2021-SciRep-mDA-Multiomics.md](./Melendez-2021-SciRep-mDA-Multiomics.md)

Use that document when you need:

- step-by-step paper-oriented RNA-seq notes
- detailed ATAC-seq notes and caveats
- accession-level context from GEO / SRA
- rationale behind parameter choices

Use this top-level `README.md` when you need:

- to understand what is in the repository
- to decide which entry point to run
- to get started quickly without reading all workflow notes

## Known Limitations

This repository is practical and reproducible, but not every discrepancy with the paper can be eliminated.

### RNA-seq

- exact DEG counts are sensitive to the upstream processing chain
- figure-level results depend strongly on ranking and filtering definitions
- the paper does not expose every implementation detail behind the published DEG tables

### ATAC-seq

- public ATAC raw data are limited to one released library per time point
- this constrains direct reproduction of paper-level differential accessibility statistics
- ATAC workflow notes are documented, but the repository is not yet centered on a full ATAC automation pipeline

### Docker

- on Apple Silicon, use `--platform linux/amd64`
- the first image build may take a long time because it solves and downloads large Conda environments
- occasional package download timeouts can happen during the first build; rerunning the build is usually sufficient

## Citation

If you use this repository, please cite the original paper:

```text
Melendez-Ramirez LY, et al. Dynamic landscape of chromatin accessibility and transcriptomic changes during differentiation of human embryonic stem cells into dopaminergic neurons. Scientific Reports. 2021.
```

Primary paper DOI:

- [10.1038/s41598-021-96263-1](https://doi.org/10.1038/s41598-021-96263-1)

Primary data accession:

- `GSE153005`
