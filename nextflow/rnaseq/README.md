# RNA-seq Nextflow pipeline

This pipeline implements the RNA-seq workflow already described in this repository:

1. ENA/local FASTQ input
2. raw FastQC + MultiQC
3. Trimmomatic PE trimming
4. trimmed FastQC + MultiQC
5. STAR genome index generation or reuse
6. STAR alignment
7. HTSeq gene counts
8. Cufflinks FPKM quantification
9. matrix aggregation with `scripts/generate_matrices.py`

## Directory

```text
nextflow/rnaseq/
├── assets/
│   ├── gse153005_rnaseq_samplesheet.tsv
│   └── local_samplesheet.example.tsv
├── main.nf
└── nextflow.config
```

## Required software

Keep the same toolchain you already use in the markdown workflow:

- `nextflow`
- `aria2c`
- `fastqc`
- `multiqc`
- `trimmomatic`
- `STAR`
- `htseq-count`
- `cufflinks`
- `python3`

## Required references

Default paths assume you launch from the repository root and keep references in `ref/`:

- `ref/GRCh38.primary_assembly.genome.fa`
- `ref/gencode.v29.annotation.gtf`
- `ref/TruSeq3-PE.fa`
- optional prebuilt STAR index at `ref/star_index/`

If `ref/star_index/` is missing or empty, the pipeline will build it.
If `ref/GRCh38.primary_assembly.genome.fa`, `ref/gencode.v29.annotation.gtf`, or
`ref/TruSeq3-PE.fa` are missing, the pipeline will download them automatically by default.

## Inputs

The samplesheet must be tab-delimited and contain:

- `sample`
- `fastq_1`
- `fastq_2`

Optional columns:

- `input_type`: `url` or `local`
- `day`
- `replicate`

If `input_type` is omitted, the pipeline infers it from `fastq_1`:

- `http://`, `https://`, `ftp://` -> `url`
- everything else -> `local`

## Run

Run from the repository root:

```bash
nextflow run nextflow/rnaseq \
  -profile standard \
  --outdir results/rnaseq_nf
```

The default samplesheet is already set to:

```text
nextflow/rnaseq/assets/gse153005_rnaseq_samplesheet.tsv
```

To use your own local FASTQ files:

```bash
nextflow run nextflow/rnaseq \
  -profile standard \
  --samplesheet nextflow/rnaseq/assets/local_samplesheet.example.tsv \
  --outdir results/rnaseq_local
```

To run on SLURM:

```bash
nextflow run nextflow/rnaseq \
  -profile slurm \
  -c nextflow/rnaseq/conf/slurm.example.config \
  --samplesheet nextflow/rnaseq/assets/gse153005_rnaseq_samplesheet.tsv \
  --outdir results/rnaseq_nf
```

The example cluster overlay is here:

```text
nextflow/rnaseq/conf/slurm.example.config
```

To run with Docker from a host that already has Nextflow:

```bash
docker build -t melendez-rnaseq:latest docker/rnaseq

nextflow run nextflow/rnaseq \
  -profile docker \
  --container_image melendez-rnaseq:latest \
  --outdir results/rnaseq_nf_docker \
  -resume
```

To run with the all-in-one container only, without installing Nextflow on the
host:

```bash
docker build -t melendez-rnaseq:latest docker/rnaseq

docker run --rm \
  -v "$PWD":/workspace \
  -w /workspace \
  melendez-rnaseq:latest run \
  --samplesheet nextflow/rnaseq/assets/test_one_sample.tsv \
  --outdir results/rnaseq_nf_test \
  -resume
```

The container also exposes helper subcommands:

```bash
docker run --rm -v "$PWD":/workspace -w /workspace melendez-rnaseq:latest help
docker run --rm -v "$PWD":/workspace -w /workspace melendez-rnaseq:latest config
docker run --rm -v "$PWD":/workspace -w /workspace melendez-rnaseq:latest nextflow -version
```

On Apple Silicon hosts, build and run the image as `linux/amd64` because
`cufflinks` is not available from Bioconda for `linux-aarch64`:

```bash
docker build --platform linux/amd64 -t melendez-rnaseq:latest docker/rnaseq

docker run --rm --platform linux/amd64 \
  -v "$PWD":/workspace \
  -w /workspace \
  melendez-rnaseq:latest run \
  --samplesheet nextflow/rnaseq/assets/test_one_sample.tsv \
  --outdir results/rnaseq_nf_test \
  -resume
```

If your host does not have Java 17+ for Nextflow, you can validate the pipeline
configuration with the official Nextflow container instead of changing the host
JDK:

```bash
docker run --rm --platform linux/amd64 \
  -v "$PWD":/workspace \
  -w /workspace \
  nextflow/nextflow:25.10.4 \
  nextflow config nextflow/rnaseq -profile docker
```

## Main parameters

```text
--samplesheet
--outdir
--container_image
--reference_dir
--auto_download_refs
--run_multiqc
--aria2c_cmd
--fastqc_cmd
--multiqc_cmd
--trimmomatic_cmd
--star_cmd
--htseq_cmd
--cufflinks_cmd
--genome_fasta
--genome_fasta_url
--gtf
--gtf_url
--star_index
--read_length
--sjdb_overhang
--publish_mode
--trimmomatic_java_heap
--trimmomatic_adapter_fa
--trimmomatic_adapter_url
--trimmomatic_clip_args
--trimmomatic_trim_args
--star_sort_ram
--star_extra
--htseq_extra
--cufflinks_extra
--python_bin
```

## Output layout

```text
results/rnaseq_nf/
├── counts/
├── fpkm_cufflinks/
├── mapped_data/
├── pipeline_info/
├── qc_raw/
├── qc_trimmed/
├── raw_data/
├── ref/
├── trimmed_data/
├── GSE153005_count_mat.txt
└── GSE153005_FPKM_mat.txt
```

## Notes

- `sample` names are used as FASTQ/download prefixes and BAM/count output prefixes.
- The bundled GSE153005 samplesheet uses `SRR` accessions as sample IDs so that `generate_matrices.py` can rename them back to `d0_r1 ... d28_r3`.
- `publish_mode` defaults to `symlink` to avoid unnecessary duplication of BAM/FASTQ files. Switch to `copy` if you want a self-contained results directory.
- The pipeline keeps the same RNA quantification logic as your existing markdown workflow: `STAR -> HTSeq counts + Cufflinks FPKM -> generate_matrices.py`.
- Missing references are auto-downloaded into `reference_dir` by default.
- Set `--auto_download_refs false` if you want strict local-only behavior.
- Set `--run_multiqc false` if your `multiqc` environment is incompatible and you want to skip only the report aggregation step.
- Tool commands are parameterized. This lets you split environments on bare metal, for example:

```bash
nextflow run nextflow/rnaseq \
  -profile standard \
  --multiqc_cmd "conda run -n multiqc_env multiqc" \
  --cufflinks_cmd "conda run -n cuff cufflinks" \
  --outdir results/rnaseq_nf \
  -resume
```

- A Docker image recipe is provided in:

```text
docker/rnaseq/
```

- The Docker image contains:
  - `nextflow=25.10.4`
- The Docker profile assumes one image with two conda environments:
  - `rnaseq`: `aria2`, `fastqc`, `multiqc`, `trimmomatic`, `STAR`, `htseq`, `python`
  - `cufflinks`: `cufflinks=2.2.1`
