# Melendez-2021-SciRep-mDA-Multiomics

## RNA-seq

This guide details the reproducible bioinformatics pipeline for the RNA-Seq data from the study: *"Dynamic landscape of chromatin accessibility and transcriptomic changes during differentiation of human embryonic stem cells into dopaminergic neurons"* (Meléndez-Ramírez et al., Scientific Reports, 2021). The original raw sequencing datasets are deposited under GEO accession **GSE153005**.

### Environment Preparation

To begin with, run all the commands in a `tmux` session is required to maintain the process. And cancel:`ulimit -n 65535`

To ensure a highly reproducible and conflict-free workflow, we use **Mamba** (a faster Conda) to manage all bioinformatics software packages and their dependencies.

Create a custom environment, and install necessary software:

```bash
mamba install multiqc fastqc cutadapt trimmomatic star htseq cufflinks=2.2.1 aria2 -y
```

(Note: The paper explicitly specifies `Cufflinks v2.2.1` for transcript assembly and FPKM calculation.)

### High-Speed Data Downloading (The ENA + aria2c Strategy)

The project includes paired-end 75-cycle RNA-Seq data across three developmental stages (D0, D14, D28) with a total of 10 samples (3, 4, and 3 replicates, respectively).

Instead of using the standard NCBI SRA Toolkit (`prefetch` and `fastq-dump` is too slow for large files), we utilize the **European Nucleotide Archive (ENA) API** combined with **`aria2c`** for accelerated downloading. The corresponding BioProject accession for GSE153005 is **PRJNA641151**.

#### 1. Fetch Metadata and Download Links from ENA

We query the ENA API to retrieve a TSV report containing the direct FTP links for the pre-converted FASTQ files:

```bash
mkdir -p raw_data && cd raw_data

curl -o PRJNA641151_filereport.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA641151&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,library_strategy&format=tsv"

# Extract RNA-Seq links and split paired-end URLs into separate lines
awk -F'\t' '$4 == "RNA-Seq" {print $2}' PRJNA641151_filereport.tsv | tr ';' '\n' > urls.txt

# Replace ftp:// with http:// for enhanced download stability
sed -i 's/ftp.sra.ebi.ac.uk/http:\/\/ftp.sra.ebi.ac.uk/g' urls.txt
```

#### 2. Execute Multi-Threaded Download with aria2c

```bash
aria2c -i urls.txt \
       -d . \
       -c \
       --max-concurrent-downloads=5 \
       -x 16 \
       -s 16 \
       --min-split-size=50M
```

#### 3. Verify the MD5 Value

```bash
# Extract MD5 and filename by using awk
awk -F'\t' '$4 == "RNA-Seq" {
    split($3, md5s, ";"); 
    split($2, urls, ";"); 
    for(i in md5s) {
        split(urls[i], parts, "/");
        print md5s[i] "  " parts[length(parts)]
    }
}' PRJNA641151_filereport.tsv > md5sums.txt
```

```bash
md5sum -c md5sums.txt
```

And make sure all pass(like this):

<img src="/Users/larryivanhan/Library/Application Support/typora-user-images/Screenshot 2026-03-02 at 21.27.38.png" alt="Screenshot 2026-03-02 at 21.27.38" style="zoom:30%;" />

### Quality Control and Adapter Trimming

According to the original study, standard quality procedures were performed using cutadapt, trimmomatic, and FASTQC. Since Trimmomatic integrates highly efficient adapter-clipping capabilities alongside quality filtering, we can achieve optimal cleaning by utilizing a comprehensive Trimmomatic pipeline.

#### 1. Pre-Trimming Quality Assessment (Pre-QC)

Before making any modifications to the raw data, we evaluate the baseline sequencing quality and adapter contamination levels using FastQC and MultiQC.

```bash
mkdir -p ../qc_raw

# Run FastQC using 20 threads for high-speed parallel processing
fastqc -t 20 *.fastq.gz -o ../qc_raw

# Aggregate all individual reports into a single comprehensive MultiQC HTML report
multiqc ../qc_raw -o ../qc_raw
```

#### 2. Trimmomatic

We filter out low-quality bases and Illumina sequencing adapters to prepare the reads for downstream mapping. The original sequencing was performed using a paired-end 75-cycle layout.

**ps:** When running Trimmomatic with multiple threads (`-threads 20`) on a high-performance server, the concurrent gzip compression can easily exhaust the default Java Heap Space, leading to a `java.lang.OutOfMemoryError`. To prevent this, we explicitly allocate 40GB of maximum heap memory to the Java Virtual Machine using the `_JAVA_OPTIONS` environment variable.

```bash
mkdir -p ../trimmed_data

# Break the Java memory limit by explicitly allocating 40GB of heap space
export _JAVA_OPTIONS="-Xmx40g"

# Batch process all paired-end samples in the raw_data directory
for f1 in *_1.fastq.gz; do
    # Extract matching reverse read and sample prefix
    f2=${f1%%_1.fastq.gz}_2.fastq.gz
    sample=${f1%%_1.fastq.gz}
    
    echo "sample: $sample ..."
    
    # Execute Trimmomatic with 20 threads per sample
    trimmomatic PE -threads 20 -phred33 \
        $f1 $f2 \
        ../trimmed_data/${sample}_1_paired.fq.gz ../trimmed_data/${sample}_1_unpaired.fq.gz \
        ../trimmed_data/${sample}_2_paired.fq.gz ../trimmed_data/${sample}_2_unpaired.fq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

echo "Done"
```

#### 3. Post-Trimming Quality Verification

After cleaning, it is essential to re-run FastQC exclusively on the successfully paired reads (`*paired.fq.gz`) to verify that adapter contents have been successfully removed and the per-base sequence quality has improved.

```bash
mkdir -p ../qc_trimmed
cd ../trimmed_data
fastqc -t 20 *_paired.fq.gz -o ../qc_trimmed
multiqc ../qc_trimmed -o ../qc_trimmed
```

### Data Processing

#### Reference Genome Preparation

According to the original study, reads were mapped to the hg38 reference genome. Furthermore, annotated genes and transcripts were derived from GENCODE v29. To ensure maximum reproducibility and fast retrieval, we utilize `aria2c` to download the primary assembly and annotation directly from the European Bioinformatics Institute (EBI) FTP servers.

```bash
mkdir -p ~/mDA_rnaseq/ref && cd ~/mDA_rnaseq/ref

# Download the hg38 genome and GENCODE v29 annotation using aria2c
aria2c -x 16 -s 16 -j 2 -c \
"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz"

aria2c -x 16 -s 16 -j 2 -c \
"ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz"

# Decompress the reference files
gunzip GRCh38.primary_assembly.genome.fa.gz gencode.v29.annotation.gtf.gz
```

#### Genome Indexing (STAR)

Before aligning the reads, a reference genome index must be generated using STAR.  Since the original sequencing was performed using a paired-end 75 cycle mode, the `--sjdbOverhang` parameter is set to 74 (calculated as read length - 1) to construct an optimal splice junction database.

```bash
mkdir -p star_index

STAR --runThreadN 60 \
     --runMode genomeGenerate \
     --genomeDir ./star_index \
     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile gencode.v29.annotation.gtf \
     --sjdbOverhang 74
```

#### Read Alignment

The trimmed paired-end reads are mapped to the hg38 genome.  We leverage STAR's high-performance parallelization and in-memory BAM sorting capabilities to eliminate I/O bottlenecks. Crucially, the `--outSAMstrandField intronMotif` parameter is included. This forces STAR to generate the `XS` attribute for spliced alignments, which is strictly required for downstream compatibility with Cufflinks.

```bash
mkdir -p ~/mDA_rnaseq/mapped_data
cd ~/mDA_rnaseq/trimmed_data

# Alignment

for f1 in *_1_paired.fq.gz; do
    f2=${f1%%_1_paired.fq.gz}_2_paired.fq.gz
    sample=${f1%%_1_paired.fq.gz}
    
    echo "Sample: $sample ..."
    
    STAR --runThreadN 50 \
         --genomeDir ~/mDA_rnaseq/ref/star_index \
         --readFilesIn $f1 $f2 \
         --readFilesCommand zcat \
         --outFileNamePrefix ../mapped_data/${sample}_ \
         --outSAMtype BAM SortedByCoordinate \
         --outBAMsortingThreadN 6 \
         --limitBAMsortRAM 100000000000 \
         --outSAMstrandField intronMotif
done
```

#### Transcriptome Quantification

The original methodology utilizes a dual-quantification approach for downstream analysis.

**1. Raw Read Counts via HTSeq**

We use HTSeq to calculate read counts for annotated genes and transcripts. These raw counts are essential for downstream differential expression analysis utilizing tools like DESeq2. The processes are sent to the background (`&`) to run concurrently, maximizing CPU utilization.

```bash
mkdir -p ~/mDA_rnaseq/counts
cd ~/mDA_rnaseq/mapped_data

# Batch execute HTSeq-count in the background
for bam in *_Aligned.sortedByCoord.out.bam; do
    sample=${bam%%_Aligned.sortedByCoord.out.bam}
    echo "counts: $sample"
    
    htseq-count -f bam -r pos -s no -t exon -i gene_id \
        $bam ~/mDA_rnaseq/ref/gencode.v29.annotation.gtf > ../counts/${sample}_counts.txt &
done

# Wait for all background HTSeq processes to finish
wait
```

**2. FPKM Calculation via Cufflinks**

Assembly of mapped reads was performed with Cufflinks v2.2.1. Furthermore, values of Fragments Per Kilobase of transcript per million Mapped reads (FPKM) were calculated for all annotated genes and transcripts. These normalized values are particularly useful for filtering lowly expressed genes and for visualization purposes.

```bash
cd ~/mDA_rnaseq/mapped_data

# Batch execute Cufflinks FPKM quantification
for bam in *_Aligned.sortedByCoord.out.bam; do
    sample=${bam%%_Aligned.sortedByCoord.out.bam}
    echo "FPKM: $sample"
    
    cufflinks -p 20 -G ~/mDA_rnaseq/ref/gencode.v29.annotation.gtf \
        -o ../fpkm_cufflinks/${sample} $bam
done
```

#### Data Aggregation and Matrix Generation

Finally, we use a custom Python script to parse the GENCODE v29 GTF file, merge individual HTSeq count files, and aggregate Cufflinks FPKM outputs into consolidated matrices. This script handles the assignment of pseudocounts and deduplication of overlapping gene IDs to ensure the final matrices are perfectly formatted for the differential expression pipeline.

```bash
# Generate the final GSE153005_count_mat.txt and GSE153005_FPKM_mat.txt matrices
python generate_matrices.py
```

## ATAC-seq

This is the corrected ATAC-seq workflow aligned to the paper and GEO processing notes for **GSE153005 / PRJNA641151**.

### Critical corrections

1. The public raw ATAC data are only three libraries:
   - `D0_ATACseq  = SRR12070840`
   - `D14_ATACseq = SRR12070841`
   - `D28_ATACseq = SRR12070842`
2. The downstream analysis is:
   - align and filter reads
   - call peaks with MACS2
   - remove ENCODE blacklist peaks
   - merge peaks within `100 bp` to generate consensus peaks
   - count reads in consensus peaks
   - run differential accessibility analysis with `DESeq2`
   - annotate peaks with `HOMER annotatePeaks`
3. The paper/GEO record states MACS2 was run with:
   - `--shift 75 --extsize 150 --nomodel --keep-dup all --call-summits`
4. The public release contains only **one ATAC library per stage**. Therefore:
   - preprocessing, peak calling, blacklist filtering, consensus peak construction, peak annotation and track generation can be reproduced
   - **DESeq2 differential accessibility from public ATAC alone is only exploratory**, because there are no biological replicates deposited for ATAC

### Environment

```bash
mamba create -n mDA_atac \
  python=3.10 \
  fastqc multiqc cutadapt trimmomatic bowtie2 samtools picard macs2 \
  bedtools htseq aria2 homer deeptools -y

mamba activate mDA_atac
```

### Working directory

```bash
mkdir -p ~/mDA_atac/{raw_data,qc_raw,qc_trimmed,trimmed_data,ref,align,peaks_raw,peaks,counts,tracks,metrics}
cd ~/mDA_atac
```

### Reference preparation

Reuse the same hg38 / GRCh38 genome FASTA used for RNA-seq. Also download the ENCODE hg38 blacklist.

```bash
mkdir -p ~/mDA_atac/ref/bowtie2_index ~/mDA_atac/ref/blacklist
cd ~/mDA_atac/ref

# Reuse the genome FASTA from the RNA-seq workflow
ln -sf ~/mDA_rnaseq/ref/GRCh38.primary_assembly.genome.fa .

# Bowtie2 index
bowtie2-build --threads 80 \
  ~/mDA_rnaseq/ref/GRCh38.primary_assembly.genome.fa \
  ~/mDA_atac/ref/bowtie2_index/GRCh38

# ENCODE hg38 blacklist
curl -L https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz \
  -o ~/mDA_atac/ref/blacklist/hg38-blacklist.v2.bed.gz
gunzip -f ~/mDA_atac/ref/blacklist/hg38-blacklist.v2.bed.gz
```

### Download raw ATAC FASTQ files

```bash
cd ~/mDA_atac/raw_data

curl -L \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA641151&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,library_strategy,sample_alias&format=tsv" \
  -o PRJNA641151_atac_filereport.tsv

awk -F'\t' 'NR==1 || $4=="ATAC-seq"' PRJNA641151_atac_filereport.tsv

awk -F'\t' 'NR>1 && $4=="ATAC-seq" {print $2}' PRJNA641151_atac_filereport.tsv | tr ';' '\n' > urls_atac.txt
sed -i 's#ftp.sra.ebi.ac.uk#http://ftp.sra.ebi.ac.uk#g' urls_atac.txt

aria2c -i urls_atac.txt \
  -d . \
  -c \
  --max-concurrent-downloads=3 \
  -x 16 \
  -s 16 \
  --min-split-size=50M

awk -F'\t' 'NR>1 && $4=="ATAC-seq" {
    split($3, md5s, ";");
    split($2, urls, ";");
    for (i in md5s) {
        split(urls[i], parts, "/");
        print md5s[i] "  " parts[length(parts)];
    }
}' PRJNA641151_atac_filereport.tsv > md5sums_atac.txt

md5sum -c md5sums_atac.txt
```

### Raw FASTQ QC

```bash
cd ~/mDA_atac/raw_data

fastqc -t 12 *.fastq.gz -o ../qc_raw
multiqc ../qc_raw -o ../qc_raw
```

### Adapter trimming and quality trimming

The paper states adapters were removed with **cutadapt** and **trimmomatic**.

```bash
cd ~/mDA_atac/raw_data
export _JAVA_OPTIONS="-Xmx40g"

for f1 in *_1.fastq.gz; do
    f2=${f1%%_1.fastq.gz}_2.fastq.gz
    sample=${f1%%_1.fastq.gz}

    cutadapt \
      -j 12 \
      -a CTGTCTCTTATACACATCT \
      -A CTGTCTCTTATACACATCT \
      -o ../trimmed_data/${sample}.cutadapt_R1.fastq.gz \
      -p ../trimmed_data/${sample}.cutadapt_R2.fastq.gz \
      "$f1" "$f2"

    trimmomatic PE -threads 12 -phred33 \
      ../trimmed_data/${sample}.cutadapt_R1.fastq.gz \
      ../trimmed_data/${sample}.cutadapt_R2.fastq.gz \
      ../trimmed_data/${sample}_1_paired.fq.gz \
      ../trimmed_data/${sample}_1_unpaired.fq.gz \
      ../trimmed_data/${sample}_2_paired.fq.gz \
      ../trimmed_data/${sample}_2_unpaired.fq.gz \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
```

### Trimmed FASTQ QC

```bash
cd ~/mDA_atac/trimmed_data

fastqc -t 12 *_paired.fq.gz -o ../qc_trimmed
multiqc ../qc_trimmed -o ../qc_trimmed
```

### Alignment with Bowtie2

The paper uses the ENCODE ATAC-seq pipeline logic and Bowtie2 parameters `-X 2000 -k 5`.

```bash
cd ~/mDA_atac/trimmed_data

for f1 in *_1_paired.fq.gz; do
    f2=${f1%%_1_paired.fq.gz}_2_paired.fq.gz
    sample=${f1%%_1_paired.fq.gz}

    bowtie2 \
      -p 32 \
      -x ~/mDA_atac/ref/bowtie2_index/GRCh38 \
      -1 "$f1" \
      -2 "$f2" \
      -X 2000 \
      -k 5 | \
      samtools sort -@ 16 -m 4G -o ../align/${sample}.coord.bam

    samtools index -@ 16 ../align/${sample}.coord.bam
done
```

### Read filtering and duplicate removal

The paper specifies removing:

- unmapped reads
- fragments with unmapped mates
- non-primary alignments
- reads failing platform QC
- low quality reads (`MAPQ < 30`)
- mitochondrial reads
- optical/PCR duplicates

Use Picard first, then filter the duplicate-marked BAM with `samtools`.

```bash
cd ~/mDA_atac/align
export _JAVA_OPTIONS="-Xmx40g"

for bam in *.coord.bam; do
    sample=${bam%%.coord.bam}

    picard MarkDuplicates \
      I=${sample}.coord.bam \
      O=${sample}.dedup.bam \
      M=../metrics/${sample}.dup_metrics.txt \
      REMOVE_DUPLICATES=true \
      VALIDATION_STRINGENCY=LENIENT

    samtools index -@ 16 ${sample}.dedup.bam

    samtools idxstats ${sample}.dedup.bam | cut -f1 | \
      grep -v -E '^(chrM|MT)$' > ${sample}.keep_chroms.txt

    samtools view -@ 16 -b \
      -f 2 \
      -F 1804 \
      -q 30 \
      ${sample}.dedup.bam \
      $(cat ${sample}.keep_chroms.txt) | \
      samtools sort -@ 16 -m 4G \
      -o ${sample}.filtered.coord.bam \
      -

    samtools index -@ 16 ${sample}.filtered.coord.bam

    test -s ${sample}.filtered.coord.bam || {
      echo "[ERROR] Missing filtered BAM for ${sample}" >&2
      exit 1
    }

    picard CollectInsertSizeMetrics \
      I=${sample}.filtered.coord.bam \
      O=../metrics/${sample}.insert_size_metrics.txt \
      H=../metrics/${sample}.insert_size_histogram.pdf
done
```

### BigWig tracks

The GEO record states bigWig files were generated with `bamCoverage` using `binSize = 50`.

```bash
cd ~/mDA_atac/align

for bam in *.filtered.coord.bam; do
    sample=${bam%%.filtered.coord.bam}

    bamCoverage \
      -b "$bam" \
      -o ../tracks/${sample}.final.bw \
      --binSize 50
done
```

### Peak calling with MACS2

Use the exact paper/GEO MACS2 arguments for the paper-faithful rerun.

```bash
cd ~/mDA_atac/align

for bam in *.filtered.coord.bam; do
    sample=${bam%%.filtered.coord.bam}

    macs2 callpeak \
      -t "$bam" \
      -n "$sample" \
      --outdir ../peaks_raw \
      -f BAM \
      -g hs \
      --shift 75 \
      --extsize 150 \
      --nomodel \
      --keep-dup all \
      --call-summits
done
```

### Remove ENCODE blacklist peaks

```bash
cd ~/mDA_atac

for npk in peaks_raw/*_peaks.narrowPeak; do
    sample=$(basename "$npk" _peaks.narrowPeak)

    bedtools intersect -v \
      -a "$npk" \
      -b ref/blacklist/hg38-blacklist.v2.bed \
      > peaks/${sample}_peaks_filtered.narrowPeak

    bedtools intersect -v \
      -a peaks_raw/${sample}_summits.bed \
      -b ref/blacklist/hg38-blacklist.v2.bed \
      > peaks/${sample}_summits_filtered.bed
done
```

### Build consensus peaks by merging within 100 bp

This matches the paper statement: *"A consensus set of unique peaks was generated merging peaks within 100 bp."*

```bash
cd ~/mDA_atac

cat peaks/*_peaks_filtered.narrowPeak | \
  cut -f1-3 | \
  sort -k1,1 -k2,2n | \
  bedtools merge -i - -d 100 > counts/consensus_peaks.merge100bp.bed

awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"cons_peak_"NR}' \
  counts/consensus_peaks.merge100bp.bed > counts/consensus_peaks.merge100bp.id.bed
```

### Optional: stage-exclusive / shared peaks

If you want stage-specific and shared peaks like the paper discussion, use the filtered peak sets, not the raw MACS2 output.

```bash
cd ~/mDA_atac

bedtools multiinter -i \
  peaks/SRR12070840_peaks_filtered.narrowPeak \
  peaks/SRR12070841_peaks_filtered.narrowPeak \
  peaks/SRR12070842_peaks_filtered.narrowPeak \
  > counts/atac_multiinter.tsv
```

### Convert consensus peaks to GTF for HTSeq counting

The paper states reads in peaks were counted using **HTSeq**.

```bash
cd ~/mDA_atac

awk 'BEGIN{OFS="\t"} {
    print $1, "consensus_peak", "peak", $2+1, $3, ".", "+", ".", \
    "gene_id \"" $4 "\"; transcript_id \"" $4 "\";"
}' counts/consensus_peaks.merge100bp.id.bed > counts/consensus_peaks.merge100bp.gtf
```

### Count reads in consensus peaks

`htseq-count` is safest on name-sorted BAM for paired-end counting.

```bash
cd ~/mDA_atac/align

for bam in *.filtered.coord.bam; do
    sample=${bam%%.filtered.coord.bam}

    samtools sort -n -@ 16 -m 4G -o ../counts/${sample}.name.bam "$bam"

    htseq-count \
      -f bam \
      -r name \
      -s no \
      -m union \
      --secondary-alignments=ignore \
      --supplementary-alignments=ignore \
      -t peak \
      -i gene_id \
      ../counts/${sample}.name.bam \
      ../counts/consensus_peaks.merge100bp.gtf \
      > ../counts/${sample}.counts.txt
done
```

### Merge per-sample peak counts into a matrix

```bash
cd ~/mDA_atac/counts

paste \
  SRR12070840.counts.txt \
  SRR12070841.counts.txt \
  SRR12070842.counts.txt | \
  awk 'BEGIN{OFS="\t"} $1 !~ /^__/ {print $1, $2, $4, $6}' \
  > consensus_peak_count_matrix.tsv

printf "peak_id\tD0\tD14\tD28\n" | cat - consensus_peak_count_matrix.tsv > tmp && mv tmp consensus_peak_count_matrix.tsv
```

### Peak annotation with HOMER

The paper uses `annotatePeaks` with **GENCODE v29** and defines promoter-TSS as `1 kb upstream + 100 bp downstream` of TSS.

```bash
cd ~/mDA_atac

annotatePeaks.pl \
  counts/consensus_peaks.merge100bp.id.bed \
  hg38 \
  -gtf ~/mDA_rnaseq/ref/gencode.v29.annotation.gtf \
  > counts/consensus_peaks.annotated.txt
```

### Differential accessibility analysis

The paper uses:

- `DESeq2`
- `|log2FC| > 1.5`
- `p-value < 0.05`

However, the public ATAC release contains only **one library per stage**. Therefore:

- you can reproduce the **consensus peak count matrix**
- you can reproduce **peak annotations, stage-specific peaks, tracks, and motif inputs**
- you **cannot faithfully reproduce the published ATAC DESeq2 statistics from public data alone**

If you later obtain replicate ATAC BAMs from the authors, then run `DESeq2` on `consensus_peak_count_matrix.tsv`.

### Motif analysis and downstream integration

Once you have:

- `counts/consensus_peaks.merge100bp.id.bed`
- `counts/consensus_peaks.annotated.txt`
- a set of OCRs of interest

you can continue with:

- motif enrichment
- promoter OCR analysis
- RNA/ATAC integration
- IGV track inspection

For example:

```bash
findMotifsGenome.pl \
  peaks_of_interest.bed \
  hg38 \
  motifs_out/ \
  -size given
```
