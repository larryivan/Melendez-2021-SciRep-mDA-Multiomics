nextflow.enable.dsl=2

params.samplesheet = params.samplesheet ?: "${projectDir}/assets/gse153005_rnaseq_samplesheet.tsv"
params.outdir = params.outdir ?: "${launchDir}/results/rnaseq_nf"
params.publish_mode = params.publish_mode ?: 'symlink'
params.reference_dir = params.reference_dir ?: "${launchDir}/ref"
params.auto_download_refs = params.auto_download_refs == null ? true : params.auto_download_refs.toString().toBoolean()

params.genome_fasta = params.genome_fasta ?: "${params.reference_dir}/GRCh38.primary_assembly.genome.fa"
params.gtf = params.gtf ?: "${params.reference_dir}/gencode.v29.annotation.gtf"
params.star_index = params.star_index ?: "${params.reference_dir}/star_index"
params.read_length = (params.read_length ?: 75) as Integer
params.sjdb_overhang = (params.sjdb_overhang ?: (params.read_length - 1)) as Integer
params.genome_fasta_url = params.genome_fasta_url ?: 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz'
params.gtf_url = params.gtf_url ?: 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz'

params.download_connections = params.download_connections ?: 8
params.download_splits = params.download_splits ?: 8
params.download_min_split = params.download_min_split ?: '50M'
params.run_multiqc = params.run_multiqc == null ? true : params.run_multiqc.toString().toBoolean()
params.aria2c_cmd = params.aria2c_cmd ?: 'aria2c'
params.fastqc_cmd = params.fastqc_cmd ?: 'fastqc'
params.multiqc_cmd = params.multiqc_cmd ?: 'multiqc'
params.trimmomatic_cmd = params.trimmomatic_cmd ?: 'trimmomatic'
params.star_cmd = params.star_cmd ?: 'STAR'
params.htseq_cmd = params.htseq_cmd ?: 'htseq-count'
params.cufflinks_cmd = params.cufflinks_cmd ?: 'cufflinks'

params.trimmomatic_java_heap = params.trimmomatic_java_heap ?: '40g'
params.trimmomatic_adapter_fa = params.trimmomatic_adapter_fa ?: "${params.reference_dir}/TruSeq3-PE.fa"
params.trimmomatic_adapter_url = params.trimmomatic_adapter_url ?: 'https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE.fa'
params.trimmomatic_clip_args = params.trimmomatic_clip_args ?: '2:30:10'
params.trimmomatic_trim_args = params.trimmomatic_trim_args ?: 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

params.star_sort_ram = (params.star_sort_ram ?: 100000000000) as Long
params.star_extra = params.star_extra ?: ''
params.htseq_extra = params.htseq_extra ?: ''
params.cufflinks_extra = params.cufflinks_extra ?: ''
params.python_bin = params.python_bin ?: 'python3'

def inferInputType(String fastq1, String explicitType) {
    if (explicitType) {
        return explicitType.trim().toLowerCase()
    }
    return fastq1 ==~ /^(https?|ftp):\/\/.*/ ? 'url' : 'local'
}

def requireReadableFile(String label, String pathString) {
    def f = new File(pathString)
    if (!f.exists()) {
        throw new IllegalArgumentException("${label} not found: ${pathString}")
    }
    if (!f.canRead()) {
        throw new IllegalArgumentException("${label} is not readable: ${pathString}")
    }
}

def readableFileExists(String pathString) {
    def f = new File(pathString)
    return f.exists() && f.canRead()
}

def useExistingStarIndex(String pathString) {
    def d = new File(pathString)
    return d.exists() && d.isDirectory() && (d.list()?.size() ?: 0) > 0
}

requireReadableFile('Samplesheet', params.samplesheet.toString())
if (!params.auto_download_refs && !readableFileExists(params.gtf.toString())) {
    requireReadableFile('GTF', params.gtf.toString())
}
if (!params.auto_download_refs && !readableFileExists(params.trimmomatic_adapter_fa.toString())) {
    requireReadableFile('Trimmomatic adapter FASTA', params.trimmomatic_adapter_fa.toString())
}
if (!useExistingStarIndex(params.star_index.toString()) && !params.auto_download_refs) {
    requireReadableFile('Genome FASTA', params.genome_fasta.toString())
}

process DOWNLOAD_FASTQ {
    tag "${sample}"
    label 'medium'

    publishDir "${params.outdir}/raw_data", mode: params.publish_mode

    input:
    tuple val(sample), val(fastq1_url), val(fastq2_url)

    output:
    tuple val(sample), path("${sample}_1.fastq.gz"), path("${sample}_2.fastq.gz"), emit: reads

    script:
    """
    ${params.aria2c_cmd} \
      --allow-overwrite=false \
      --auto-file-renaming=false \
      --continue=true \
      --file-allocation=none \
      --max-connection-per-server=${params.download_connections} \
      --split=${params.download_splits} \
      --min-split-size=${params.download_min_split} \
      -o ${sample}_1.fastq.gz \
      "${fastq1_url}"

    ${params.aria2c_cmd} \
      --allow-overwrite=false \
      --auto-file-renaming=false \
      --continue=true \
      --file-allocation=none \
      --max-connection-per-server=${params.download_connections} \
      --split=${params.download_splits} \
      --min-split-size=${params.download_min_split} \
      -o ${sample}_2.fastq.gz \
      "${fastq2_url}"
    """
}

process DOWNLOAD_GENOME_FASTA {
    label 'medium'

    publishDir "${params.reference_dir}", mode: 'copy'

    input:
    val genome_fasta_url

    output:
    path("GRCh38.primary_assembly.genome.fa"), emit: fasta

    script:
    """
    ${params.aria2c_cmd} \
      --allow-overwrite=false \
      --auto-file-renaming=false \
      --continue=true \
      --file-allocation=none \
      --max-connection-per-server=${params.download_connections} \
      --split=${params.download_splits} \
      --min-split-size=${params.download_min_split} \
      -o genome.fa.gz \
      "${genome_fasta_url}"

    gzip -dc genome.fa.gz > GRCh38.primary_assembly.genome.fa
    """
}

process DOWNLOAD_GTF {
    label 'medium'

    publishDir "${params.reference_dir}", mode: 'copy'

    input:
    val gtf_url

    output:
    path("gencode.v29.annotation.gtf"), emit: gtf

    script:
    """
    ${params.aria2c_cmd} \
      --allow-overwrite=false \
      --auto-file-renaming=false \
      --continue=true \
      --file-allocation=none \
      --max-connection-per-server=${params.download_connections} \
      --split=${params.download_splits} \
      --min-split-size=${params.download_min_split} \
      -o gencode.v29.annotation.gtf.gz \
      "${gtf_url}"

    gzip -dc gencode.v29.annotation.gtf.gz > gencode.v29.annotation.gtf
    """
}

process DOWNLOAD_TRIMMOMATIC_ADAPTER {
    label 'light'

    publishDir "${params.reference_dir}", mode: 'copy'

    input:
    val adapter_url

    output:
    path("TruSeq3-PE.fa"), emit: adapter

    script:
    """
    ${params.aria2c_cmd} \
      --allow-overwrite=false \
      --auto-file-renaming=false \
      --continue=true \
      --file-allocation=none \
      --max-connection-per-server=${params.download_connections} \
      --split=1 \
      -o TruSeq3-PE.fa \
      "${adapter_url}"
    """
}

process FASTQC_RAW {
    tag "${sample}"
    label 'light'

    publishDir "${params.outdir}/qc_raw/fastqc", mode: params.publish_mode

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path("*_fastqc.zip"), emit: zips
    path("*_fastqc.html"), emit: html

    script:
    """
    ${params.fastqc_cmd} \
      --threads ${task.cpus} \
      --outdir . \
      ${read1} ${read2}
    """
}

process MULTIQC_RAW {
    label 'light'

    publishDir "${params.outdir}/qc_raw", mode: params.publish_mode

    input:
    path fastqc_zips

    output:
    path("multiqc_raw.html"), emit: html
    path("multiqc_raw_data"), emit: data

    script:
    """
    ${params.multiqc_cmd} \
      --force \
      --filename multiqc_raw.html \
      --outdir . \
      .
    """
}

process TRIMMOMATIC_PE {
    tag "${sample}"
    label 'heavy'

    publishDir "${params.outdir}/trimmed_data", mode: params.publish_mode

    input:
    tuple val(sample), path(read1), path(read2), path(adapter_fa)

    output:
    tuple val(sample), path("${sample}_1_paired.fq.gz"), path("${sample}_2_paired.fq.gz"), emit: paired
    path("${sample}_1_unpaired.fq.gz"), emit: unpaired_r1
    path("${sample}_2_unpaired.fq.gz"), emit: unpaired_r2

    script:
    """
    export _JAVA_OPTIONS="-Xmx${params.trimmomatic_java_heap}"

    ${params.trimmomatic_cmd} PE \
      -threads ${task.cpus} \
      -phred33 \
      ${read1} \
      ${read2} \
      ${sample}_1_paired.fq.gz \
      ${sample}_1_unpaired.fq.gz \
      ${sample}_2_paired.fq.gz \
      ${sample}_2_unpaired.fq.gz \
      ILLUMINACLIP:${adapter_fa}:${params.trimmomatic_clip_args} \
      ${params.trimmomatic_trim_args}
    """
}

process FASTQC_TRIMMED {
    tag "${sample}"
    label 'light'

    publishDir "${params.outdir}/qc_trimmed/fastqc", mode: params.publish_mode

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path("*_fastqc.zip"), emit: zips
    path("*_fastqc.html"), emit: html

    script:
    """
    ${params.fastqc_cmd} \
      --threads ${task.cpus} \
      --outdir . \
      ${read1} ${read2}
    """
}

process MULTIQC_TRIMMED {
    label 'light'

    publishDir "${params.outdir}/qc_trimmed", mode: params.publish_mode

    input:
    path fastqc_zips

    output:
    path("multiqc_trimmed.html"), emit: html
    path("multiqc_trimmed_data"), emit: data

    script:
    """
    ${params.multiqc_cmd} \
      --force \
      --filename multiqc_trimmed.html \
      --outdir . \
      .
    """
}

process STAR_GENOME_INDEX {
    label 'heavy'

    publishDir "${params.outdir}/ref", mode: params.publish_mode

    input:
    path genome_fasta
    path gtf

    output:
    path("star_index"), emit: index

    script:
    """
    mkdir -p star_index

    ${params.star_cmd} \
      --runThreadN ${task.cpus} \
      --runMode genomeGenerate \
      --genomeDir star_index \
      --genomeFastaFiles ${genome_fasta} \
      --sjdbGTFfile ${gtf} \
      --sjdbOverhang ${params.sjdb_overhang}
    """
}

process STAR_ALIGN {
    tag "${sample}"
    label 'heavy'

    publishDir "${params.outdir}/mapped_data", mode: params.publish_mode

    input:
    tuple val(sample), path(read1), path(read2), path(star_index)

    output:
    tuple val(sample), path("${sample}_Aligned.sortedByCoord.out.bam"), emit: bam
    path("${sample}_Log.final.out"), emit: log_final
    path("${sample}_Log.out"), emit: log_out
    path("${sample}_Log.progress.out"), emit: log_progress
    path("${sample}_SJ.out.tab"), emit: sj

    script:
    """
    ${params.star_cmd} \
      --runThreadN ${task.cpus} \
      --genomeDir ${star_index} \
      --readFilesIn ${read1} ${read2} \
      --readFilesCommand zcat \
      --outFileNamePrefix ${sample}_ \
      --outSAMtype BAM SortedByCoordinate \
      --outBAMsortingThreadN ${task.cpus} \
      --limitBAMsortRAM ${params.star_sort_ram} \
      --outSAMstrandField intronMotif \
      ${params.star_extra}
    """
}

process HTSEQ_COUNT {
    tag "${sample}"
    label 'medium'

    publishDir "${params.outdir}/counts", mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(gtf)

    output:
    path("${sample}_counts.txt"), emit: counts

    script:
    """
    ${params.htseq_cmd} \
      -f bam \
      -r pos \
      -s no \
      -m union \
      -t exon \
      -i gene_id \
      ${params.htseq_extra} \
      ${bam} \
      ${gtf} \
      > ${sample}_counts.txt
    """
}

process CUFFLINKS_FPKM {
    tag "${sample}"
    label 'heavy'

    publishDir "${params.outdir}/fpkm_cufflinks", mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(gtf)

    output:
    path("${sample}"), emit: sample_dir

    script:
    """
    mkdir -p ${sample}

    ${params.cufflinks_cmd} \
      -p ${task.cpus} \
      -G ${gtf} \
      -o ${sample} \
      ${params.cufflinks_extra} \
      ${bam}
    """
}

process BUILD_MATRICES {
    label 'medium'

    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path gtf
    path count_files
    path fpkm_dirs
    path matrix_script

    output:
    path("GSE153005_count_mat.txt"), emit: count_matrix
    path("GSE153005_FPKM_mat.txt"), emit: fpkm_matrix

    script:
    """
    mkdir -p ref counts fpkm_cufflinks

    ln -s "${gtf}" ref/gencode.v29.annotation.gtf

    for f in ${count_files}; do
      ln -s "\$f" "counts/\$(basename "\$f")"
    done

    for d in ${fpkm_dirs}; do
      ln -s "\$d" "fpkm_cufflinks/\$(basename "\$d")"
    done

    ln -s "${matrix_script}" generate_matrices.py

    ${params.python_bin} generate_matrices.py
    """
}

workflow {
    Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sample = row.sample?.toString()?.trim()
            def fastq1 = row.fastq_1?.toString()?.trim()
            def fastq2 = row.fastq_2?.toString()?.trim()
            def inputType = inferInputType(fastq1 ?: '', row.input_type?.toString())

            if (!sample) {
                throw new IllegalArgumentException("Each samplesheet row must define 'sample'")
            }
            if (!fastq1 || !fastq2) {
                throw new IllegalArgumentException("Samplesheet row for ${sample} is missing fastq_1 or fastq_2")
            }
            if (!(inputType in ['local', 'url'])) {
                throw new IllegalArgumentException("Unsupported input_type '${inputType}' for sample ${sample}")
            }

            tuple(sample, inputType, fastq1, fastq2)
        }
        .set { samples_ch }

    local_reads_ch = samples_ch
        .filter { sample, inputType, fastq1, fastq2 -> inputType == 'local' }
        .map { sample, inputType, fastq1, fastq2 ->
            tuple(sample, file(fastq1, checkIfExists: true), file(fastq2, checkIfExists: true))
        }

    url_reads_ch = samples_ch
        .filter { sample, inputType, fastq1, fastq2 -> inputType == 'url' }
        .map { sample, inputType, fastq1, fastq2 -> tuple(sample, fastq1, fastq2) }

    downloaded_reads = DOWNLOAD_FASTQ(url_reads_ch)
    raw_reads_ch = local_reads_ch.mix(downloaded_reads.reads)

    raw_fastqc = FASTQC_RAW(raw_reads_ch)
    if (params.run_multiqc) {
        MULTIQC_RAW(raw_fastqc.zips.collect())
    }

    if (readableFileExists(params.trimmomatic_adapter_fa.toString())) {
        adapter_ch = Channel.value(file(params.trimmomatic_adapter_fa, checkIfExists: true))
    } else {
        adapter_dl = DOWNLOAD_TRIMMOMATIC_ADAPTER(Channel.value(params.trimmomatic_adapter_url))
        adapter_ch = adapter_dl.adapter
    }

    trim_input_ch = raw_reads_ch.combine(adapter_ch).map { sample, read1, read2, adapter_fa ->
        tuple(sample, read1, read2, adapter_fa)
    }
    trimmed = TRIMMOMATIC_PE(trim_input_ch)

    trimmed_fastqc = FASTQC_TRIMMED(trimmed.paired)
    if (params.run_multiqc) {
        MULTIQC_TRIMMED(trimmed_fastqc.zips.collect())
    }

    if (readableFileExists(params.gtf.toString())) {
        gtf_ch = Channel.value(file(params.gtf, checkIfExists: true))
    } else {
        gtf_dl = DOWNLOAD_GTF(Channel.value(params.gtf_url))
        gtf_ch = gtf_dl.gtf
    }

    matrix_script_ch = Channel.value(file("${projectDir}/../../scripts/generate_matrices.py", checkIfExists: true))

    if (useExistingStarIndex(params.star_index.toString())) {
        star_index_ch = Channel.value(file(params.star_index, checkIfExists: true))
    } else {
        if (readableFileExists(params.genome_fasta.toString())) {
            fasta_ch = Channel.value(file(params.genome_fasta, checkIfExists: true))
        } else {
            fasta_dl = DOWNLOAD_GENOME_FASTA(Channel.value(params.genome_fasta_url))
            fasta_ch = fasta_dl.fasta
        }
        indexed = STAR_GENOME_INDEX(fasta_ch, gtf_ch)
        star_index_ch = indexed.index
    }

    aligned_input_ch = trimmed.paired.combine(star_index_ch).map { sample, read1, read2, star_index ->
        tuple(sample, read1, read2, star_index)
    }
    aligned = STAR_ALIGN(aligned_input_ch)

    counts_input_ch = aligned.bam.combine(gtf_ch).map { sample, bam, gtf ->
        tuple(sample, bam, gtf)
    }
    fpkm_input_ch = aligned.bam.combine(gtf_ch).map { sample, bam, gtf ->
        tuple(sample, bam, gtf)
    }

    counted = HTSEQ_COUNT(counts_input_ch)
    quantified = CUFFLINKS_FPKM(fpkm_input_ch)

    BUILD_MATRICES(
        gtf_ch,
        counted.counts.collect(),
        quantified.sample_dir.collect(),
        matrix_script_ch
    )
}
