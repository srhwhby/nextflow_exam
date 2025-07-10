#!/usr/bin/env nextflow

params.reads_dir = "/mnt/c/Users/User/Documents/Nextflow/nf-scripts/data/samples/"
params.reads_glob = "*_R{1,2}*.{fq,fastq}" // define how to pick only paired-end reads files with global pattern; optionally add: {.gz,}.
params.out = "results/"
params.accession = "NC_000913.3"

workflow {
    download_ref(params.accession)

    raw_reads_ch = Channel.fromFilePairs("${params.reads_dir}/${params.reads_glob}", flat: true, checkIfExists: true)
    // raw_reads_ch.view() // check channel output object/structure
    fastqc(raw_reads_ch) // [S1, /mnt/c/Users/User/Documents/Nextflow/nf-scripts/data/samples/S1_R1.fq, /mnt/c/Users/User/Documents/Nextflow/nf-scripts/data/samples/S1_R2.fq]
    /*
    fastp(raw_reads_ch)
    bwa_align(fastp.out)
    call_variants(bwa_align.out)
    multiqc(params.out)
    */
}

// sample_ID = 'S1'
// reads = [ '/path/to/S1_R1.fq', '/path/to/S1_R2.fq' ]

process download_ref {

    label 'download_ref'
    publishDir params.out, mode: 'copy'

    input:
        val accession
    
    output:
        path "${accession}.fa"
    
    script:
        """
        if [ ! -f "${accession}.fa" ]; then
            esearch -db nucleotide -query "$accession" | \
            efetch -format fasta > "${accession}.fa"
        else
            echo "${accession}.fa already exists, skipping download"
        fi
        """
}

process fastqc {

    label 'fastqc'
    publishDir params.out, mode: 'copy' 

    input:
    tuple val(sample_ID), path(reads)

    output:
    path "*.zip", emit: zip // output channel zip - only use emit when output is channeled further.
    path "*.html", emit: html // output channel html

    script:
    """
    fastqc ${reads.join(' ')} --outdir .
    mv ${reads[0].getBaseName().split('_')[0]}_R1_fastqc.zip ${sample_ID}.zip
    mv ${reads[0].getBaseName().split('_')[0]}_R1_fastqc.html ${sample_ID}.html
    """
}


// for testing use nextflow run main.nf -resume - this will cache the output and does not hash it to the work dir.
// nextflow log - displays all the runs with names
// nextflow clean -before <run name> -n - deletes every run in the work dir before the named one.
    // -n tells you which hash dir will be removed, so you need to confirm it by
    // nextflow clean -before <run name> -f finally.