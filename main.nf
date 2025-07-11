#!/usr/bin/env nextflow

params.reads_raw = "hepatitis"
// params.reads_glob = "*.{fa,fasta}"
params.out = "results"
params.accession_ref = "M21012"
params.combined_fasta = "combined_fasta.fasta"
params.alignment_mafft = "alignment_mafft.fasta"
params.trimming_trimal_report = "trimming_trimal.html"
params.trimming_trimal_fasta = "trimming_trimal.fasta"

process download_ref {

    conda 'bioconda::entrez-direct=24.0'
    publishDir params.reads_raw, mode: 'copy'

    input:
        val accession_ref
    
    output:
        path "${accession_ref}.fa"
    
    script:
        """
        if [ ! -f "${accession_ref}.fa" ]; then
            esearch -db nucleotide -query "$accession_ref" | \
            efetch -format fasta > "${accession_ref}.fa"
        else
            echo "${accession_ref}.fa already exists, skipping download"
        fi
        """
}

process combine_fasta {

    publishDir params.out, mode: 'copy'

    input:
    val reads_raw_channel //path reads_raw_channel //"${params.reads_raw}/${reads_glob}".collect()

    output:
    path "${params.combined_fasta}"

    script:
    """
    cat ${reads_raw_channel.join(' ')} > ${params.combined_fasta}
    """
}

process align_mafft {

    conda 'bioconda::mafft=7.525'
    publishDir params.out, mode: 'copy'

    input:
    path "${params.combined_fasta}" // via combine_fasta.out

    output:
    path "${params.alignment_mafft}"

    script:
    """
    mafft --auto --thread -1 "${params.combined_fasta}" > "${params.alignment_mafft}"
    """

}

process trimming_trimal {

    conda 'bioconda::trimal=1.5.0'
    publishDir params.out, mode: 'copy'

    input:
    path "${params.alignment_mafft}" // via align_mafft.out

    output:
    path "${params.trimming_trimal_report}"
    path "${params.trimming_trimal_fasta}"

    script:
    """
    trimal -in "${params.alignment_mafft}" -out "${params.trimming_trimal_fasta}" -automated1 -htmlout "${params.trimming_trimal_report}" -fasta
    """

}

workflow {
    download_ref_ch = download_ref(params.accession_ref)

    reads_raw_channel = Channel.fromPath("${params.reads_raw}/*.fasta", type:'file')
                                .mix(download_ref_ch)
                                .collect()
    reads_raw_channel.view()
    combine_fasta(reads_raw_channel)
    
    //combine_fasta_channel = Channel.fromPath("${params.out}/${params.combined_fasta}")
    
    align_mafft(combine_fasta.out)

    //alignment_channel = Channel.fromPath("${params.out}/${params.alignment_mafft}")
    
    trimming_trimal(align_mafft.out)
    
}