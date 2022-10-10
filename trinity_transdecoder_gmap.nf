#!/usr/bin/env nextflow

//params.reads = "$baseDir/data/external/trimmed_reads/*_R{1,2}*.fq"
params.reads = "$baseDir/data/external/other_calanus/*_{1,2}_val_{1,2}.fq"
params.assembly = "$baseDir/data/external/assembly.fasta"

//Channel
//    .fromPath("$baseDir/data/external/*.fasta")
//    .set { assemblies}
//
/*
* Create the `read_pairs_ch` channel that emits tuples containing three elements:
* the pair ID, the first read-pair file and the second read-pair file
*/

Channel
    .fromFilePairs( params.reads, flat: true )
    .set { samples_ch }


// Run with 2.13.2? 

process Trinity {
    cpus 16
    executor 'slurm'
    scratch true
    time '240h'
    clusterOptions = '-A snic2021-5-453 -C mem256GB'
    module 'bioinfo-tools:trinity'
    tag "$sampleId"

    input:
    set sampleId, file(left), file(right) from samples_ch
    
    
    publishDir "results/trinity/$sampleId", mode: 'copy'
    
    output:
    path '*.fasta' into trinity

    """
    Trinity \
    --CPU ${task.cpus} \
    --max_memory 220G \
    --seqType fq \
    --left $left \
    --right $right \
    --full_cleanup \
    --output ${sampleId}.trinity
    """
}

//trinity_transcripts = Channel.fromPath("$baseDir/results/trinity/**/*.fasta")

// 2.13.2 gives error

process longest_transcripts {
    cpus 1
    executor 'slurm'
    scratch true
    time '15min'
    clusterOptions = '-A snic2021-5-453'
    tag "$transcriptome"
    module 'bioinfo-tools:trinity/2.11.0'

    input:
    path transcriptome from trinity
    
    publishDir "data/interim/longest_isoforms/", mode: 'copy'
    
    output:
    tuple path(transcriptome), path("*.longest_isoforms.fasta") into (longest_isoforms, remember_name)


    """
    module purge
    module load bioinfo-tools trinity/2.11.0
    get_longest_isoform_seq_per_trinity_gene.pl ${transcriptome} > ${transcriptome}.longest_isoforms.fasta
    """
}

process get_orfs {
    cpus 1
    executor 'slurm'
    scratch true
    time '24h'
    clusterOptions = '-A snic2021-5-453'
    tag "$transcriptome"
    module 'bioinfo-tools:TransDecoder'
    tag "$longest_isoform"

    input:
    tuple path(name), path(longest_isoform) from longest_isoforms

    publishDir "$baseDir/data/interim/transdecoder/", mode: 'copy'

    output:
    tuple path(longest_isoform), path("**/longest_orfs.pep"), path("**/longest_orfs.cds"), path("**/longest_orfs.gff3"), path("**/base_freqs.dat"), path(longest_isoform) into (lo_pfam, lo_sprot, lo_predict)

    """
    TransDecoder.LongOrfs \
    -t ${longest_isoform}
    """
}
Channel.fromPath("/sw/data/blast_databases/uniprot_sprot.*").into{blastdb;pass_files}

process diamond_prep {
    cpus 1
    executor 'slurm'
    scratch true
    time '1h'
    clusterOptions = '-A snic2021-5-453'
    module 'bioinfo-tools:diamond'
   //stageInMode 'copy'

    input:
    path x from blastdb.collect()

    output:
    path "uniprot_sprot.*" into prep_db

    """
    diamond prepdb -d ./uniprot_sprot
    """

}

process search_sprot {
    cpus 16
    executor 'slurm'
    scratch true
    time '24h'
    clusterOptions = '-A snic2021-5-453'
    tag "$longest_orf"
    module 'bioinfo-tools:diamond'
    stageInMode 'copy'
    tag "$name"

    input:
    tuple val(name), path(pep), path(cds), path(gff), path(freqs), path(lo) from lo_sprot
    path x from prep_db.collect()
    path y from pass_files.collect()

    output:
    file "blastp_out" into sprot
    
    //publishDir "$baseDir/data/interim/sprot/", mode: 'copy'

    """
    diamond blastp \
        --query ${pep} \
        --db ./uniprot_sprot \
        --max-target-seqs 1 \
        --outfmt 6 \
        --evalue 1e-5 \
        --threads ${task.cpus} \
        --out blastp_out
    """
}

process search_pfam {
    cpus 8
    executor 'slurm'
    scratch true
    time '24h'
    clusterOptions = '-A snic2021-5-453'
    tag "$longest_orf"
    module 'bioinfo-tools:hmmer'
    //stageInMode 'copy' 
    tag "$name"

    input:
    tuple val(name), path(pep), path(cds), path(gff), path(freqs), path(lo)  from lo_pfam

    output:
    file "pfam_out" into pfam

    //publishDir "$baseDir/data/interim/pfam/", mode: 'copy'

    """
     hmmsearch \
        --cpu ${task.cpus} \
        --domtblout pfam_out \
        /sw/data/Pfam/31.0/Pfam-A.hmm \
        ${pep}
    """

}

process transdecoder_predict {
    cpus 1
    executor 'slurm'
    scratch true
    time '1h'
    clusterOptions = '-A snic2021-5-453'
    tag "$longest_orf"
    module 'bioinfo-tools:TransDecoder'
    //stageInMode 'copy' 
    tag "$name"

    input:
    tuple val(name), path(pep), path(cds), path(gff), path(freqs), path(lo)  from lo_predict
    file pfam
    file sprot


    output:
    path "*.transdecoder.*" into predicted_transcripts
    file "*.transdecoder.cds" into predicted_cds

    publishDir "$baseDir/data/interim/transdecoder/", mode: 'copy'

    """
    TransDecoder.Predict -t ${pep} \
        -t ${name} \
        --single_best_only \
        --output_dir \$PWD \
        --retain_pfam_hits ${pfam} \
        --retain_blastp_hits ${sprot}
    """

}

process gmap_index {
    conda '/domus/h1/fele1522/.local/bin/mambaforge/envs/gmap2'
    cpus 8
    executor 'slurm'
    scratch true
    time '240h'
    clusterOptions = '-A snic2021-5-453'

    input:
    path assembly from param.assembly

    output:
    path '*' into gmap_index

    """
    module purge
    gmap_build -s none -d ${assembly}.index --dir \$PWD ${assembly}
    """
}

process gmap_map { 
    conda '/domus/h1/fele1522/.local/bin/mambaforge/envs/gmap'
    cpus 16
    executor 'slurm'
    scratch true
    time '48h'
    clusterOptions = '-A snic2021-5-453'

    input:
    //path genome from gmap_index
    tuple path(cdna_file), path(genome) from predicted_cds.combine(gmap_index)

    publishDir "$baseDir/results/gmap/", mode: 'copy'
    output:
    path '*.gff3'

    """
    gmapl --dir \$PWD -f 2 -d ${genome} ${cdna_file} -t ${task.cpus} > ${genome}.${cdna_file}.gff3
    """
}
