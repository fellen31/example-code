#!/usr/bin/env nextflow

params.in = "/crex/proj/uppstore2018098/nobackup/felix/projects/calanus/results/polish/pilon/bases/flye_andreas.medaka.pilon_no_ncbi.fasta"
params.bam = "$baseDir/data/external/assembly.fasta.minimap2.sam.no_secondary.sorted.bam"
params.bamlist = "$baseDir/bams.list"
params.bai = "$baseDir/data/external/assembly.fasta.minimap2.sam.no_secondary.sorted.bam.bai"
params.out = "out"
params.chunkSize = 500

Channel
    .fromPath(params.in)
    .splitFasta(by :params.chunkSize, file: true)
    .into { fasta_ch;fasta_ch2;fasta_ch3}

/*

process preprocess {
    
    executor 'slurm'
    cpus 1
    scratch true
    time '1h'
    clusterOptions = '-A snic2021-5-453' 
    queue 'core'
    module 'bioinfo-tools:samtools'
    tag "${x}"
    //stageInMode 'copy'
    

    input:
    path x from fasta_ch
    path 'bam' from params.bam
    path 'bam.bai' from params.bai

    output:
    tuple path(x), path("*.bam"), path("*.bai"), path ("*.fai") into (preprecessed_ch)
    //file 'split_bam' into bam_ch
    //path "${x}.fai" into fai_ch
    //file 'split_bam.bai' into bai_ch

    """
    samtools faidx ${x}
    cat ${x}.fai | cut -f1 | xargs samtools view -b bam > split.bam
    samtools index split.bam 
    """

}
*/

/*
process pmd {
    
    executor 'slurm'
    cpus 20
    scratch true
    time '24h'
    clusterOptions = '-A snic2021-5-453' 
    queue 'core'
    tag "${x}"
    stageInMode 'copy'
    //errorStrategy 'ignore'
    

    input:
    tuple path(x), path(bam), path(bai), path(fai) from preprecessed_ch
    //path x from fasta_ch2
    //path 'bam' from bam_ch
    //path 'bam.bai' from bai_ch
    //path 'fasta.fai' from fai_ch
    publishDir 'results'
    output:
    path '*'

    """
    INPUT_DIR=`pwd`
    OUTPUT_DIR=`pwd`
    OUTPUT_PREFIX=${x}
    THREADS=${task.cpus}
    BAM=${bam}
    REF=${x}
    singularity exec --nv --bind /usr/lib/locale/ \
    /crex/proj/uppstore2018098/nobackup/felix/projects/calanus/data/external/assembly/flye/pepper_deepvariant_r0.8.sif \
    run_pepper_margin_deepvariant call_variant \
    -b "\${INPUT_DIR}/\${BAM}" \
    -f "\${INPUT_DIR}/\${REF}" \
    -o "\${OUTPUT_DIR}" \
    -p "\${OUTPUT_PREFIX}" \
    -t \${THREADS} \
    --ont_r9_guppy5_sup 
    """

}*/
//polished_ch

/*

process clair3 {
    
    executor 'slurm'
    cpus 16
    scratch true
    time '24h'
    clusterOptions = '-A snic2021-5-453' 
    queue 'core'
    tag "${x}"
    stageInMode 'copy'
    //errorStrategy 'ignore'
    

    input:
    tuple path(x), path(bam), path(bai), path(fai) from preprecessed_ch
    publishDir 'results/clair3/', mode: 'copy'
    output:
    path '*.gz'

    """
    INPUT_DIR=`pwd`
    OUTPUT_DIR=`pwd`
    THREADS=${task.cpus}
    BAM=${bam}
    REF=${x}
    MODEL_NAME="r941_prom_sup_g5014"
    singularity exec /crex/proj/uppstore2018098/nobackup/felix/projects/calanus/data/external/assembly/flye/clair3_latest.sif \
    /opt/bin/run_clair3.sh \
        --bam_fn=\${INPUT_DIR}/\${BAM} \
        --ref_fn=\${INPUT_DIR}/\${REF} \
        --threads=\${THREADS} \
        --platform="ont" \
        --model_path="/opt/models/\${MODEL_NAME}" \
        --output=\${OUTPUT_DIR} --include_all_ctgs --enable_long_indel 
        for file in *.gz; do mv \$file ${x}.\$file; done 
    """

}

*/

process fb {
    
    executor 'slurm'
    cpus 2
    scratch true
    time '24h'
    clusterOptions = '-A snic2021-5-453' 
    queue 'core'
    tag "${x}"
    //stageInMode 'copy'
    module 'bioinfo-tools:samtools:freebayes'
    //errorStrategy 'ignore'
    

    input:
    path x from fasta_ch3
    //path 'fasta' from params.in
    path 'bam' from params.bamlist
    publishDir 'results/freebayes/', mode: 'copy'
    output:
    path '*.vcf'

    """
    samtools faidx ${x}
    cat ${x}.fai | awk -v OFS='\\t' '{print \$1,0,\$2}' > region.bed
    freebayes -f ${x} -L bam -T 0.01 -O -E 0 -t region.bed --use-best-n-alleles 5 > ${x}.vcf
    """

}

// .collectFile(name: params.out, storeDir: "$baseDir/results/")


