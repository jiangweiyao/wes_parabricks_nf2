#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_fastqs = "/home/ubuntu/WES/fastq/*_{R1_001,R2_001}.fastq.gz"
fastq_files = Channel.fromFilePairs(params.input_fastqs, type: 'file').first()
params.reference_folder = "/home/ubuntu/WES/GRCh38_reference_genome"
reference_folder = file(params.reference_folder)
params.reference_name = "GRCh38_full_analysis_set_plus_decoy_hla.fa"
params.out = "parabricks_wes_output"
params.baits_int = "/home/ubuntu/WES/panel-00005_withChr.baits.grch38.interval_list"
baits_int = file(params.baits_int)
params.target_int = "/home/ubuntu/WES/panel-00005_withChr.targets.grch38.interval_list"
target_int = file(params.target_int)
params.target_bed = "/home/ubuntu/WES/panel-00005_withChr.targets.grch38.bed"
target_bed = file(params.target_bed)

workflow {
    fq2bam(fastq_files, reference_folder, params.reference_name)
    metrics(fq2bam.out,reference_folder, params.reference_name, baits_int, target_int, target_bed)
}

process fq2bam {

    //errorStrategy 'ignore'
    //publishDir params.out, mode: 'copy', overwrite: true
    cpus 32
    memory '120 GB'
    input:
    tuple val(name), file(fastq) 
    file(reference_folder)
    val(reference_name)

    output:
    tuple val("${name}"), path("${name}_output/*")

    """
    mkdir ${name}_output
    pbrun fq2bam \
      --ref ${reference_folder}/${reference_name} \
      --in-fq ${fastq} "@RG\\tID:s${name}\\tLB:lib${name}\\tPL:SING\\tSM:${name}\\tPU:FC1.${name}" \
      --bwa-options="-M" \
      --out-duplicate-metrics ${name}_output/${name}_dupemetrics.txt \
      --optical-duplicate-pixel-distance 300 \
      --out-bam ${name}_output/mark_dups_${name}.bam \
      --logfile ${name}_output/fq2bam.log
    """
}

process metrics {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true
    cpus 8
    memory '61 GB'
    input:
    tuple val(name), file(fq2bam_input)
    file(reference_folder)
    val(reference_name)
    file(baits_int)
    file(target_int)
    file(target_bed)

    output:
    file("${name}_metric/*")

    """
    mkdir ${name}_metric

    gatk --java-options "-Xmx16G" CollectHsMetrics \
         --INPUT mark_dups_${name}.bam \
         --OUTPUT ${name}_metric/${name}_hs_metrics.txt \
         --REFERENCE_SEQUENCE ${reference_folder}/${reference_name} \
         --BAIT_INTERVALS ${baits_int} \
        --TARGET_INTERVALS ${target_int} \
        --MINIMUM_BASE_QUALITY 10 \
        --MINIMUM_MAPPING_QUALITY 10 \
        --PER_TARGET_COVERAGE ${name}_metric/${name}_hs_metrics_per_target.txt \
        --COVERAGE_CAP 1000 \
        --METRIC_ACCUMULATION_LEVEL ALL_READS \
        --NEAR_DISTANCE 250 \
        --CLIP_OVERLAPPING_READS true 



    samtools sort -@ 4 -n -u mark_dups_${name}.bam | fgbio ClipBam -Xmx60g -i /dev/stdin \
         --output=clipped.bam \
         --ref=${reference_folder}/${reference_name} \
         --clipping-mode=Hard \
         --clip-overlapping-reads=true

    samtools sort -@ 8 clipped.bam > clipped_sorted.bam
    samtools index clipped_sorted.bam

    bedtools nuc \
         -fi ${reference_folder}/${reference_name} \
         -bed ${target_bed} \
         -seq > panel_targetlevel_nucstats.txt

    samtools bedcov \
         -Q 0 \
         panel_targetlevel_nucstats.txt \
         clipped_sorted.bam > tmp.txt

    cat <(paste <(cat panel_targetlevel_nucstats.txt | head -n 1) <(echo "14_coverage")) <(cat tmp.txt) > ${name}_metric/${name}_panel_targetlevel_nucstats_covstats.txt


    """
}
