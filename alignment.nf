// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
params.fastq = "/data/sample.fastq.gz"
params.reference_fa = "data/GRCh38_full_analysis_set_plus_decoy_hla.fa"
params.alignment_threads = 8
params.merge_threads = 8
params.reads_per_split = 1000

process SplitFastq {
    time '1h'
    memory '2 GB'

  input:
    tuple (path(fastq), val(reads_per_split))
  output:
    path "output/*", emit: fastqs
  script:
    """
    ont_alignment split-fastq --fastq_file $fastq --outdir output/ --reads_per_split $reads_per_split
    """
}


process MiniMap{
  cpus params.alignment_threads
  time '1h'
  memory '6 GB'

  input:
    tuple path(fastq), path(reference_fa)
  output:
    path "sorted.bam"
  script:
    """
    minimap2 -y -t $task.cpus -ax map-ont $reference_fa $fastq > minimap.sam
    samtools view -bSh minimap.sam > minimap.bam
    samtools sort minimap.bam > sorted.bam
    """
}

process MergeBams{
  cpus params.alignment_threads
  time '1h'
  memory '6 GB'

  input:
    path(bams, stageAs: "?/*")
  output:
    path "merged.bam"
  script:
    """
        ont_alignment merge-bams --output merged.bam --tempdir temp --ncores $task.cpus --bam ${bams.join(' --bam ')}
    """

}


process MarkDuplicates{
  time '1h'
  memory '6 GB'

  input:
    path(bam)
  output:
    tuple path("duplicates_marked.bam"), path("duplicates_marked.bam.bai")
  script:
    """
        picard -Xmx4G -Xms4G MarkDuplicates \
        INPUT=${bam} \
        OUTPUT=duplicates_marked.bam \
        METRICS_FILE=markduplicates_metrics.txt \
        REMOVE_DUPLICATES=False \
        ASSUME_SORTED=True \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=tempdir \
        MAX_RECORDS_IN_RAM=150000
        samtools index duplicates_marked.bam
    """
}


process MapulaCount{
  time '1h'
  memory '6 GB'

  input:
    tuple path(bam), path(bai), path(reference)
  output:
    path "mapula_count.csv"
  script:
    """
      mapula count ${bam}  -r ${reference} -n mapula_count
    """
}



workflow {

    Channel.fromPath( params.fastq )
            | map { tuple( it, params.reads_per_split ) }
            | set { sample_fastq_file }

    Channel.fromPath( params.reference_fa) | set {reference_ch}
    Channel.fromPath( params.reference_fa) | set {reference_ch_2}


    SplitFastq(sample_fastq_file) | flatten | combine(reference_ch) | MiniMap | collect | MergeBams | MarkDuplicates | combine(reference_ch) | MapulaCount

}