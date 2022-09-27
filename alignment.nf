// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
params.samplesheet = "/data/sample.csv"
params.reference_fa = "data/GRCh38_full_analysis_set_plus_decoy_hla.fa"
params.threads = 8
params.merge_threads = 8
params.reads_per_split = 1000
params.min_qual = 10
params.min_length = 10000

process SplitFastq {
    time '24h'
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
  cpus params.threads
  time '6h'
  memory '6 GB'

  input:
    tuple val(sample_id), val(lane_id), path(fastq), path(reference_fa)
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
  cpus params.merge_threads
  time '96h'
  memory '16 GB'

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
  time '48h'
  memory '16 GB'

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
  time '48h'
  memory '16 GB'

  input:
    tuple path(bam), path(bai), path(reference)
  output:
    path "mapula_count.csv"
  script:
    """
      mapula count ${bam}  -r ${reference} -n mapula_count
    """
}

process NanoPlot{
  time '48h'
  memory '16 GB'
  cpus params.threads

  input:
    tuple path(bam), path(bai), val(min_qual), val(min_length)
  output:
    path "tempdir/nanoplot_NanoPlot-report.html"
  script:
    """
      NanoPlot --bam ${bam} --maxlength 40000 -o tempdir --prefix nanoplot_ -t $task.cpus --minqual $min_qual --minlength $min_length
    """


}

workflow {
    Channel.fromPath( params.reference_fa) | set {reference_ch}

    Channel
        .fromPath( params.samplesheet )
        .splitCsv( header: true, sep: ',' )
        .map { row -> tuple( row.sample_id, row.lane_id, file(row.fastq)) }
        .combine(reference_ch)
        .set{sample_fastq_ch}

    MiniMap(sample_fastq_ch) | collect | MergeBams | MarkDuplicates | set{bam_ch}

    bam_ch | combine([params.min_qual]) | combine([params.min_length]) | NanoPlot

    bam_ch | combine(reference_ch) | MapulaCount




}