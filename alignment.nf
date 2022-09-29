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
  memory '16 GB'

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

process UnmappedStats{
  time '6h'
  memory '16 GB'

  input:
    path(sorted_bam)
  output:
    path "unmapped_stats.csv"
  script:
    """
      bamtools split -in ${sorted_bam} -mapped
      bedtools bamtofastq -i *UNMAPPED.bam -fq temp.unmapped.fq
      fastcat -s temp.unmapped.fq -r unmapped_stats.csv -x temp.unmapped.fq >> temp.uncompressed.fastq
    """
}


process MergeBams{
  cpus params.merge_threads
  time '96h'
  memory '16 GB'

  input:
    path(bams, stageAs: "?/*")
  output:
    tuple path("merged.bam"), path("merged.bam.bai")
  script:
    """
        ont_alignment merge-bams --output merged.bam --tempdir temp --ncores $task.cpus --bam ${bams.join(' --bam ')}
        samtools index merged.bam
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


process MapulaMerge{
  time '48h'
  memory '16 GB'

  input:
    path(mapula_json, stageAs: "?/*")
  output:
    path "merged_mapula.json"
  script:
    """
      mapula merge $mapula_json -f all -n merged_mapula
    """
}


process MapulaCount{
  time '48h'
  memory '16 GB'

  input:
    tuple path(bam), path(reference)
  output:
    path "mapula_count.json"
  script:
    """
      mapula count ${bam}  -r ${reference} -n mapula_count -f json
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


process MergeUnmappedStats{
  time '48h'
  memory '16 GB'
  input:
    path(stats, stageAs: "?/*")
  output:
    path "merged_stats.csv"
  script:
    """
      ont_alignment merge_unmapped_stats --stats ${stats.join(' --stats ')} --output merged_stats.csv
    """
}



workflow align {
    take:
        reference_ch
        sample_sheet_ch

    main:

        sample_sheet_ch
        | splitCsv( header: true, sep: ',' )
        | map { row -> tuple( row.sample_id, row.lane_id, file(row.fastq)) }
        | combine(reference_ch)
        | set{sample_fastq_ch}

        MiniMap(sample_fastq_ch) | set {aligned_bams}

        aligned_bams | UnmappedStats | collect | MergeUnmappedStats

        aligned_bams | combine(reference_ch) | MapulaCount | set{mapula_jsons}

        aligned_bams | collect | MergeBams

        mapula_jsons | collect | MapulaMerge

    emit:
        MergeBams.out

}

workflow postprocess_filter {

    take:
        bam_ch
    main:
        bam_ch | combine([params.min_qual]) | combine([params.min_length]) | NanoPlot

}

workflow postprocess {

    take:
        bam_ch
    main:
        bam_ch | combine([0]) | combine([0]) | NanoPlot

}


workflow{
        Channel.fromPath( params.reference_fa) | set {reference_ch}
        Channel.fromPath( params.samplesheet) | set {samplesheet_ch}

        align(reference_ch, samplesheet_ch)

        postprocess(align.out)
        postprocess_filter(align.out)

}