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
    path "tempdir/nanoplot_NanoPlot-report.html", emit: html
    path "tempdir/nanoplot_NanoStats.txt", emit: stats
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
      ont_alignment merge-unmapped-stats --stats_files ${stats.join(' --stats_files ')} --output merged_stats.csv
    """
}


process AlignmentReport{
  time '48h'
  memory '16 GB'
  input:
    tuple file(mapula_json), file(unmapped_stats), file(reference), val(sample_name)
  output:
    path('alignment_report.html')
  script:
    """
    ont_alignment plot-alignment-report --mapula_json ${mapula_json} --report_html alignment_report.html \
      --reference ${reference}  --sample_name ${sample_name} --unmapped_stats ${unmapped_stats}
    """
}


process BamFastqc{
  time '48h'
  memory '16 GB'
  input:
    tuple file(bam), file(bai)
  output:
    path('output/*fastqc.html'), emit: fastqc_html
    path('output/*fastqc.zip'), emit: fastqc_zip
  script:
    """
      mkdir output
      fastqc --threads 1 -f bam -o output ${bam}
    """
}

process SamtoolsFlagstat{
  time '48h'
  memory '16 GB'
  input:
    tuple file(bam), file(bai)
  output:
    path('flagstat.txt')
  script:
    """
      samtools flagstat ${bam} > flagstat.txt
    """
}



process MultiQc{
  time '48h'
  memory '16 GB'
  input:
    file(fastqc_zip)
    file(nanoplot_html)
    file(flagstat)
  output:
    path('output/multiqc.html')
  script:
    """
      mkdir output
      multiqc -o output -n multiqc .
    """

}

process IgvtoolsCount{
  time '48h'
  memory '16 GB'
  input:
    tuple file(bam), file(reference)
  output:
    file("merged.bam.tdf")
  script:
    """
      igvtools count ${bam} merged.bam.tdf ${reference}
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

        aligned_bams | UnmappedStats | collect | MergeUnmappedStats | set {unmapped_stats_merged}

        aligned_bams | combine(reference_ch) | MapulaCount | set{mapula_jsons}

        aligned_bams | collect | MergeBams | set {bamfile}

        mapula_jsons | collect | MapulaMerge | set {mapula_json_merged}

        mapula_json_merged | combine(unmapped_stats_merged) | combine(reference_ch) | combine(["SA123"]) | AlignmentReport

        bamfile | combine([0]) | combine([0]) | NanoPlot

        bamfile | combine(reference_ch) | IgvtoolsCount

        bamfile | BamFastqc

        bamfile | SamtoolsFlagstat

        MultiQc(BamFastqc.out.fastqc_zip, NanoPlot.out.stats, SamtoolsFlagstat.out)

    emit:
        MergeBams.out

}


workflow{
        Channel.fromPath( params.reference_fa) | set {reference_ch}
        Channel.fromPath( params.samplesheet) | set {samplesheet_ch}

        align(reference_ch, samplesheet_ch)

        postprocess(align.out)
        postprocess_filter(align.out)

}