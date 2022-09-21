"""Console script for csverve."""

import click

from ont_alignment import utils


@click.group()
def cli():
    pass


@cli.command()
@click.option('--fastq_file', required=True, help='fastq file')
@click.option('--outdir', required=True, help='output dir')
@click.option('--reads_per_split', default=1e7, help='number of reads per split file')
def split_fastq(
        fastq_file, outdir, reads_per_split
):
    utils.split_fastq(fastq_file, outdir, reads_per_split=reads_per_split)


@cli.command()
@click.option('--bam', multiple=True, required=True, help='fastq file')
@click.option('--output', required=True, help='fastq file')
@click.option('--tempdir', required=True, help='output dir')
@click.option('--ncores', default=8, help='number of cores')
def merge_bams(
        bam, output, tempdir, ncores=8
):
    utils.merge_bams(bam, output, tempdir, ncores=ncores)
