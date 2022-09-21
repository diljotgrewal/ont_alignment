"""Console script for csverve."""

import click

from ont_alignment import utils


@click.group()
def cli():
    pass


@cli.command()
@click.option('--fastq_file', multiple=True, required=True, help='fastq file')
@click.option('--outdir', required=True, help='output dir')
@click.option('--reads_per_split', default=1e7, help='number of reads per split file')
def split_fastq(
        fastq_file, outdir, num_reads_per_split
):
    utils.split_fastq(fastq_file, outdir, num_reads_per_split)
