import errno
import gzip
import itertools
import os
import subprocess
from subprocess import Popen, PIPE


def _run_cmd(cmd, output=None):
    stdout = PIPE
    if output:
        stdout = open(output, "w")

    p = Popen(cmd, stdout=stdout, stderr=PIPE)

    cmdout, cmderr = p.communicate()
    retc = p.returncode

    cmdout = cmdout.decode('utf-8')
    cmderr = cmderr.decode('utf-8')

    if output:
        stdout.close()

    console_out = '-' * 30 + ' cmdout ' + '-' * 30 + '\n'
    console_out += '\n'.join(cmdout.split('\n')) + '\n'
    console_out += '-' * 70

    console_err = '-' * 30 + ' cmderr ' + '-' * 30 + '\n'
    console_err += '\n'.join(cmderr.split('\n')) + '\n'
    console_err += '-' * 70

    if retc:
        raise Exception("command failed.\n {}\n {}".format(console_out, console_err))

    print(console_out)
    print(console_err)


def _makedirs(directory, isfile=False):
    if isfile:
        directory = os.path.dirname(directory)
        if not directory:
            return

    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def _chunks(bamfiles, numcores):
    output = []
    for i in range(0, len(bamfiles), numcores):
        output.append(bamfiles[i:i + numcores])
    return output


def _get_merge_command(bams, output, ncores=1):
    if len(bams) == 1:
        command = ['cp', bams[0], output]
    else:
        command = ['sambamba', 'merge', '-t', str(ncores), output]
        command.extend(bams)

    return command


def _build_shell_script(command, tag, tempdir):
    outfile = os.path.join(tempdir, "{}.sh".format(tag))
    with open(outfile, 'w') as scriptfile:
        scriptfile.write("#!/bin/bash\n")
        if isinstance(command, list) or isinstance(command, tuple):
            command = ' '.join(map(str, command)) + '\n'
        scriptfile.write(command)
    return outfile


def _run_in_gnu_parallel(commands, tempdir, ncores):
    _makedirs(tempdir)

    scriptfiles = []

    for tag, command in enumerate(commands):
        scriptfiles.append(_build_shell_script(command, tag, tempdir))

    parallel_outfile = os.path.join(tempdir, "commands.txt")
    with open(parallel_outfile, 'w') as outfile:
        for scriptfile in scriptfiles:
            outfile.write("sh {}\n".format(scriptfile))

    subprocess.run(['parallel', '--jobs', str(ncores)], stdin=open(parallel_outfile))


def merge_bams(infiles, outfile, tempdir, ncores):
    assert len(infiles) > 0

    if len(infiles) < ncores * 2:
        _run_cmd(_get_merge_command(infiles, outfile, ncores=ncores))
        return

    chunked_infiles = _chunks(list(infiles), ncores)

    commands = []
    outputs = []
    for i, chunk in enumerate(chunked_infiles):
        chunk_tempdir = os.path.join(tempdir, str(i))
        _makedirs(chunk_tempdir)
        output = os.path.join(chunk_tempdir, 'merged.bam')
        outputs.append(output)
        commands.append(_get_merge_command(chunk, output))

    parallel_temp_dir = os.path.join(tempdir, 'gnu_parallel_temp')
    _run_in_gnu_parallel(commands, parallel_temp_dir, ncores)

    command = _get_merge_command(outputs, outfile, ncores=ncores)
    _run_cmd(command)


def split_fastq(in_filename, outdir, reads_per_split=1e7):
    """ Split a fastq file.
    """

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if in_filename.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open

    out_filenames = lambda x: os.path.join(outdir, f'{x}.fastq.gz')

    with opener(in_filename, 'rt') as in_file:
        file_number = 0
        out_file = None
        out_file_read_count = None
        try:
            for name, seq, comment, qual in itertools.zip_longest(*[in_file] * 4):
                if out_file is None or out_file_read_count == reads_per_split:
                    if out_file is not None:
                        out_file.close()
                    out_file = open(out_filenames(file_number), 'wt')
                    out_file_read_count = 0
                    file_number += 1
                out_file.write(name)
                out_file.write(seq)
                out_file.write(comment)
                out_file.write(qual)
                out_file_read_count += 1
        finally:
            if out_file is not None:
                out_file.close()


def merge_unmapped_stats(stats_files, merged_data):
    with open(merged_data, 'wt') as writer:

        header = None
        for infile in stats_files:
            with open(infile, 'rt') as reader:
                if header is None:
                    header = reader.readline()
                    assert header.startswith('read_id')
                    writer.write(header)
                else:
                    assert header == reader.readline()
                for line in reader:
                    writer.write(line)
