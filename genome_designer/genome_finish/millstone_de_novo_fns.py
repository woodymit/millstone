import subprocess
import os
import shutil
import pickle
import re

from Bio import SeqIO
from django.conf import settings
import pysam

from genome_finish import __path__ as gf_path_list
from settings import SAMTOOLS_BINARY
from settings import BASH_PATH
from utils.bam_utils import clipping_stats

GENOME_FINISH_PATH = gf_path_list[0]
VELVETH_BINARY = settings.TOOLS_DIR + '/velvet/velveth'
VELVETG_BINARY = settings.TOOLS_DIR + '/velvet/velvetg'


def get_clipped_reads(bam_filename, output_filename, clipping_threshold=None):
    if clipping_threshold is None:
        stats = clipping_stats(bam_filename, sample_size=10000)
        clipping_threshold = int(stats['mean'] + stats['std'])

    print "get_clipped_reads:\n\tCLIPPING THRESHOLD:", clipping_threshold, "\n"

    cmd = ' | '.join([
            '{samtools} view -h {bam_filename}',
            '{extract_clipped_script} -i stdin -t {clipping_threshold}',
            '{samtools} view -Sb -']).format(
                    samtools=SAMTOOLS_BINARY,
                    bam_filename=bam_filename,
                    clipping_threshold=clipping_threshold,
                    extract_clipped_script=os.path.join(
                            GENOME_FINISH_PATH,
                            'extractClippedReads.py'))

    with open(output_filename, 'w') as fh:
        subprocess.check_call(cmd, stdout=fh, shell=True, executable=BASH_PATH)

    # sort the split reads, overwrite the old file
    subprocess.check_call(
            [SAMTOOLS_BINARY, 'sort', output_filename] +
            [os.path.splitext(output_filename)[0]])


def get_match_counts(bam_filename):
    cmd = ' | '.join([
            '{samtools} view -h {bam_filename}',
            '{assess_alignment_script} -i stdin']).format(
                    samtools=SAMTOOLS_BINARY,
                    bam_filename=bam_filename,
                    assess_alignment_script=os.path.join(
                        GENOME_FINISH_PATH,
                        'assess_alignment.py'))

    output = subprocess.check_output(cmd, shell=True, executable=BASH_PATH)
    return pickle.loads(output)


def get_insertion_location(bam_path):

    BAM_CSOFT_CLIP = 4
    BAM_CHARD_CLIP = 5
    clip_codes = [BAM_CSOFT_CLIP, BAM_CHARD_CLIP]

    samfile = pysam.AlignmentFile(bam_path)
    contig_length = 0

    left_read = None
    right_read = None
    for read in samfile:
        contig_length = max(contig_length, read.query_length)
        if read.cigartuples[0][0] not in clip_codes:
            left_read = read
        elif read.cigartuples[-1][0] not in clip_codes:
            right_read = read

    if left_read is None:
        error_string = 'Expected homology on left end of contig not found'
    elif right_read is None:
        error_string = 'Expected homology on right end of contig not found'
    else:
        left_start_pos = left_read.reference_start
        right_end_pos = right_read.reference_end
        chromosome_seqrecord_id = samfile.getrname(left_read.reference_id)

        if left_start_pos > right_end_pos:
            error_string = ('Contig alignment to reference suggests a ' +
                    'structural variant more complex than a simple insertion')
        elif contig_length < right_end_pos - left_start_pos:
            error_string = ('Contig not long enough to span putative ' +
                    'insertion region, suggesting insertion is longer than ' +
                    'assembled contig')
        else:
            locations_dict = {
                    'chromosome_seqrecord_id': chromosome_seqrecord_id,
                    'left_end': left_start_pos,
                    'right_end': right_end_pos,
            }
            return locations_dict

    return {'error_string': error_string}


def make_sliced_fasta(fasta_path, seqrecord_id, left_index, right_index,
        output_fasta_path=None):
    with open(fasta_path, 'r') as fh:
        seqrecord = SeqIO.parse(fh, 'fasta').next()
        while(seqrecord.id != seqrecord_id):
            seqrecord = SeqIO.parse(fh, 'fasta').next()

    seqrecord = seqrecord[left_index:right_index]

    if output_fasta_path:
        fasta_prefix, fasta_ext = os.path.splitext(output_fasta_path)
        sliced_fasta_path = output_fasta_path
    else:
        fasta_prefix, fasta_ext = os.path.splitext(fasta_path)
        sliced_fasta_path = (fasta_prefix + '.split_' + str(left_index) +
                '_to_' + str(right_index) + fasta_ext)
        if os.path.exists(sliced_fasta_path):
            raise Exception("Sliced fasta already exists with path:",
                    sliced_fasta_path)

    with open(sliced_fasta_path, 'w') as fh:
        SeqIO.write(seqrecord, fh, 'fasta')

    return sliced_fasta_path


def get_local_contig_placement(ref_fasta, contig_fasta):
    cmd = ' '.join([
            os.path.join(settings.TOOLS_DIR, 'seqan/place_contig'),
            '-ref', ref_fasta,
            '-read', contig_fasta
        ])

    output_raw = subprocess.check_output(cmd, shell=True, executable=BASH_PATH)
    output_raw = output_raw.strip('\n')
    output_list = output_raw.split('\n')

    assert len(output_list) == 3, ('output generated by place_contig script ' +
            'not in expected format of 3 newline delimited integers.  ' +
            'output was split as:' + str(output_list))

    return {
        'ref_split_pos': int(output_list[0]),
        'contig_start_pos': int(output_list[1]),
        'contig_end_pos': int(output_list[2])
    }


def get_unmapped_reads(bam_filename, output_filename):
    cmd = '{samtools} view -h -b -f 0x4 {bam_filename}'.format(
            samtools=SAMTOOLS_BINARY,
            bam_filename=bam_filename)
    with open(output_filename, 'w') as output_fh:
        subprocess.call(
                cmd, stdout=output_fh, shell=True, executable=BASH_PATH)


def get_split_reads(bam_filename, output_filename):
    """Isolate split reads from a sample alignment.

    This uses a python script supplied with Lumppy, that is run as a
    separate process.

    NOTE THAT THIS SCRIPT ONLY WORKS WITH BWA MEM.
    """

    # Use lumpy bwa-mem split read script to pull out split reads.
    filter_split_reads = ' | '.join([
            '{samtools} view -h {bam_filename}',
            'python {lumpy_bwa_mem_sr_script} -i stdin',
            '{samtools} view -Sb -']).format(
                    samtools=SAMTOOLS_BINARY,
                    bam_filename=bam_filename,
                    lumpy_bwa_mem_sr_script=
                            settings.LUMPY_EXTRACT_SPLIT_READS_BWA_MEM)

    try:
        with open(output_filename, 'w') as fh:
            subprocess.check_call(
                    filter_split_reads, stdout=fh, shell=True,
                    executable=BASH_PATH)

        # sort the split reads, overwrite the old file
        subprocess.check_call(
                [SAMTOOLS_BINARY, 'sort', output_filename,
                 os.path.splitext(output_filename)[0]])

    except subprocess.CalledProcessError:
        raise Exception('Exception caught in split reads generator, ' +
                        'perhaps due to no split reads')


def _parse_sam_line(line):
    parts = line.split()
    return {
        'read_id': parts[0],
        'flags': parts[1]
    }


def add_paired_mates(input_sam_path, source_bam_filename, output_sam_path):
    """Creates a file at output_sam_path that contains all the reads in
    input_sam_path, as well as all of their paired mates.

    The resulting sam is not sorted and may contain duplicates. Clients
    should filter elsewhere.

    TODO: This is currently nasty inefficient. Is there a better way to do
    this, e.g. leverage indexing somehow?

    TODO: This is potentially memory-overwhelming. Need to think about
    the de novo assembly feature more carefully before making it user-facing.
    """
    # Strategy:
    # 1. Copy input to output, to preserve existing reads.
    # 2. Load all ids into a dictionary, storing information about whether
    #    both pairs are already in the output.
    #    NOTE: This step could potentially overwhelm memory for large datasets.
    # 3. Loop through the bam file, checking against the id dictionary to
    #    see whether or not we want to keep that data.

    # 1. Copy input to output, to preserve existing reads.
    shutil.copyfile(input_sam_path, output_sam_path)

    # 2. Create dictionary of read ids.
    read_id_to_flags_map = {}
    with open(input_sam_path) as fh:
        for line in fh:
            if re.match('@', line):
                continue
            parsed_line = _parse_sam_line(line)
            read_id = parsed_line['read_id']
            flags_value = parsed_line['flags']
            if read_id not in read_id_to_flags_map:
                read_id_to_flags_map[read_id] = []
            if flags_value not in read_id_to_flags_map[read_id]:
                read_id_to_flags_map[read_id].append(flags_value)
            assert 1 <= len(read_id_to_flags_map[read_id]) <= 2

    # 3. Loop through bam file, appending lines whose pairs aren't already
    #    present.
    with open(output_sam_path, 'a') as output_fh:
        # Use samtools to read the bam.
        read_bam_cmd = [
            SAMTOOLS_BINARY,
            'view',
            source_bam_filename
        ]
        samtools_proc = subprocess.Popen(read_bam_cmd, stdout=subprocess.PIPE)

        for sam_line in samtools_proc.stdout:
            # Use flags as unique identifier for each read (i.e. which of the
            # pairs for a read id that we are looking at.
            parsed_line = _parse_sam_line(sam_line)
            read_id = parsed_line['read_id']
            flags = parsed_line['flags']
            if read_id not in read_id_to_flags_map:
                continue
            observed_flags = read_id_to_flags_map[read_id]
            if flags not in observed_flags:
                output_fh.write(sam_line)
                # Update flags in case we encounter same one again.
                observed_flags.append(flags)


def run_velvet(reads, output_dir, opt_dict=None):
    default_opts = {
        'velveth': {
            'hash_length': 21
        },
    }

    if not opt_dict:
        opt_dict = default_opts

    if 'hash_length' not in opt_dict['velveth']:
        raise Exception('If passing an option dictionary, an output_dir_name key \
        must exist in the velveth sub-dictionary')

    cmd = ' '.join([
            VELVETH_BINARY,
            output_dir,
            str(opt_dict['velveth']['hash_length'])] +
            ['-' + key + ' ' + str(opt_dict['velveth'][key])
             for key in opt_dict['velveth']
             if key not in ['hash_length']] +
            ['-bam', reads])
    subprocess.check_call(cmd, shell=True, executable=BASH_PATH)

    cmd = ' '.join([
            VELVETG_BINARY,
            output_dir] +
            ['-' + key + ' ' + str(opt_dict['velvetg'][key])
             for key in opt_dict['velvetg']])
    subprocess.check_call(cmd, shell=True, executable=BASH_PATH)
