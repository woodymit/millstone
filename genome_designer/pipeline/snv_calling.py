"""
Functions for calling SNPs.
"""

import os
import re
import subprocess

from celery import task

from main.models import clean_filesystem_location
from main.models import Dataset
from main.models import ensure_exists_0775_dir
from main.models import get_dataset_with_type
from read_alignment import get_insert_size
from main.s3 import project_files_needed
from variant_effects import run_snpeff
from scripts.vcf_parser import parse_alignment_group_vcf
from scripts.jbrowse_util import add_vcf_track
from settings import TOOLS_DIR

# TODO: These VCF types should be set somewhere else. snpeff_util and vcf_parser
# also use them, but where should they go? settings.py seems logical, but it
# cannot import from models.py... -dbg

# Dataset type to use for snp calling.
VCF_DATASET_TYPE = Dataset.TYPE.VCF_FREEBAYES
# Dataset type to use for snp annotation.
VCF_ANNOTATED_DATASET_TYPE = Dataset.TYPE.VCF_FREEBAYES_SNPEFF
# Dataset type for results of finding SVs.
VCF_PINDEL_TYPE = Dataset.TYPE.VCF_PINDEL
VCF_DELLY_TYPE = Dataset.TYPE.VCF_DELLY

# Returns a dictionary of common parameters required for all variant callers
# (freebayes, pindel, delly)
def get_common_tool_params(alignment_group):
    alignment_type = Dataset.TYPE.BWA_ALIGN
    return {
            'alignment_group': alignment_group,
            'alignment_type': alignment_type,
            'fasta_ref': _get_fasta_ref(alignment_group),
            'output_dir': _create_output_dir(alignment_group),
            'bam_files': _find_valid_bam_files(alignment_group, alignment_type),
            }

# Returns a tuple of variant tools params to pass into find_variants_with_tool
def get_variant_tool_params():
    return (
            ('freebayes', Dataset.TYPE.VCF_FREEBAYES, run_freebayes),
            ('pindel', Dataset.TYPE.VCF_PINDEL, run_pindel),
            ('delly', Dataset.TYPE.VCF_DELLY, run_delly),
            )

def _get_fasta_ref(alignment_group):
    # Grab the reference genome fasta for the alignment.
    return get_dataset_with_type(
            alignment_group.reference_genome,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

def _create_output_dir(alignment_group):
    # Prepare a directory to put the output files.
    # We'll put them in /projects/<project_uid>/alignment_groups/vcf/<variant tool>/
    #     <alignment_type>.vcf
    # We'll save these for now, maybe it's not necessary later.
    vcf_dir = os.path.join(alignment_group.get_model_data_dir(), 'vcf')
    ensure_exists_0775_dir(vcf_dir)
    return vcf_dir

def _find_valid_bam_files(alignment_group, alignment_type):
    """ Returns a list of tuples (bam file path, corresponding sample) """
    sample_alignment_list = (
            alignment_group.experimentsampletoalignment_set.all())

    # Filter out mis-aligned files.
    # TODO: Should we show in the UI that some alignments failed and are
    # being skipped?
    def _is_successful_alignment(sample_alignment):
        bam_dataset = get_dataset_with_type(sample_alignment, alignment_type)
        return bam_dataset.status == Dataset.STATUS.READY
    sample_alignment_list = [sample_alignment for sample_alignment in
            sample_alignment_list if _is_successful_alignment(sample_alignment)]

    if len(sample_alignment_list) == 0:
        raise Exception('No successful alignments, Freebayes cannot proceed.')

    # Get handles for each of the bam files.
    def _get_bam_location(sample_alignment):
        bam_dataset = get_dataset_with_type(sample_alignment, alignment_type)
        return bam_dataset.get_absolute_location()
    bam_files = [(_get_bam_location(sample_alignment), sample_alignment.experiment_sample)
            for sample_alignment in sample_alignment_list if sample_alignment]

    # Keep only valid bam_files
    valid_bam_files = []
    for bam_file, sample in bam_files:
        if bam_file is None:
            continue
        if not os.stat(bam_file).st_size > 0:
            continue
        valid_bam_files.append((bam_file, sample))
    assert len(valid_bam_files) == len(sample_alignment_list), (
            "Expected %d bam files, but found %d" % (
                    len(sample_alignment_list), len(bam_files)))
    return valid_bam_files

@task
@project_files_needed
def find_variants_with_tool(alignment_group, variant_params):
    common_params = get_common_tool_params(alignment_group)
    tool_name, vcf_dataset_type, tool_function = variant_params

    # Create subdirectory for this tool
    tool_dir = os.path.join(common_params['output_dir'], tool_name)
    ensure_exists_0775_dir(tool_dir)
    vcf_output_filename = os.path.join(tool_dir, common_params['alignment_type'] + '.vcf')

    # Run the tool
    tool_function(common_params['fasta_ref'],
            common_params['bam_files'], tool_dir, vcf_output_filename)

    # Add dataset
    # If a Dataset already exists, delete it, might have been a bad run.
    existing_set = Dataset.objects.filter(
            type=vcf_dataset_type,
            label=vcf_dataset_type,
            filesystem_location=clean_filesystem_location(vcf_output_filename),
    )
    if len(existing_set) > 0:
        existing_set[0].delete()

    dataset = Dataset.objects.create(
            type=vcf_dataset_type,
            label=vcf_dataset_type,
            filesystem_location=clean_filesystem_location(vcf_output_filename),
    )
    alignment_group.dataset_set.add(dataset)

    # Do the following only for freebayes; right now just special if condition
    if tool_name == 'freebayes':
        # For now, automatically run snpeff if a genbank annotation is available.
        # If no annotation, then skip it, and pass the unannotated vcf type.
        if alignment_group.reference_genome.is_annotated():
            run_snpeff(alignment_group, Dataset.TYPE.BWA_ALIGN)
            vcf_dataset_type = VCF_ANNOTATED_DATASET_TYPE
        else:
            vcf_dataset_type = VCF_DATASET_TYPE

    # Tabix index and add the VCF track to Jbrowse
    add_vcf_track(alignment_group.reference_genome, alignment_group,
        vcf_dataset_type)

    # Parse the resulting vcf
    parse_alignment_group_vcf(alignment_group, vcf_dataset_type)


def run_freebayes(fasta_ref, bam_files, vcf_output_dir, vcf_output_filename):
    """Run freebayes using the bam alignment files keyed by the alignment_type
    for all Genomes of the passed in ReferenceGenome.

    NOTE: If a Genome doesn't have a bam alignment file with this
    alignment_type, then it won't be used.
    """
    vcf_dataset_type = VCF_DATASET_TYPE

    # Build up the bam part of the freebayes binary call.
    bam_part = []
    for bam_file, sample in bam_files:
        bam_part.append('--bam')
        bam_part.append(bam_file)

    # Build the full command and execute it for all bam files at once.
    full_command = (['%s/freebayes/freebayes' %  TOOLS_DIR] + bam_part + [
        '--fasta-reference', fasta_ref,
        '--pvar', '0.001',
        '--ploidy', '2',
        '--min-alternate-fraction', '.3',
        '--hwe-priors-off',
        '--binomial-obs-priors-off',
        '--use-mapping-quality',
        '--min-base-quality', '25',
        '--min-mapping-quality', '30'
    ])

    with open(vcf_output_filename, 'w') as fh:
        subprocess.check_call(full_command, stdout=fh)


def _filter_small_variants(vcf_file, cutoff):
    """Go through each line of vcf, and remove small structural variants"""
    vcf_file_tmp = vcf_file + '.tmp'
    with open(vcf_file_tmp, 'w') as fh:
        for line in open(vcf_file):
            match = re.search('SVLEN=(-?[0-9]+);', line)
            # Check if SVLEN > cutoff
            if not match or abs(int(match.group(1))) > cutoff:
                fh.write(line)

    # move temporary file back to vcf_file path
    subprocess.check_call(['mv', vcf_file_tmp, vcf_file])


def run_pindel(fasta_ref, bam_files, vcf_output_dir, vcf_output_filename):
    """Run pindel to find SVs."""
    vcf_dataset_type = VCF_PINDEL_TYPE

    if not os.path.isdir('%s/pindel' % TOOLS_DIR):
        raise Exception('Pindel is not installed. Aborting.')

    # Create pindel config file
    pindel_config = os.path.join(vcf_output_dir, 'pindel_config.txt')
    with open(pindel_config, 'w') as fh:
        for bam_file, sample in bam_files:
            insert_size = get_insert_size(bam_file)
            fh.write('%s %s %s\n' % (bam_file, insert_size, sample.uid))

    # Build the full pindel command.
    pindel_root = vcf_output_filename[:-4]  # get rid of .vcf extension
    subprocess.check_call(['%s/pindel/pindel' % TOOLS_DIR,
        '-f', fasta_ref,
        '-i', pindel_config,
        '-c', 'ALL',
        '-o', pindel_root
    ])

    # convert all different structural variant types to vcf
    subprocess.check_call(['%s/pindel/pindel2vcf' % TOOLS_DIR,
        '-P', pindel_root,
        '-r', fasta_ref,
        '-R', 'name',
        '-d', 'date'
    ])
    _filter_small_variants(vcf_output_filename, 10)


def run_delly(fasta_ref, bam_files, vcf_output_dir, vcf_output_filename):
    """Run delly to find SVs."""
    vcf_dataset_type = VCF_DELLY_TYPE

    if not os.path.isdir('%s/delly' % TOOLS_DIR):
        raise Exception('Delly is not installed. Aborting.')

    delly_root = vcf_output_filename[:-4]  # get rid of .vcf extension
    transformations = ['DEL', 'DUP', 'INV']
    vcf_outputs = map(lambda transformation:
            '%s_%s.vcf' % (delly_root, transformation), transformations)

    # Rename bam files, because Delly uses the name of the file as the sample uid.
    # Use cp instead of mv, because other sv callers will be reading from the
    #   original bam file.
    new_bam_files = []
    for bam_file, sample in bam_files:
        new_bam_file = os.path.join(os.path.dirname(bam_file), sample.uid + '.bam')
        subprocess.check_call(['cp', bam_file, new_bam_file])
        subprocess.check_call(['cp', bam_file + '.bai', new_bam_file + '.bai'])
        new_bam_files.append(new_bam_file)

    # run delly for each type of transformation
    for transformation, vcf_output in zip(transformations, vcf_outputs):
        # not checked_call, because delly errors if it doesn't find any SVs
        subprocess.call(['%s/delly/delly' % TOOLS_DIR,
            '-t', transformation,
            '-o', vcf_output,
            '-g', fasta_ref] + new_bam_files)

    # combine the separate vcfs for each transformation
    vcf_outputs = filter(lambda file: os.path.exists(file), vcf_outputs)
    if vcf_outputs:
        with open(vcf_output_filename, 'w') as fh:
            subprocess.check_call(['vcf-concat'] + vcf_outputs, stdout=fh)
    else:
        # hack: create empty vcf
        subprocess.check_call(['touch', delly_root])
        subprocess.check_call(['%s/pindel/pindel2vcf' % TOOLS_DIR,
            '-p', delly_root,  # TODO does this work?
            '-r', fasta_ref,
            '-R', 'name',
            '-d', 'date'
        ])

    # Delete temporary renamed bam files
    for bam_file in new_bam_files:
        subprocess.check_call(['rm', bam_file])
        subprocess.check_call(['rm', bam_file + '.bai'])
