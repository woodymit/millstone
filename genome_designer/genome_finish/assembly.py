import datetime
import pickle
import os
import re

from Bio import SeqIO

from genome_finish.millstone_de_novo_fns import add_paired_mates
from genome_finish.millstone_de_novo_fns import get_clipped_reads
from genome_finish.millstone_de_novo_fns import get_unmapped_reads
from genome_finish.millstone_de_novo_fns import run_velvet
from main.models import Contig
from main.models import Dataset
from main.model_utils import get_dataset_with_type
from pipeline.read_alignment import align_with_bwa_mem
from pipeline.read_alignment import get_insert_size_mean_and_stdev
from utils.bam_utils import concatenate_bams
from utils.bam_utils import make_bam
from utils.bam_utils import make_sam
from utils.bam_utils import rmdup
from utils.bam_utils import sort_bam_by_name
from utils.import_util import add_dataset_to_entity
from utils.import_util import prepare_ref_genome_related_datasets

#DEBUG imports
from main.models import AlignmentGroup
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from utils.bam_utils import index_bam
from utils.bam_utils import sort_bam_by_coordinate
from utils.jbrowse_util import add_bam_file_track

# Default args for velvet assembly
VELVET_COVERAGE_CUTOFF = 10
VELVET_HASH_LENGTH = 21

DEFAULT_VELVET_OPTS = {
    'velveth': {
        'hash_length': VELVET_HASH_LENGTH,
        'shortPaired': ''
    },
    'velvetg': {
        'read_trkg': 'yes',
        'cov_cutoff': VELVET_COVERAGE_CUTOFF
    }
}


def generate_contigs(sample_alignment,
        sv_read_classes={}, input_velvet_opts={}, contig_label_base='',
        overwrite=False):

    contig_label_base = '' # Force for now

    # Grab reference genome fasta path, ensure indexed
    reference_genome = sample_alignment.alignment_group.reference_genome
    ref_fasta_dataset = reference_genome.dataset_set.get_or_create(
            type=Dataset.TYPE.REFERENCE_GENOME_FASTA)[0]
    prepare_ref_genome_related_datasets(reference_genome, ref_fasta_dataset)

    # Make data_dir directory to house genome_finishing files
    assembly_dir = os.path.join(
            sample_alignment.get_model_data_dir(),
            'assembly')

    # Make assembly directory if it does not exist
    if not os.path.exists(assembly_dir):
        os.mkdir(assembly_dir)

    data_dir_counter = 0
    data_dir = os.path.join(assembly_dir, str(data_dir_counter))
    while(os.path.exists(data_dir)):
        data_dir_counter += 1
        data_dir = os.path.join(assembly_dir, str(data_dir_counter))
    os.mkdir(data_dir)

    # # Retrieve bwa mem .bam alignment if exists otherwise generate it
    # if not sample_alignment.dataset_set.filter(
    #         type=Dataset.TYPE.BWA_ALIGN).exists():
    #     add_dataset_to_entity(
    #             sample_alignment,
    #             'sample_alignment_for_assembly',
    #             Dataset.TYPE.BWA_ALIGN)
    #     align_with_bwa_mem(
    #             sample_alignment.alignment_group,
    #             sample_alignment,
    #             project=reference_genome.project)

    # Get a bam of sorted SV indicants with pairs
    sv_indicants_bam = get_sv_indicating_reads(sample_alignment,
            sv_read_classes, overwrite=overwrite)

    velvet_opts = DEFAULT_VELVET_OPTS

    # Find insertion metrics
    ins_length, ins_length_sd = get_insert_size_mean_and_stdev(
            sample_alignment)
    velvet_opts['velvetg']['ins_length'] = ins_length
    velvet_opts['velvetg']['ins_length_sd'] = ins_length_sd

    for shallow_key in ['velveth', 'velvetg']:
        if shallow_key in input_velvet_opts:
            for deep_key in input_velvet_opts[shallow_key]:
                velvet_opts[shallow_key][deep_key] = (
                        input_velvet_opts[shallow_key][deep_key])

    # Perform velvet assembly
    contig_files = assemble_with_velvet(
            data_dir, velvet_opts, sv_indicants_bam,
            sample_alignment,
            contig_label_base)

    return contig_files


def get_sv_indicating_reads(sample_alignment, input_sv_indicant_classes={},
        overwrite=False):

    sv_indicant_keys = [
            Dataset.TYPE.BWA_CLIPPED,
            Dataset.TYPE.BWA_SPLIT,
            Dataset.TYPE.BWA_UNMAPPED,
            Dataset.TYPE.BWA_DISCORDANT
    ]

    sv_indicant_class_to_filename_suffix = {
            Dataset.TYPE.BWA_CLIPPED: 'clipped',
            Dataset.TYPE.BWA_SPLIT: 'split',
            Dataset.TYPE.BWA_UNMAPPED: 'unmapped',
            Dataset.TYPE.BWA_DISCORDANT: 'discordant'
    }

    sv_indicant_class_to_generator = {
            Dataset.TYPE.BWA_CLIPPED: get_clipped_reads,
            Dataset.TYPE.BWA_UNMAPPED: get_unmapped_reads
    }

    default_sv_indicant_classes = {
            Dataset.TYPE.BWA_CLIPPED: True,
            Dataset.TYPE.BWA_SPLIT: True,
            Dataset.TYPE.BWA_UNMAPPED: True,
            Dataset.TYPE.BWA_DISCORDANT: True
    }
    default_sv_indicant_classes.update(input_sv_indicant_classes)

    # Grab alignment bam file-path
    alignment_bam = get_dataset_with_type(
            sample_alignment,
            Dataset.TYPE.BWA_ALIGN).get_absolute_location()

    # Get SV indicating reads
    sv_bams_list = []
    alignment_file_prefix = os.path.join(
            sample_alignment.get_model_data_dir(),
            'bwa_align')

    # Helper function for getting sv read datasets
    def _get_or_create_sv_dataset(key):
        dataset_query = sample_alignment.dataset_set.filter(type=key)

        if dataset_query.exists():
            assert len(dataset_query) == 1
            return dataset_query[0]
        else:
            dataset_path = '.'.join([
                    alignment_file_prefix,
                    sv_indicant_class_to_filename_suffix[key],
                    'bam'
                    ])
            generator = sv_indicant_class_to_generator[key]
            generator(alignment_bam, dataset_path)
            return add_dataset_to_entity(
                    sample_alignment,
                    sv_indicant_class_to_filename_suffix[key],
                    key,
                    filesystem_location=dataset_path)

    # Aggregate SV indicants
    for key in sv_indicant_keys:
        if default_sv_indicant_classes[key]:
            dataset = _get_or_create_sv_dataset(key)
            sv_bams_list.append(dataset.get_absolute_location())


    compilation_prefix = '.'.join([
            alignment_file_prefix,
            '_'.join([sv_indicant_class_to_filename_suffix[k] for k in
                    sorted(sv_indicant_keys)
                    if default_sv_indicant_classes[k]])
            ])

    SV_indicants_with_pairs_bam = compilation_prefix + '.with_pairs.bam'
    if os.path.exists(SV_indicants_with_pairs_bam) and not overwrite:
        print ('WARNING: Requested SV indicants bam file: ' +
                SV_indicants_with_pairs_bam +
                ' already exists and will be returned by this function.  ' +
                'To overwrite this file pass the keyword overwrite=True')
        return SV_indicants_with_pairs_bam
    if overwrite:
        print ('WARNING: overwrite is True, so SV read bam datasets ' +
                'are being overwritten')

    # Aggregate SV indicants
    SV_indicants_bam = compilation_prefix + '.bam'
    concatenate_bams(
            sv_bams_list,
            SV_indicants_bam)

    # Remove duplicates
    print 'removing duplicates'
    SV_indicants_no_dups_bam = compilation_prefix + '.no_dups.bam'
    rmdup(SV_indicants_bam, SV_indicants_no_dups_bam)

    # Convert SV indicants bam to sam
    SV_indicants_sam = compilation_prefix + '.no_dups.sam'
    make_sam(SV_indicants_no_dups_bam, SV_indicants_sam)

    # Add mate pairs to SV indicants sam
    print 'adding mate pairs'
    SV_indicants_with_pairs_sam = compilation_prefix + '.with_pairs.sam'
    add_paired_mates(
            SV_indicants_sam, alignment_bam, SV_indicants_with_pairs_sam)

    # Make bam of SV indicants w/mate pairs
    make_bam(SV_indicants_with_pairs_sam, SV_indicants_with_pairs_bam)

    # Sort for velvet assembly
    print 'sorting by name'
    sort_bam_by_name(SV_indicants_with_pairs_bam)

    return SV_indicants_with_pairs_bam


def assemble_with_velvet(data_dir, velvet_opts, sv_indicants_bam,
        sample_alignment, contig_label_base=''):

    timestamp = str(datetime.datetime.now())
    contig_number_pattern = re.compile('^NODE_(\d+)_')

    reference_genome = sample_alignment.alignment_group.reference_genome

    contig_files = []

    # Write sv_indicants filename and velvet options to file
    assembly_metadata_fn = os.path.join(data_dir, 'metadata.txt')
    with open(assembly_metadata_fn, 'w') as fh:
        assembly_metadata = {
            'sv_indicants_bam': sv_indicants_bam,
            'velvet_opts': velvet_opts
        }
        pickle.dump(assembly_metadata, fh)

    # Run velvet assembly on SV indicants
    run_velvet(
            sv_indicants_bam,
            data_dir,
            velvet_opts)

    # Collect resulting contigs fasta
    contigs_fasta = os.path.join(data_dir, 'contigs.fa')
    contig_files.append(contigs_fasta)

    for seq_record in SeqIO.parse(contigs_fasta, 'fasta'):

        contig_label = '_'.join(
                [contig_label_base, seq_record.description])

        # Create an insertion model for the contig
        contig = Contig.objects.create(
                label=contig_label,
                parent_reference_genome=reference_genome,
                experiment_sample_to_alignment=(
                        sample_alignment))

        contig.metadata['coverage'] = float(
                seq_record.description.rsplit('_', 1)[1])
        contig.metadata['timestamp'] = timestamp
        contig.metadata['node_number'] = int(
                contig_number_pattern.findall(seq_record.description)[0])
        contig.metadata['assembly_dir'] = data_dir

        contig.ensure_model_data_dir_exists()
        dataset_path = os.path.join(contig.get_model_data_dir(),
                'fasta.fa')

        with open(dataset_path, 'w') as fh:
            SeqIO.write([seq_record], fh, 'fasta')

        add_dataset_to_entity(
                contig,
                'contig_fasta',
                Dataset.TYPE.REFERENCE_GENOME_FASTA,
                filesystem_location=dataset_path)

    return contig_files


def add_bam_track(reference_genome, bam_file, label):
    print 'add_bam_file_track entered'
    ag = AlignmentGroup.objects.create(
            reference_genome=reference_genome,
            label=label)
    es = ExperimentSample.objects.create(
            project=reference_genome.project,
            label=label)
    esta = ExperimentSampleToAlignment.objects.create(
            alignment_group=ag,
            experiment_sample=es)
    coordinate_sorted_bam = (os.path.splitext(bam_file)[0] +
            '.coordinate_sorted.bam')
    print 'created related models'
    sort_bam_by_coordinate(bam_file, coordinate_sorted_bam)
    print 'sorted bam'
    index_bam(coordinate_sorted_bam)
    print 'indexed bam'
    add_dataset_to_entity(esta,
                label,
                Dataset.TYPE.BWA_ALIGN,
                filesystem_location=coordinate_sorted_bam)
    print 'about to add bam file track'
    add_bam_file_track(reference_genome, esta,
            Dataset.TYPE.BWA_ALIGN)
    print 'added bam file track'
