import datetime
import os
import re

from Bio import SeqIO

from genome_finish.millstone_de_novo_fns import add_paired_mates
from genome_finish.millstone_de_novo_fns import get_clipped_reads
from genome_finish.millstone_de_novo_fns import get_unmapped_reads
from genome_finish.millstone_de_novo_fns import get_split_reads
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
VELVET_KMER_LIST = [21]


def generate_contigs(experiment_sample_to_alignment, contig_label_base):

    # Grab reference genome fasta path
    reference_genome = (
            experiment_sample_to_alignment.alignment_group.reference_genome)
    ref_fasta_dataset = reference_genome.dataset_set.get_or_create(
            type=Dataset.TYPE.REFERENCE_GENOME_FASTA)[0]
    prepare_ref_genome_related_datasets(reference_genome, ref_fasta_dataset)

    # Make data_dir directory to house genome_finishing files
    genome_finishing_dir = os.path.join(
            experiment_sample_to_alignment.get_model_data_dir(),
            'genome_finishing')

    # Make data_dir directory if it does not exist
    if not os.path.exists(genome_finishing_dir):
        os.mkdir(genome_finishing_dir)

    data_dir = os.path.join(genome_finishing_dir, '0')
    data_dir_counter = 0
    while(os.path.exists(data_dir)):
        data_dir_counter += 1
        data_dir = os.path.join(genome_finishing_dir, str(data_dir_counter))
    os.mkdir(data_dir)

    # Retrieve bwa mem .bam alignment if exists otherwise generate it
    if not experiment_sample_to_alignment.dataset_set.filter(
            type=Dataset.TYPE.BWA_ALIGN).exists():
        add_dataset_to_entity(
                experiment_sample_to_alignment,
                'sample_alignment_for_assembly',
                Dataset.TYPE.BWA_ALIGN)
        align_with_bwa_mem(
                experiment_sample_to_alignment.alignment_group,
                experiment_sample_to_alignment,
                project=reference_genome.project)

    # Get a bam of sorted SV indicants with pairs
    sv_indicants_bam = get_sv_indicating_reads(experiment_sample_to_alignment)
            # input_sv_indicant_types={'clipped': False})

    # Find insertion metrics
    ins_length, ins_length_sd = get_insert_size_mean_and_stdev(
            experiment_sample_to_alignment)

    velvet_opts = {
        'velveth': {
            'shortPaired': ''
        },
        'velvetg': {
            'read_trkg': 'yes',
            'ins_length': ins_length,
            'ins_length_sd': ins_length_sd,
            'cov_cutoff': VELVET_COVERAGE_CUTOFF,
        }
    }

    # Perform velvet assembly
    contig_files = assemble_with_velvet(
            data_dir, velvet_opts, sv_indicants_bam,
            reference_genome, experiment_sample_to_alignment,
            contig_label_base)

    return contig_files


def get_sv_indicating_reads(sample_alignment, input_sv_indicant_classes={}):

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

    # Aggregate SV indicants
    SV_indicants_bam = compilation_prefix + '.bam'
    if os.path.exists(SV_indicants_bam):
        raise Exception(SV_indicants_bam + ' already exists')

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
    SV_indicants_with_pairs_bam = compilation_prefix + '.with_pairs.bam'
    make_bam(SV_indicants_with_pairs_sam, SV_indicants_with_pairs_bam)

    # Sort for velvet assembly
    print 'sorting by name'
    sort_bam_by_name(SV_indicants_with_pairs_bam)

    return SV_indicants_with_pairs_bam


def assemble_with_velvet(data_dir, velvet_opts, sv_indicants_bam,
        reference_genome, experiment_sample_to_alignment,
        contig_label_base, velvet_dir_prefix=''):

    timestamp = str(datetime.datetime.now())
    contig_number_pattern = re.compile('^NODE_(\d+)_')

    contig_files = []
    kmer_list = VELVET_KMER_LIST
    for kmer_length in kmer_list:
        # Set hash length argument for velveth
        velvet_opts['velveth']['hash_length'] = kmer_length

        # Run velvet assembly on SV indicants
        velvet_dir = os.path.join(data_dir, velvet_dir_prefix + 'velvet_k' + str(kmer_length))
        run_velvet(
                sv_indicants_bam,
                velvet_dir,
                velvet_opts)

        # Collect resulting contigs fasta
        contigs_fasta = os.path.join(velvet_dir, 'contigs.fa')
        contig_files.append(contigs_fasta)

        for seq_record in SeqIO.parse(contigs_fasta, 'fasta'):

            contig_label = '_'.join(
                    [contig_label_base, seq_record.description])

            # Create an insertion model for the contig
            contig = Contig.objects.create(
                    label=contig_label,
                    parent_reference_genome=reference_genome,
                    experiment_sample_to_alignment=(
                            experiment_sample_to_alignment))

            contig.metadata['coverage'] = float(
                    seq_record.description.rsplit('_', 1)[1])

            contig.metadata['timestamp'] = timestamp

            contig.metadata['node_number'] = int(contig_number_pattern.findall(
                    seq_record.description)[0])

            contig.metadata['assembly_dir'] = velvet_dir

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
