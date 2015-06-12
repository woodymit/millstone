from Bio import SeqIO
import os
import sys

# Setup Django environment. -- DEBUG
sys.path.append(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), '../'))
os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'
    
from django.conf import settings

from main.model_utils import get_dataset_with_type
from main.models import Dataset
from main.models import AlignmentGroup
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import ReferenceGenome
from pipeline.read_alignment import align_with_bwa_mem
from utils import convert_seqrecord_to_fastq
from utils.import_util import add_dataset_to_entity

def find_insertion_region(reference_genome, contig_seqrecord):
    alignment_group = AlignmentGroup.objects.create(
            reference_genome=reference_genome,
            label='contig_alignment')

    contig_sample = ExperimentSample.objects.create(
            project=reference_genome.project,
            label='contig_sample')

    # Convert inserted sequence to fastq for alignment
    fastq_path = os.path.join(
            contig_sample.get_model_data_dir(), 'insertion.fq')
    convert_seqrecord_to_fastq(contig_seqrecord, fastq_path)

    # Add fastq dataset to Experiment Sample
    add_dataset_to_entity(
            contig_sample, 'contig_fastq', Dataset.TYPE.FASTQ1,
            filesystem_location=fastq_path)

    sample_to_alignment = ExperimentSampleToAlignment.objects.create(
            alignment_group=alignment_group,
            experiment_sample=contig_sample)

    align_with_bwa_mem(
            alignment_group,
            sample_to_alignment,
            project=reference_genome.project)

    alignment_bam_path = get_dataset_with_type(
            sample_to_alignment, Dataset.TYPE.BWA_ALIGN
                    ).get_absolute_location()

    return alignment_bam_path


def test():
    reference_genome = ReferenceGenome.objects.get(
            project__title='gf_testing_2',
            label='insertion_ref')

    contig_ref = ReferenceGenome.objects.get(
            project__title='gf_testing_2',
            label='insertion_ref_insertion_reads_de_novo_contigs')

    contig_fasta_path = get_dataset_with_type(
            contig_ref,
            Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

    with open(contig_fasta_path, 'r') as fh:
        contig_seqrecord = SeqIO.parse(fh, 'fasta').next()
        return find_insertion_region(reference_genome, contig_seqrecord)


if __name__ == '__main__':
    test()


# def place_insertion(reference_genome, contig_reference_genomes,
#       new_genome_label):

#   AlignmentGroup.objects.create(
#           reference_genome=reference_genome,
#           label='contig_alignment')
#   align_with_bwa_mem(
#                 experiment_sample_to_alignment.alignment_group,
#                 experiment_sample_to_alignment,
#                 project=reference_genome.project)
