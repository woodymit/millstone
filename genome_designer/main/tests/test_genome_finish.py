"""
Tests for genome finishing features
"""
import os
import tempfile

from Bio import SeqIO
from django.conf import settings
from django.contrib.auth import authenticate
from django.contrib.auth.models import User
from django.http.request import HttpRequest
from django.test import Client
from django.test import TestCase

from genome_finish.insertion_placement import place_contig
from genome_finish.millstone_de_novo_fns import get_local_contig_placement
from genome_finish.millstone_de_novo_fns import get_match_counts
from main.model_utils import get_dataset_with_type
from main.models import AlignmentGroup
from main.models import Contig
from main.models import Dataset
from main.models import ExperimentSample
from main.models import ExperimentSampleToAlignment
from main.models import Project
import main.xhr_handlers as xhr_handlers
from pipeline.pipeline_runner import run_pipeline
from utils import convert_fasta_to_fastq
from utils import generate_safe_filename_prefix_from_label
from utils.bam_utils import minimal_bwa_align
from utils.import_util import add_dataset_to_entity
from utils.import_util import import_reference_genome_from_local_file



TEST_USERNAME = 'testuser'
TEST_PASSWORD = 'password'
TEST_EMAIL = 'test@example.com'
TEST_PROJECT_NAME = 'testModels_project'
TEST_REF_GENOME_NAME = 'mg1655_partial'
TEST_REF_GENOME_PATH = os.path.join(
        settings.PWD,
        'test_data/full_vcf_test_set/mg1655_tolC_through_zupT.gb')

TEST_FASTA_1_PATH = os.path.join(
    settings.PWD,
    'test_data/genome_finish_test/random_fasta_1.fa')
TEST_FASTA_2_PATH = os.path.join(
    settings.PWD,
    'test_data/genome_finish_test/random_fasta_2.fa')
INS_1KB_REF_GENOME_PATH = os.path.join(
        settings.PWD, 'test_data/genome_finish_test/ins_1kb.fa')
INS_1KB_FQ_1_PATH = os.path.join(
        settings.PWD, 'test_data/genome_finish_test/ins_1kb.1.fq')
INS_1KB_FQ_2_PATH = os.path.join(
        settings.PWD, 'test_data/genome_finish_test/ins_1kb.2.fq')
INS_1KB_INSERTION_SEQUENCE_PATH = os.path.join(
        settings.PWD,
        'test_data/genome_finish_test/ins_1kb_insertion.fa')
INS_1KB_CONTIG_FASTA_PATH = os.path.join(
        settings.PWD,
        'test_data/genome_finish_test/ins_1kb_contig.fa')
INS_1KB_TRANSFORMED_FASTA_PATH = os.path.join(
        settings.PWD,
        'test_data/genome_finish_test/ins_1kb_transformed.fa')
INSERTION_LENGTH = 1000


class TestContigAssembly(TestCase):

    def setUp(self):
        # Useful models.
        self.user = User.objects.create_user(
            TEST_USERNAME, password=TEST_PASSWORD, email=TEST_EMAIL)
        self.project = Project.objects.create(
            owner=self.user.get_profile(), title='Test Project')

        # Fake web browser client used to make requests.
        self.client = Client()
        self.client.login(username=TEST_USERNAME, password=TEST_PASSWORD)

    def _perform_assembly(self, read_fastqs, reference_genome_path,
                          reference_genome_format=None):
        file_extension_to_genome_format = {
            '.fa': 'fasta',
            '.fasta': 'fasta',
            '.gb': 'genbank',
            '.genbank': 'genbank',
            '.gbk': 'genbank'
        }
        if reference_genome_format is None:
            reference_genome_format = file_extension_to_genome_format[
                    os.path.splitext(reference_genome_path)[1]]

        reference_genome = import_reference_genome_from_local_file(
                self.project, 'test_ref',
                reference_genome_path, reference_genome_format)

        # Make ExperimentSample
        experiment_sample = ExperimentSample.objects.create(
                project=self.project,
                label='test_sample')

        add_dataset_to_entity(
                experiment_sample, 'reads_fq_1', Dataset.TYPE.FASTQ1,
                filesystem_location=read_fastqs[0])
        if len(read_fastqs) == 2:
            add_dataset_to_entity(
                    experiment_sample, 'reads_fq_2', Dataset.TYPE.FASTQ2,
                    filesystem_location=read_fastqs[1])

        # Make Alignment group
        alignment_group = AlignmentGroup.objects.create(
            reference_genome=reference_genome)

        # Make resulting ExperimentSampleToAlignment
        experiment_sample_to_alignment = (
                ExperimentSampleToAlignment.objects.create(
                        alignment_group=alignment_group,
                        experiment_sample=experiment_sample))

        request_data = {
                'sampleAlignmentUid': experiment_sample_to_alignment.uid
        }

        request = HttpRequest()
        request.GET = request_data
        request.method = 'GET'
        request.user = self.user

        authenticate(username=TEST_USERNAME, password=TEST_PASSWORD)
        self.assertTrue(request.user.is_authenticated())

        initial_contig_uids = [c.uid for c in Contig.objects.filter(
                parent_reference_genome=reference_genome,
                experiment_sample_to_alignment=experiment_sample_to_alignment
        )]

        xhr_handlers.generate_contigs(request)

        contigs = Contig.objects.filter(
                parent_reference_genome=reference_genome,
                experiment_sample_to_alignment=experiment_sample_to_alignment
                        ).exclude(uid__in=initial_contig_uids)

        # Assert contigs were generated
        for c in contigs:
            self.assertTrue(c.num_bases > 0)

        return contigs

    def test_generate_contigs_with_existing_alignment(self):
        """Tests contig generation on samples with preexisting alignment
        """
        # Make ReferenceGenome
        ins_1kb_ref = import_reference_genome_from_local_file(
                self.project, 'ins_1kb', INS_1KB_REF_GENOME_PATH,
                'fasta', move=False)

        # Make ExperimentSample
        ins_1kb_reads = ExperimentSample.objects.create(
            project=self.project,
            label='ins_1kb_reads')
        add_dataset_to_entity(
                ins_1kb_reads, 'reads_fq_1', Dataset.TYPE.FASTQ1,
                filesystem_location=INS_1KB_FQ_1_PATH)
        add_dataset_to_entity(
                ins_1kb_reads, 'reads_fq_2', Dataset.TYPE.FASTQ2,
                filesystem_location=INS_1KB_FQ_2_PATH)

        # Run alignment of reads to reference
        alignment_group_label = 'fq_reads_pipeline_alignment'
        ref_genome = ins_1kb_ref
        sample_list = [ins_1kb_reads]
        alignment_group, _, _ = run_pipeline(
                alignment_group_label, ref_genome, sample_list,
                perform_variant_calling=False, alignment_options={})

        # Get resulting ExperimentSampleToAlignment
        reads_align = ExperimentSampleToAlignment.objects.get(
            alignment_group=alignment_group,
            experiment_sample=ins_1kb_reads)

        request_data = {
            'sampleAlignmentUid': reads_align.uid
        }

        print 'request_data:', request_data

        request = HttpRequest()
        request.GET = request_data
        request.method = 'GET'
        request.user = self.user

        authenticate(username=TEST_USERNAME, password=TEST_PASSWORD)
        self.assertTrue(request.user.is_authenticated())

        initial_contig_uids = [c.uid for c in Contig.objects.filter(
                parent_reference_genome=ins_1kb_ref,
                experiment_sample_to_alignment=reads_align)]

        xhr_handlers.generate_contigs(request)

        contigs = Contig.objects.filter(
                parent_reference_genome=ins_1kb_ref,
                experiment_sample_to_alignment=reads_align).exclude(
                        uid__in=initial_contig_uids)

        # Assert contigs were generated
        self.assertEqual(contigs.count(), 1)
        self.assertTrue(contigs[0].num_bases > 0)

    def test_1kb_insertion_detection(self):
        """ Tests ability to detect and assemble 1kb insertion
        """

        contigs = self._perform_assembly(
                [INS_1KB_FQ_1_PATH, INS_1KB_FQ_2_PATH],
                INS_1KB_REF_GENOME_PATH)

        # Only one contig, corresponding to the single 1kb insertion, should
        # be created
        self.assertTrue(contigs.count() == 1)
        contig = contigs[0]

        # Create Experiment Sample for the inserted sequence
        insertion_sequence_sample = ExperimentSample.objects.create(
                project=self.project, label='insertion_sequence')

        # Convert inserted sequence to fastq for alignment
        fastq_path = os.path.join(
                insertion_sequence_sample.get_model_data_dir(),
                'insertion.fq')
        convert_fasta_to_fastq(INS_1KB_INSERTION_SEQUENCE_PATH, fastq_path)

        # Add fastq dataset to Experiment Sample
        add_dataset_to_entity(
                insertion_sequence_sample, 'insertion_fastq',
                Dataset.TYPE.FASTQ1, filesystem_location=fastq_path)

        # Select insertion sequence fastq and contigs fasta for bwa alignment
        reads = [get_dataset_with_type(
                insertion_sequence_sample, Dataset.TYPE.FASTQ1
                ).get_absolute_location()]
        ref_fasta = get_dataset_with_type(
                contig, Dataset.TYPE.REFERENCE_GENOME_FASTA
                ).get_absolute_location()
        data_dir = os.path.join(
                contig.get_model_data_dir(), 'alignment_data')
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)
        alignment_bam = minimal_bwa_align(reads, ref_fasta, data_dir)

        # Find the maximum number of matches to the insertion captured by
        # any contig
        match_counts = get_match_counts(alignment_bam)
        max_matches = 0
        for k in match_counts:
            if max(match_counts[k]) > max_matches:
                max_matches = max(match_counts[k])

        # Ensure thorough coverage of the insertion
        INSERTION_COVERAGE_CUTOFF = 0.95
        max_cov_fraction = float(max_matches)/INSERTION_LENGTH
        print 'Contigs covered %.3f of the %d base inserted sequence\n' % (
                max_cov_fraction, INSERTION_LENGTH)
        assert max_cov_fraction >= INSERTION_COVERAGE_CUTOFF, (
                'The maximum fraction of the insertion captured by any ' +
                'contig was: %.3f, less than the passing cutoff: %.3f ' % (
                        max_cov_fraction, INSERTION_COVERAGE_CUTOFF))


class TestContigPlacement(TestCase):

    def setUp(self):
        # Useful models.
        self.user = User.objects.create_user(
            TEST_USERNAME, password=TEST_PASSWORD, email=TEST_EMAIL)
        self.project = Project.objects.create(
            owner=self.user.get_profile(), title='Test Project')

        # Fake web browser client used to make requests.
        self.client = Client()
        self.client.login(username=TEST_USERNAME, password=TEST_PASSWORD)

    def _make_temp_file(self, label, extension):
        if not os.path.exists(settings.TEMP_FILE_ROOT):
            os.mkdir(settings.TEMP_FILE_ROOT)
        _, temp_file_path = tempfile.mkstemp(
                suffix=('_' + generate_safe_filename_prefix_from_label(label) +
                        extension),
                dir=settings.TEMP_FILE_ROOT)
        return temp_file_path

    def _seqan_place_contig_string_input(self, ref_string, contig_string):

        ref_temp_fasta = self._make_temp_file('ref', '.fa')
        contig_temp_fasta = self._make_temp_file('contig', '.fa')

        with open(ref_temp_fasta, 'w') as fh:
            fh.write('\n'.join(['>ref', ref_string]))

        with open(contig_temp_fasta, 'w') as fh:
            fh.write('\n'.join(['>contig', contig_string]))

        return get_local_contig_placement(
                ref_temp_fasta, contig_temp_fasta)

    def test_seqan_place_contig(self):

        ref_string =    'AGCATGTTAGATAAGATAG'    +    'CTGTGCTAGTAGGCAGTCAGCGCCAT'
        contig_string = 'AGCATGTTAGATAAGATAGCCCCCCCCCCCCTGTGCTAGTAGGCAGTCAGCGCCAT'

        local_contig_placement = self._seqan_place_contig_string_input(
                ref_string, contig_string)

        self.assertEqual(local_contig_placement['contig_start_pos'], 19)
        self.assertEqual(local_contig_placement['contig_end_pos'], 30)

    def test_1kb_insertion_placement(self):
        new_reference_genome_label = 'insertion_incorporated'

        reference_genome = import_reference_genome_from_local_file(
                self.project, 'test_ref',
                INS_1KB_REF_GENOME_PATH, 'fasta')

        with open(INS_1KB_CONTIG_FASTA_PATH, 'r') as fh:
            contig_seqrecord = SeqIO.parse(fh, 'fasta').next()

        placed_contig_ref_genome = place_contig(
                reference_genome, contig_seqrecord,
                new_reference_genome_label)

        # Verify expected transformation
        placed_contig_fasta = get_dataset_with_type(
                placed_contig_ref_genome,
                Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

        with open(placed_contig_fasta, 'r') as fh:
            placed_contig_seqrecord = SeqIO.parse(fh, 'fasta').next()

        with open(INS_1KB_TRANSFORMED_FASTA_PATH, 'r') as fh:
            transformed_seqrecord = SeqIO.parse(fh, 'fasta').next()

        self.assertEqual(str(placed_contig_seqrecord.seq),
                str(transformed_seqrecord.seq))
