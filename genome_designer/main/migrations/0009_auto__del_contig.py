# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Deleting model 'Contig'
        db.delete_table(u'main_contig')


    def backwards(self, orm):
        # Adding model 'Contig'
        db.create_table(u'main_contig', (
            ('parent_reference_genome', self.gf('django.db.models.fields.related.ForeignKey')(related_name='+', to=orm['main.ReferenceGenome'])),
            ('experiment_sample_to_alignment', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['main.ExperimentSampleToAlignment'])),
            (u'referencegenome_ptr', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['main.ReferenceGenome'], unique=True, primary_key=True)),
        ))
        db.send_create_signal(u'main', ['Contig'])


    models = {
        u'auth.group': {
            'Meta': {'object_name': 'Group'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '80'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'})
        },
        u'auth.permission': {
            'Meta': {'ordering': "(u'content_type__app_label', u'content_type__model', u'codename')", 'unique_together': "((u'content_type', u'codename'),)", 'object_name': 'Permission'},
            'codename': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['contenttypes.ContentType']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'auth.user': {
            'Meta': {'object_name': 'User'},
            'date_joined': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'email': ('django.db.models.fields.EmailField', [], {'max_length': '75', 'blank': 'True'}),
            'first_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Group']", 'symmetrical': 'False', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'}),
            'username': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        u'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'main.alignmentgroup': {
            'Meta': {'object_name': 'AlignmentGroup'},
            'aligner': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'alignment_options': ('main.custom_fields.PostgresJsonField', [], {'default': '\'{"skip_het_only": false, "call_as_haploid": false}\''}),
            'dataset_set': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['main.Dataset']", 'null': 'True', 'blank': 'True'}),
            'end_time': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '256', 'blank': 'True'}),
            'reference_genome': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.ReferenceGenome']"}),
            'start_time': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'NOT_STARTED'", 'max_length': '40'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'3c5e28cf'", 'unique': 'True', 'max_length': '8'})
        },
        u'main.chromosome': {
            'Meta': {'object_name': 'Chromosome'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'num_bases': ('django.db.models.fields.BigIntegerField', [], {}),
            'reference_genome': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.ReferenceGenome']"}),
            'seqrecord_id': ('django.db.models.fields.CharField', [], {'default': "'chrom_1'", 'max_length': '256'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'5831f7c8'", 'unique': 'True', 'max_length': '8'})
        },
        u'main.dataset': {
            'Meta': {'object_name': 'Dataset'},
            'filesystem_idx_location': ('django.db.models.fields.CharField', [], {'max_length': '512', 'blank': 'True'}),
            'filesystem_location': ('django.db.models.fields.CharField', [], {'max_length': '512', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'status': ('django.db.models.fields.CharField', [], {'default': "'READY'", 'max_length': '40'}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '40'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'77fd6f04'", 'unique': 'True', 'max_length': '8'})
        },
        u'main.experimentsample': {
            'Meta': {'object_name': 'ExperimentSample'},
            'children': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'parents'", 'symmetrical': 'False', 'through': u"orm['main.ExperimentSampleRelation']", 'to': u"orm['main.ExperimentSample']"}),
            'data': ('main.custom_fields.PostgresJsonField', [], {}),
            'dataset_set': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['main.Dataset']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'project': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.Project']"}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'5074348b'", 'unique': 'True', 'max_length': '8'})
        },
        u'main.experimentsamplerelation': {
            'Meta': {'object_name': 'ExperimentSampleRelation'},
            'child': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'child_relationships'", 'to': u"orm['main.ExperimentSample']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'parent': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'parent_relationships'", 'to': u"orm['main.ExperimentSample']"}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'57fe9081'", 'unique': 'True', 'max_length': '8'})
        },
        u'main.experimentsampletoalignment': {
            'Meta': {'object_name': 'ExperimentSampleToAlignment'},
            'alignment_group': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.AlignmentGroup']"}),
            'dataset_set': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['main.Dataset']", 'null': 'True', 'blank': 'True'}),
            'experiment_sample': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.ExperimentSample']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'4cdc40c0'", 'unique': 'True', 'max_length': '8'})
        },
        u'main.project': {
            'Meta': {'object_name': 'Project'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'owner': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.UserProfile']"}),
            's3_backed': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'title': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'79997ce0'", 'unique': 'True', 'max_length': '8'})
        },
        u'main.referencegenome': {
            'Meta': {'object_name': 'ReferenceGenome'},
            'dataset_set': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['main.Dataset']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_materialized_variant_view_valid': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'metadata': ('main.custom_fields.PostgresJsonField', [], {}),
            'project': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.Project']"}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'093688ce'", 'unique': 'True', 'max_length': '8'}),
            'variant_key_map': ('main.custom_fields.PostgresJsonField', [], {})
        },
        u'main.region': {
            'Meta': {'object_name': 'Region'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'reference_genome': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.ReferenceGenome']"}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '40'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'53125830'", 'unique': 'True', 'max_length': '8'})
        },
        u'main.regioninterval': {
            'Meta': {'object_name': 'RegionInterval'},
            'end': ('django.db.models.fields.BigIntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'region': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.Region']"}),
            'start': ('django.db.models.fields.BigIntegerField', [], {})
        },
        u'main.s3file': {
            'Meta': {'object_name': 'S3File'},
            'bucket': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'created_at': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'key': ('django.db.models.fields.CharField', [], {'max_length': '200'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '200', 'null': 'True'})
        },
        u'main.savedvariantfilterquery': {
            'Meta': {'object_name': 'SavedVariantFilterQuery'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'owner': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.UserProfile']"}),
            'text': ('django.db.models.fields.TextField', [], {}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'745883ec'", 'unique': 'True', 'max_length': '8'})
        },
        u'main.userprofile': {
            'Meta': {'object_name': 'UserProfile'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'bd8f5678'", 'unique': 'True', 'max_length': '8'}),
            'user': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['auth.User']", 'unique': 'True'})
        },
        u'main.variant': {
            'Meta': {'object_name': 'Variant'},
            'chromosome': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.Chromosome']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'position': ('django.db.models.fields.BigIntegerField', [], {}),
            'ref_value': ('django.db.models.fields.TextField', [], {}),
            'reference_genome': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.ReferenceGenome']"}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '40'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'909e894e'", 'unique': 'True', 'max_length': '8'})
        },
        u'main.variantalternate': {
            'Meta': {'object_name': 'VariantAlternate'},
            'alt_value': ('django.db.models.fields.TextField', [], {}),
            'data': ('main.custom_fields.PostgresJsonField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_primary': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'c0412015'", 'unique': 'True', 'max_length': '8'}),
            'variant': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.Variant']", 'null': 'True'})
        },
        u'main.variantcallercommondata': {
            'Meta': {'object_name': 'VariantCallerCommonData'},
            'alignment_group': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.AlignmentGroup']"}),
            'data': ('main.custom_fields.PostgresJsonField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'source_dataset': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.Dataset']"}),
            'variant': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.Variant']"})
        },
        u'main.variantevidence': {
            'Meta': {'object_name': 'VariantEvidence'},
            'data': ('main.custom_fields.PostgresJsonField', [], {}),
            'experiment_sample': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.ExperimentSample']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'573bc825'", 'unique': 'True', 'max_length': '8'}),
            'variant_caller_common_data': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.VariantCallerCommonData']"}),
            'variantalternate_set': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['main.VariantAlternate']", 'symmetrical': 'False'})
        },
        u'main.variantset': {
            'Meta': {'object_name': 'VariantSet'},
            'dataset_set': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['main.Dataset']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'label': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'reference_genome': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.ReferenceGenome']"}),
            'uid': ('django.db.models.fields.CharField', [], {'default': "'d359360b'", 'unique': 'True', 'max_length': '8'}),
            'variants': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['main.Variant']", 'null': 'True', 'through': u"orm['main.VariantToVariantSet']", 'blank': 'True'})
        },
        u'main.varianttovariantset': {
            'Meta': {'object_name': 'VariantToVariantSet'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'sample_variant_set_association': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['main.ExperimentSample']", 'null': 'True', 'blank': 'True'}),
            'variant': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.Variant']"}),
            'variant_set': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['main.VariantSet']"})
        }
    }

    complete_apps = ['main']