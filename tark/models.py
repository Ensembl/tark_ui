# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Make sure each ForeignKey has `on_delete` set to the desired behavior.
#   * Remove `managed = False` lines if you wish to allow Django to create, modify, and delete the table
# Feel free to rename the models, but don't rename db_table values or field names.
from __future__ import unicode_literals

from django.db import models
from django.apps import apps

from tark.fields import ChecksumField, SequenceField

import pprint

FEATURE_CHILD = {
    'gene': 'transcript',
    'transcript': 'translation'
}

FEATURE_TYPES = ['gene', 'transcript', 'exon', 'translation']

class feature(models.Model):
    @classmethod
    def fetch_feature(cls, fetch_children=False, filter_pk=True, **kwargs):
        features = cls.objects.filter(**kwargs)
        
        if not features:
            return None
        
        features_ary = []
        
        for feature in features:
            feature_obj = feature.to_dict(filter_pk=filter_pk)
            
            if fetch_children:
                child_key, child_features = feature.fetch_children(filter_pk=filter_pk)
                if child_features:
                    feature_obj[child_key] = child_features
                    
                if cls.__name__ == 'transcript':
                    exons_ary = []
                    for exon_transcript in feature.exons.all().order_by('exon_order'):
                        exon_obj = exon_transcript.exon.to_dict()
                        exons_ary.append(exon_obj)
                        
                    if exons_ary:
                        feature_obj['exon'] = exons_ary
                    
        
            features_ary.append(feature_obj)

        return features_ary
    
    
    def to_dict(self, filter_pk=True):
        feature_obj = {}
        for field in self._meta.fields:
            if field.name == 'session':
                continue
            elif (field.name == self._meta.pk.name) and filter_pk:
                continue
            elif field.get_internal_type() == "ForeignKey":
                feature_obj[field.name] = str(getattr(self, field.name))

                if field.name == 'seq_checksum':
#                    seq = sequence.objects.get(pk=feature_obj[field.name])
                    seq = sequence.fetch_sequence(feature_obj[field.name], **{'feature_type': type(self).__name__})
                    feature_obj['sequence'] = seq.sequence
            else:
                feature_obj[field.name] = getattr(self, field.name)
                
        return feature_obj
        
    
    def fetch_children(self, fetch_children=True, filter_pk=True):
        if self.__class__.__name__ not in FEATURE_CHILD:
            return None, None
        
        child_name = FEATURE_CHILD[self.__class__.__name__]
        
        child_model = apps.get_model('tark', child_name)
        filter_name = "{}_id".format(self.__class__.__name__)
        return child_name, child_model.fetch_feature(fetch_children=fetch_children, filter_pk=filter_pk, **{ filter_name: getattr(self, filter_name) })

    def __str__(self):
        return "{}.{}".format(self.stable_id, self.stable_id_version) if self.stable_id_version else "{}".format(self.stable_id)

    class Meta:
        abstract = True


class assembly(models.Model):
    assembly_id = models.AutoField(primary_key=True)
    genome = models.ForeignKey('genome', models.DO_NOTHING, blank=True, null=True)
    assembly_name = models.CharField(max_length=128, blank=True, null=True)
    assembly_accession = models.CharField(max_length=32, blank=True, null=True)
    assembly_version = models.IntegerField(blank=True, null=True)
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    def __str__(self):
        return "{}:{}.{}".format(self.assembly_name, self.assembly_accession, str(self.assembly_version))

    class Meta:
        managed = False
        db_table = 'assembly'
        unique_together = (('assembly_name', 'assembly_version'),)


class exon(feature):
    exon_id = models.AutoField(primary_key=True)
    stable_id = models.CharField(max_length=64)
    stable_id_version = models.IntegerField()
    assembly = models.ForeignKey(assembly, models.DO_NOTHING, blank=True, null=True)
    loc_start = models.IntegerField(blank=True, null=True)
    loc_end = models.IntegerField(blank=True, null=True)
    loc_strand = models.IntegerField(blank=True, null=True)
    loc_region = models.CharField(max_length=42, blank=True, null=True)
    loc_checksum = ChecksumField(max_length=20, blank=True, null=True)
    exon_checksum = ChecksumField(unique=True, max_length=20, blank=True, null=True)
    seq_checksum = models.ForeignKey('sequence', models.DO_NOTHING, db_column='seq_checksum', blank=True, null=True)
    session = models.ForeignKey('session', models.DO_NOTHING, blank=True, null=True)

#    @classmethod
#    def fetch_set(cls, transcript_id):
#        exon_transcript.objects.filter(transcript)

    class Meta:
        managed = False
        db_table = 'exon'


class exontranscript(models.Model):
    exon_transcript_id = models.AutoField(primary_key=True)
    transcript = models.ForeignKey('transcript', models.DO_NOTHING, blank=True, null=True, related_name='exons')
    exon = models.ForeignKey(exon, models.DO_NOTHING, blank=True, null=True, related_name='transcript')
    exon_order = models.SmallIntegerField(blank=True, null=True)
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'exon_transcript'


class gene(feature):
    gene_id = models.AutoField(primary_key=True)
    stable_id = models.CharField(max_length=64)
    stable_id_version = models.IntegerField()
    assembly = models.ForeignKey(assembly, models.DO_NOTHING, blank=True, null=True)
    loc_start = models.IntegerField(blank=True, null=True)
    loc_end = models.IntegerField(blank=True, null=True)
    loc_strand = models.IntegerField(blank=True, null=True)
    loc_region = models.CharField(max_length=42, blank=True, null=True)
    gene_checksum = ChecksumField(unique=True, max_length=20, blank=True, null=True)
    session = models.ForeignKey('session', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'gene'


class genenames(models.Model):
    gene_names_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=32, blank=True, null=True)
    gene = models.ForeignKey(gene, models.DO_NOTHING, blank=True, null=True, related_name='gene_names')
    assembly = models.ForeignKey(assembly, models.DO_NOTHING, blank=True, null=True)
    source = models.CharField(max_length=32, blank=True, null=True)
    primary_id = models.IntegerField(blank=True, null=True)
    session = models.ForeignKey('session', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'gene_names'


class genome(models.Model):
    genome_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=128, blank=True, null=True)
    tax_id = models.IntegerField(blank=True, null=True)
    session = models.ForeignKey('session', models.DO_NOTHING, blank=True, null=True)

    def __str__(self):
        return "{} ({})".format(self.name, self.tax_id)

    class Meta:
        managed = False
        db_table = 'genome'
        unique_together = (('name', 'tax_id'),)


class operon(feature):
    operon_id = models.AutoField(primary_key=True)
    stable_id = models.CharField(max_length=64, blank=True, null=True)
    stable_id_version = models.IntegerField(blank=True, null=True)
    assembly = models.ForeignKey(assembly, models.DO_NOTHING, blank=True, null=True)
    loc_start = models.IntegerField(blank=True, null=True)
    loc_end = models.IntegerField(blank=True, null=True)
    loc_strand = models.IntegerField(blank=True, null=True)
    loc_region = models.CharField(max_length=42, blank=True, null=True)
    operon_checksum = ChecksumField(max_length=20, blank=True, null=True)
    seq_checksum = models.ForeignKey('sequence', models.DO_NOTHING, db_column='seq_checksum', blank=True, null=True)
    session = models.ForeignKey('session', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'operon'


class operontranscript(models.Model):
    operon_transcript_id = models.AutoField(primary_key=True)
    stable_id = models.CharField(max_length=64, blank=True, null=True)
    stable_id_version = models.IntegerField(blank=True, null=True)
    operon = models.ForeignKey(operon, models.DO_NOTHING, blank=True, null=True)
    transcript = models.ForeignKey('transcript', models.DO_NOTHING, blank=True, null=True)
    session = models.ForeignKey('session', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'operon_transcript'


class releaseset(models.Model):
    release_id = models.AutoField(primary_key=True)
    shortname = models.CharField(max_length=24, blank=True, null=True)
    description = models.CharField(max_length=256, blank=True, null=True)
    assembly = models.ForeignKey(assembly, models.DO_NOTHING, blank=True, null=True)
    release_date = models.DateField(blank=True, null=True)
    session = models.ForeignKey('session', models.DO_NOTHING, blank=True, null=True)
    release_checksum = ChecksumField(max_length=20, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'release_set'
        unique_together = (('shortname', 'assembly'),)


class releasetag(models.Model):
    feature_id = models.PositiveIntegerField()
    feature_type = models
    feature_id = models.IntegerField()
    feature_type = models.IntegerField()
    release = models.ForeignKey(releaseset, models.DO_NOTHING, related_name='features')
    session = models.ForeignKey('session', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'release_tag'
        unique_together = (('feature_id', 'feature_type', 'release'),)


class sequence(models.Model):
    seq_checksum = ChecksumField(primary_key=True, max_length=20)
#    sequence = models.TextField(blank=True, null=True)
    sequence = SequenceField(blank=True, null=True)
    session = models.ForeignKey('session', models.DO_NOTHING, blank=True, null=True)

    @property
    def seq(self):
        return str(self.sequence)

    @classmethod
    def fetch_sequence(cls, checksum, **kwargs):
        seq = sequence.objects.get(pk=checksum)
        
        if 'feature_type' in kwargs:
            seq.sequence.alphabet = SequenceField.alphabet_type(kwargs['feature_type'])
                
        return seq


    def __str__(self):
        return "{}".format(self.seq_checksum)

    class Meta:
        managed = False
        db_table = 'sequence'


class session(models.Model):
    session_id = models.AutoField(primary_key=True)
    client_id = models.CharField(max_length=128, blank=True, null=True)
    start_date = models.DateTimeField(blank=True, null=True)
    status = models.CharField(max_length=45, blank=True, null=True)

    def __str__(self):
        return "Session:{}.Client: {}".format(self.session_id, self.client_id)

    class Meta:
        managed = False
        db_table = 'session'


class tag(models.Model):
    transcript = models.ForeignKey('transcript', models.DO_NOTHING)
    tagset = models.ForeignKey('tagset', models.DO_NOTHING)
    session = models.ForeignKey(session, models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'tag'
        unique_together = (('transcript', 'tagset'),)


class tagset(models.Model):
    tagset_id = models.AutoField(primary_key=True)
    shortname = models.CharField(max_length=45, blank=True, null=True)
    description = models.CharField(max_length=255, blank=True, null=True)
    version = models.CharField(max_length=20, blank=True, null=True)
    is_current = models.IntegerField(blank=True, null=True)
    session = models.ForeignKey(session, models.DO_NOTHING, blank=True, null=True)
    tagset_checksum = ChecksumField(max_length=20, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'tagset'
        unique_together = (('shortname', 'version'),)


class transcript(feature):
    transcript_id = models.AutoField(primary_key=True)
    stable_id = models.CharField(max_length=64)
    stable_id_version = models.IntegerField()
    assembly = models.ForeignKey(assembly, models.DO_NOTHING, blank=True, null=True)
    loc_start = models.IntegerField(blank=True, null=True)
    loc_end = models.IntegerField(blank=True, null=True)
    loc_strand = models.IntegerField(blank=True, null=True)
    loc_region = models.CharField(max_length=42, blank=True, null=True)
    loc_checksum = ChecksumField(max_length=20, blank=True, null=True)
    exon_set_checksum = ChecksumField(max_length=20, blank=True, null=True)
    transcript_checksum = ChecksumField(unique=True, max_length=20, blank=True, null=True)
    seq_checksum = models.ForeignKey(sequence, models.DO_NOTHING, db_column='seq_checksum', blank=True, null=True)
    gene = models.ForeignKey(gene, models.DO_NOTHING, blank=True, null=True, related_name='transcripts')
    session = models.ForeignKey(session, models.DO_NOTHING, blank=True, null=True)

#    def to_dict(self, **kwargs):
#        feature_obj = super(transcript, self).to_dict(self, **kwargs)
#        
#        if kwargs.get('fetch_children', False):
#            exons_ary = []
#            for exon_transcript in feature.exons.all().order_by('exon_order'):
#                exon_obj = exon_transcript.exon.to_dict()
#                exons_ary.append(exon_obj)
#                        
#            if exons_ary:
#                feature_obj['exon'] = exons_ary
#    
#        return feature_obj

    class Meta:
        managed = False
        db_table = 'transcript'


class translation(feature):
    translation_id = models.AutoField(primary_key=True)
    stable_id = models.CharField(max_length=64)
    stable_id_version = models.IntegerField()
    assembly = models.ForeignKey(assembly, models.DO_NOTHING, blank=True, null=True)
    transcript = models.ForeignKey(transcript, models.DO_NOTHING, blank=True, null=True, related_name='translations')
    loc_start = models.IntegerField(blank=True, null=True)
    loc_end = models.IntegerField(blank=True, null=True)
    loc_strand = models.IntegerField(blank=True, null=True)
    loc_region = models.CharField(max_length=42, blank=True, null=True)
    loc_checksum = ChecksumField(max_length=20, blank=True, null=True)
    translation_checksum = ChecksumField(unique=True, max_length=20, blank=True, null=True)
    seq_checksum = models.ForeignKey(sequence, models.DO_NOTHING, db_column='seq_checksum', blank=True, null=True)
    session = models.ForeignKey(session, models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'translation'
