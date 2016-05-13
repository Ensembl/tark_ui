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

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation

from tark.fields import ChecksumField, SequenceField
from tark.tark_exceptions import AssemblyNotFound, ReleaseNotFound,\
    FeatureNotFound
import hashlib
import pprint
from tark.lib.mapper import Mapper

FEATURE_TYPES = {'gene': 'Gene', 
                 'transcript': 'Transcript', 
                 'exon': 'Exon',
                 'translation': 'Translation'
                 }

FEATURE_LOOKUP = {'Gene': 1,
                  'Transcript': 2,
                  'Exon': 3,
                  'Translation': 4,
                  'Operon': 5}

SPECIES_NAMES = {'human': 'GCA_000001405'}

class FeatureQuerySet(models.query.QuerySet):
    def dict_iterator(self, **kwargs):          
        for feature in self.iterator():
            feature_obj = feature.to_dict(**kwargs)
            
            yield feature_obj
    
    def seq_iterator(self, format=None, **kwargs):
        format = format

        for feature in self.all():
            if not feature.has_sequence:
                continue
            
            if format:
                yield SeqRecord(feature.seq,
                                id=feature.stable_id_versioned,
                                description=feature.location).format(format)
            else:
                yield SeqRecord(feature.seq,
                                id=feature.stable_id_versioned,
                                description=feature.location)
        
    def to_dict(self, **kwargs):
        features_ary = []

        for feature in self.all():
            feature_obj = feature.to_dict(**kwargs)
            
            features_ary.append(feature_obj)
            
        return features_ary

    def release(self, release_set, exclude=False):
        feature_id = FEATURE_LOOKUP[self.model.__name__]
        set_table_name = 'release_tag_' + str(release_set.release_id)
        
        if exclude:
            release_tags = Releasetag.objects.filter(release_id=release_set, feature_type=feature_id)
            pk_column = self.model._meta.pk.name + '__in'
            return self.exclude(**{pk_column: release_tags})
        else:
            return self.extra( where=[set_table_name + ".feature_type=%s", 
                                      set_table_name + ".feature_id = " + self.model._meta.db_table + '.' + self.model._meta.pk.name, 
                                      set_table_name + ".release_id='%s'"], 
                                params=[feature_id, release_set.release_id], 
                                tables=['release_tag` as `' + set_table_name])

    def tag(self, tag_set, exclude=False):
        set_table_name = 'tag_' + str(tag_set.tagset_id)
        if exclude:
            settags = Tag.objects.filter(tagset_id=tag_set)
            return self.exclude(transcript_id__in=settags)
        else:
            return self.extra( where=[set_table_name + ".transcript_id = " + self.model._meta.db_table + '.' + self.model._meta.pk.name, 
                                      set_table_name + ".tagset_id='%s'"], 
                              params=[tag_set.tagset_id], 
                              tables=['tag` as `' + set_table_name] )

    def build_filters(self, **kwargs):
        filter = {}
        
        if 'assembly' in kwargs:
            if type(kwargs['assembly']) == Assembly:
                assembly = kwargs['assembly']
            else:
                assembly = Assembly.fetch_by_accession(kwargs['assembly'])

            filter.update({'assembly_id': assembly})

        if filter:
            return self.filter(**filter)
        else:
            return self

    def by_stable_id(self, stable_id, **kwargs):
        print "trying {}".format(stable_id)
        try:
            split_stable_id = stable_id.split('.')
            # Will this be a problem for null version stable_ids?
            if len(split_stable_id) == 2:
                feature = self.build_filters(**kwargs).filter(stable_id=split_stable_id[0], stable_id_version=split_stable_id[1])
                
                if feature:
                    return feature
        except Exception as e:
            print e        

        return self.build_filters(**kwargs).filter(stable_id=stable_id)     

    def fetch_by_name(self, name, **kwargs):
        if self.model.__name__ == 'Gene':
            return self.filter(genenames__name=name).build_filters(**kwargs)
        elif self.model.__name__ == 'Transcript':
            return self.filter(gene__genenames__name=name).build_filters(**kwargs)
        
        raise FeatureNotFound("Feature {} not found".format(name))
                    
    def checksum(self):
        checksums = []
        checksum_field = self.model._meta.db_table + '_checksum'
        
        for feature in self.order_by('pk').all():
            checksum = getattr(feature, checksum_field)
            checksums.append(checksum)
            
        set_checksum = hashlib.sha1(':'.join(checksums)).hexdigest()

        return set_checksum


class FeatureManager(models.Manager):
    def get_queryset(self):
        return FeatureQuerySet(self.model, using=self._db)
    
    def release(self, release_set):
        # Fetch the id of the feature type
        feature_id = FEATURE_LOOKUP[self.model.__name__]
        return self.get_queryset().release(release_set)
#        return self.get_queryset().extra( where=["release_tag.feature_type=%s", 
#                                                 "release_tag.feature_id = " + self.model._meta.db_table + '.' + self.model._meta.pk.name, 
#                                                 "release_tag.release_id='%s'"], 
#                                         params=[feature_id, release_set.release_id], 
#                                         tables=['release_tag'] )

    def tag(self, tag_set):
        # Fetch the id of the feature type
        return self.get_queryset.tag(tag_set)
#        return self.get_queryset().extra( where=["tag.transcript_id = " + self.model._meta.db_table + '.' + self.model._meta.pk.name, 
#                                                 "tag.tagset_id='%s'"], 
#                                         params=[tag_set.tagset_id], 
#                                         tables=['tag'] )

    def by_stable_id(self, stable_id, **kwargs):
        return self.get_queryset().by_stable_id(stable_id, **kwargs)

    def fetch_by_name(self, name, **kwargs):
        return self.get_queryset().fetch_by_name(name, **kwargs)
    
#    def dict_iterator(self, **kwargs):
#        for stable_id in self.stable_ids:
#            featureset = self.get_queryset().by_stable_id(stable_id)
#            for feature in featureset.to_dict( **kwargs ):
#                yield feature 

#    def seq_iterator(self, format=None, **kwargs):
#        format = format

#        for stable_id in self.stable_ids:
#            featureset = self.get_queryset().by_stable_id(stable_id, **kwargs)
#            for feature in featureset.all():
#                if not feature.has_sequence:
#                    continue
#            
#                if format:
#                    yield SeqRecord(feature.seq,
#                                    id=feature.stable_id_versioned,
#                                    description=feature.location).format(format)
#                else:
#                    yield SeqRecord(feature.seq,
#                                    id=feature.stable_id_versioned,
#                                    description=feature.location)

    def by_stable_ids(self, stable_ids, **kwargs):
        self.stable_ids = stable_ids
        return self
#        return self.get_queryset().by_stable_ids(stable_ids, **kwargs)


class Feature(models.Model):
    objects = FeatureManager()

    @property
    def location(self):
        loc = "{}:{}:{}:{}:{}".format(str(self.assembly), self.loc_region, self.loc_start, self.loc_end, self.loc_strand)
        
        return loc

    @property
    def stable_id_versioned(self):
        return str(self)
    
    @property
    def has_sequence(self):
        return True if hasattr(self, 'seq_checksum') else False
        
    @property
    def seq(self, **kwargs):
        seq = Sequence.objects.get(pk=self.seq_checksum)

        return seq.sequence

    @classmethod
    def build_filters(cls, **kwargs):
        filter = {}
        
        if 'assembly' in kwargs:
            if type(kwargs['assembly']) == Assembly:
                assembly = kwargs['assembly']
            else:
                assembly = Assembly.fetch_by_accession(kwargs['assembly'])
            filter.update({'assembly_id': assembly})
            
        return filter

    def dict_iterator(self, **kwargs):
        return self.to_dict(**kwargs)
            
    def to_dict(self, **kwargs):
        feature_obj = {}
        filter_pk = kwargs.get('filter_pk', True)
        for field in self._meta.fields:
            if field.name == 'session':
                continue
            elif (field.name == self._meta.pk.name) and filter_pk:
                continue
            elif field.get_internal_type() == "ForeignKey":
                feature_obj[field.name] = str(getattr(self, field.name))

                if field.name == 'seq_checksum' and not kwargs.get('skip_sequence', False):

                    feature_obj['seq_checksum'] = str(getattr(self, field.name))
                    seq = Sequence.fetch_sequence(feature_obj[field.name], **{'feature_type': type(self).__name__})
                    feature_obj['sequence'] = seq.sequence
            else:
                feature_obj[field.name] = getattr(self, field.name)
                
        return feature_obj

    @property
    def feature_type(self):
        return type(self).__name__.lower()
    
    def feature_location(self, start, end):
        return FeatureLocation(start, end, strand=self.loc_strand, ref=self.stable_id, ref_db=self.feature_type)

    def __contains__(self, item):
        if type(item) == FeatureLocation:
            if self.loc_strand != item.strand:
                raise Exception("Stands of feature {} does not match given location {}".format(self, item))
            
            # start < end, always, so strand doesn't matter
            if item.start >= self.loc_start and item.end <= self.loc_end:
                return True
            
        return False
                     

    def __str__(self):
        return "{}.{}".format(self.stable_id, self.stable_id_version) if self.stable_id_version else "{}".format(self.stable_id)

    class Meta:
        abstract = True

class VariationFeature():
    def __init__(self, feature=None, alt_seq=None):
        if not feature:
            raise Exception("No feature given")
        
        self.feature = feature
        if alt_seq:
            self.alt_seq = alt_seq

class Assembly(models.Model):
    assembly_id = models.AutoField(primary_key=True)
    genome = models.ForeignKey('Genome', models.DO_NOTHING, blank=True, null=True)
    assembly_name = models.CharField(max_length=128, blank=True, null=True)
    assembly_accession = models.CharField(max_length=32, blank=True, null=True)
    assembly_version = models.IntegerField(blank=True, null=True)
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    @classmethod
    def fetch_by_accession(cls, accession):
        try:            
            split_accession = accession.rsplit('.', 1)
            if len(split_accession) == 2:
                assembly = Assembly.objects.get(assembly_accession=split_accession[0], assembly_version=split_accession[1])
                return assembly
            
            assembly = Assembly.objects.filter(assembly_accession=accession)
            if assembly:
                return assembly
        except Exception as e:
            pass
            
        raise AssemblyNotFound("Accession " + accession + " not found")

    @classmethod
    def fetch_by_name(cls, name):
        accession = SPECIES_NAMES[name.lower()]
        
        if not accession:
            raise AssemblyNotFound("Assembly for " + name + " not found")

        return cls.fetch_by_accession(accession)

    def __str__(self):
        return "{}:{}.{}".format(self.assembly_name, self.assembly_accession, str(self.assembly_version))

    class Meta:
        managed = False
        db_table = 'assembly'
        unique_together = (('assembly_name', 'assembly_version'),)


class Exon(Feature):
    exon_id = models.AutoField(primary_key=True)
    stable_id = models.CharField(max_length=64)
    stable_id_version = models.IntegerField()
    assembly = models.ForeignKey(Assembly, models.DO_NOTHING, blank=True, null=True)
    loc_start = models.IntegerField(blank=True, null=True)
    loc_end = models.IntegerField(blank=True, null=True)
    loc_strand = models.IntegerField(blank=True, null=True)
    loc_region = models.CharField(max_length=42, blank=True, null=True)
    loc_checksum = ChecksumField(max_length=20, blank=True, null=True)
    exon_checksum = ChecksumField(unique=True, max_length=20, blank=True, null=True)
    seq_checksum = models.ForeignKey('Sequence', models.DO_NOTHING, db_column='seq_checksum', blank=True, null=True)
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

#    @classmethod
#    def fetch_set(cls, transcript_id):
#        exon_transcript.objects.filter(transcript)

    class Meta:
        managed = False
        db_table = 'exon'


class Exontranscript(models.Model):
    exon_transcript_id = models.AutoField(primary_key=True)
    transcript = models.ForeignKey('Transcript', models.DO_NOTHING, blank=True, null=True, related_name='exons')
    exon = models.ForeignKey(Exon, models.DO_NOTHING, blank=True, null=True, related_name='transcript')
    exon_order = models.SmallIntegerField(blank=True, null=True)
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'exon_transcript'


class Gene(Feature):
    gene_id = models.AutoField(primary_key=True)
    stable_id = models.CharField(max_length=64)
    stable_id_version = models.IntegerField()
    assembly = models.ForeignKey(Assembly, models.DO_NOTHING, blank=True, null=True)
    loc_start = models.IntegerField(blank=True, null=True)
    loc_end = models.IntegerField(blank=True, null=True)
    loc_strand = models.IntegerField(blank=True, null=True)
    loc_region = models.CharField(max_length=42, blank=True, null=True)
    gene_checksum = ChecksumField(unique=True, max_length=20, blank=True, null=True)
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    def to_dict(self, **kwargs):
        feature_obj = super(Gene, self).to_dict(**kwargs)
        
        if kwargs.get('expand', False):
            children = Transcript.objects.filter(gene_id=self.gene_id).to_dict(**kwargs)
            if children:
                feature_obj['translation'] = children
    
        return feature_obj

    @property
    def mapper(self):
        if not hasattr(self, '_mapper'):
            self._mapper = Mapper(self)
            
        return self._mapper

    class Meta:
        managed = False
        db_table = 'gene'


class Genenames(models.Model):
    gene_names_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=32, blank=True, null=True)
    gene = models.ForeignKey(Gene, models.DO_NOTHING, blank=True, null=True, related_name='genenames')
    assembly = models.ForeignKey(Assembly, models.DO_NOTHING, blank=True, null=True)
    source = models.CharField(max_length=32, blank=True, null=True)
    primary_id = models.IntegerField(blank=True, null=True)
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'gene_names'


class Genome(models.Model):
    genome_id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=128, blank=True, null=True)
    tax_id = models.IntegerField(blank=True, null=True)
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    def __str__(self):
        return "{} ({})".format(self.name, self.tax_id)

    class Meta:
        managed = False
        db_table = 'genome'
        unique_together = (('name', 'tax_id'),)


class Operon(Feature):
    operon_id = models.AutoField(primary_key=True)
    stable_id = models.CharField(max_length=64, blank=True, null=True)
    stable_id_version = models.IntegerField(blank=True, null=True)
    assembly = models.ForeignKey(Assembly, models.DO_NOTHING, blank=True, null=True)
    loc_start = models.IntegerField(blank=True, null=True)
    loc_end = models.IntegerField(blank=True, null=True)
    loc_strand = models.IntegerField(blank=True, null=True)
    loc_region = models.CharField(max_length=42, blank=True, null=True)
    operon_checksum = ChecksumField(max_length=20, blank=True, null=True)
    seq_checksum = models.ForeignKey('Sequence', models.DO_NOTHING, db_column='seq_checksum', blank=True, null=True)
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'operon'


class Operontranscript(models.Model):
    operon_transcript_id = models.AutoField(primary_key=True)
    stable_id = models.CharField(max_length=64, blank=True, null=True)
    stable_id_version = models.IntegerField(blank=True, null=True)
    operon = models.ForeignKey(Operon, models.DO_NOTHING, blank=True, null=True)
    transcript = models.ForeignKey('Transcript', models.DO_NOTHING, blank=True, null=True)
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'operon_transcript'


class Releaseset(models.Model):
    release_id = models.AutoField(primary_key=True)
    shortname = models.CharField(max_length=24, blank=True, null=True)
    description = models.CharField(max_length=256, blank=True, null=True)
    assembly = models.ForeignKey(Assembly, models.DO_NOTHING, blank=True, null=True)
    release_date = models.DateField(blank=True, null=True)
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)
    release_checksum = ChecksumField(max_length=20, blank=True, null=True)

    @classmethod
    def fetch_set(cls, **kwargs):
        filter = {}
        if 'assembly' in kwargs:
            assembly = Assembly.fetch_by_accession(kwargs['assembly'])
            filter['assembly'] = assembly
        
        if 'tag' in kwargs:
            filter['shortname'] = kwargs['tag']
            
        if 'description' in kwargs:
            filter['description__contains'] = kwargs['description']
   
        try:
            return Releaseset.objects.get(**filter)
        
        except Exception as e:
            raise ReleaseNotFound("Unique release not found, trying specifying one or more of assembly, tag or description")
        

    class Meta:
        managed = False
        db_table = 'release_set'
        unique_together = (('shortname', 'assembly'),)


class Releasetag(models.Model):
    feature_id = models.PositiveIntegerField(primary_key=True)
#    feature_type = models
#    feature_id = models.IntegerField()
    feature_type = models.IntegerField()
    release = models.ForeignKey(Releaseset, models.DO_NOTHING, related_name='features')
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'release_tag'
        unique_together = (('feature_id', 'feature_type', 'release'),)


class Sequence(models.Model):
    seq_checksum = ChecksumField(primary_key=True, max_length=20)
#    sequence = models.TextField(blank=True, null=True)
    sequence = SequenceField(blank=True, null=True)
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    @property
    def seq(self):
        return str(self.sequence)

    @classmethod
    def fetch_sequence(cls, checksum, **kwargs):
        seq = Sequence.objects.get(pk=checksum)
        
        if 'feature_type' in kwargs:
            seq.sequence.alphabet = SequenceField.alphabet_type(kwargs['feature_type'])

        return seq


    def __str__(self):
        return "{}".format(self.seq_checksum)

    class Meta:
        managed = False
        db_table = 'sequence'


class Session(models.Model):
    session_id = models.AutoField(primary_key=True)
    client_id = models.CharField(max_length=128, blank=True, null=True)
    start_date = models.DateTimeField(blank=True, null=True)
    status = models.CharField(max_length=45, blank=True, null=True)

    def __str__(self):
        return "Session:{}.Client: {}".format(self.session_id, self.client_id)

    class Meta:
        managed = False
        db_table = 'session'


class Tag(models.Model):
    transcript = models.ForeignKey('Transcript', models.DO_NOTHING, primary_key=True)
    tagset = models.ForeignKey('Tagset', models.DO_NOTHING, related_name='tags')
    session = models.ForeignKey(Session, models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'tag'
        unique_together = (('transcript', 'tagset'),)


class Tagset(models.Model):
    tagset_id = models.AutoField(primary_key=True)
    shortname = models.CharField(max_length=45, blank=True, null=True)
    description = models.CharField(max_length=255, blank=True, null=True)
    version = models.CharField(max_length=20, blank=True, null=True)
    is_current = models.IntegerField(blank=True, null=True)
    session = models.ForeignKey(Session, models.DO_NOTHING, blank=True, null=True)
    tagset_checksum = ChecksumField(max_length=20, blank=True, null=True)

    @classmethod
    def fetch_set(cls, **kwargs):
        filter = {}
                
        if 'tag' in kwargs:
            filter['shortname'] = kwargs['tag']

        if 'version' in kwargs:
            filter['version'] = kwargs['version']
            
        if 'description' in kwargs:
            filter['description__contains'] = kwargs['description']
   
        return Tagset.objects.get(**filter)

    class Meta:
        managed = False
        db_table = 'tagset'
        unique_together = (('shortname', 'version'),)


class Transcript(Feature):
    transcript_id = models.AutoField(primary_key=True)
    stable_id = models.CharField(max_length=64)
    stable_id_version = models.IntegerField()
    assembly = models.ForeignKey(Assembly, models.DO_NOTHING, blank=True, null=True)
    loc_start = models.IntegerField(blank=True, null=True)
    loc_end = models.IntegerField(blank=True, null=True)
    loc_strand = models.IntegerField(blank=True, null=True)
    loc_region = models.CharField(max_length=42, blank=True, null=True)
    loc_checksum = ChecksumField(max_length=20, blank=True, null=True)
    exon_set_checksum = ChecksumField(max_length=20, blank=True, null=True)
    transcript_checksum = ChecksumField(unique=True, max_length=20, blank=True, null=True)
    seq_checksum = models.ForeignKey(Sequence, models.DO_NOTHING, db_column='seq_checksum', blank=True, null=True)
    gene = models.ForeignKey(Gene, models.DO_NOTHING, blank=True, null=True, related_name='transcripts')
    session = models.ForeignKey(Session, models.DO_NOTHING, blank=True, null=True)

    def is_coding(self):
        if self.translations.exists():
            return True
        
        return False

    @property
    def mapper(self):
        return self.gene.mapper.transcript(self)

    def to_dict(self, **kwargs):
        feature_obj = super(Transcript, self).to_dict(**kwargs)
       
        if kwargs.get('expand', False):
            children = Translation.objects.filter(transcript_id=self.transcript_id).to_dict(**kwargs)
            if children:
                feature_obj['translation'] = children
            
            exons_ary = []
            for exon_transcript in self.exons.all().order_by('exon_order'):
#                print exon_transcript.exon.exon_id
                exon_obj = exon_transcript.exon.to_dict()
                exons_ary.append(exon_obj)
                        
            if exons_ary:
                feature_obj['exon'] = exons_ary
    
        return feature_obj

    class Meta:
        managed = False
        db_table = 'transcript'


class Translation(Feature):
    translation_id = models.AutoField(primary_key=True)
    stable_id = models.CharField(max_length=64)
    stable_id_version = models.IntegerField()
    assembly = models.ForeignKey(Assembly, models.DO_NOTHING, blank=True, null=True)
    transcript = models.ForeignKey(Transcript, models.DO_NOTHING, blank=True, null=True, related_name='translations')
    loc_start = models.IntegerField(blank=True, null=True)
    loc_end = models.IntegerField(blank=True, null=True)
    loc_strand = models.IntegerField(blank=True, null=True)
    loc_region = models.CharField(max_length=42, blank=True, null=True)
    loc_checksum = ChecksumField(max_length=20, blank=True, null=True)
    translation_checksum = ChecksumField(unique=True, max_length=20, blank=True, null=True)
    seq_checksum = models.ForeignKey(Sequence, models.DO_NOTHING, db_column='seq_checksum', blank=True, null=True)
    session = models.ForeignKey(Session, models.DO_NOTHING, blank=True, null=True)

    @property
    def mapper(self):
        return self.transcript.mapper

    class Meta:
        managed = False
        db_table = 'translation'
