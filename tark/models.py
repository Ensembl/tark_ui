from __future__ import unicode_literals

from django.db import models
from django.apps import apps
from django.db.models import Q
import operator
from django.conf import settings

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation
from Bio.Seq import Seq

from .fields import ChecksumField, SequenceField, HGNCField, GeneSetField
from .exceptions import FeatureNotFound, IncompatibleFeatureType, AssemblyNotFound, ReleaseNotFound, FeatureNonCoding
import hashlib
import copy
import pprint
from .lib.mapper import Mapper
from .lib.seqfetcher import SeqFetcher

from __builtin__ import int, True
from django.db.models.fields.related import ForeignKey

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

SPECIES_NAMES = {'human': 'homo_sapiens'}

class FeatureQuerySet(models.query.QuerySet):
    def dict_iterator(self, **kwargs):          
        for feature in self.iterator():
            feature_obj = feature.to_dict(**kwargs)
            
            yield feature_obj
    
    def seq_iterator(self, format=None, **kwargs):
        format = format

        for feature in self.all():
            if not feature.has_sequence:
                if settings.DEBUG:
                    print "no sequence: {}".format(feature)
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

    def release(self, release_set=None, exclude=False, **kwargs):
#        print "release_set: {}".format(release_set)
        # Allow chaining with an empty release
        if not release_set:
            return self

        if type(release_set) != Releaseset:
            release_set = Releaseset.fetch_set(tag=release_set, **kwargs)
        
        # Get the class for the release_tag for this particular feature type
        release_tag = self.model.releasetag()
        feature_id = FEATURE_LOOKUP[self.model.__name__]
        # We need to alias the table name in case we're joining with release tags multiple times
        set_table_name = 'release_tag_' + str(release_set.release_id)
        
        if exclude:
            release_tags = release_tag.objects.filter(release_id=release_set)
#            release_tags = Releasetag.objects.filter(release_id=release_set, feature_type=feature_id)
            pk_column = self.model._meta.pk.name + '__in'
            return self.exclude(**{pk_column: release_tags})
        else:
            release_tags = release_tag.objects.filter(release_id=release_set)
            pk_column = self.model._meta.pk.name + '__in'
            return self.filter(**{pk_column: release_tags}).annotate(_default_release=models.Value(release_set.shortname, output_field=models.CharField()))
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
                              tables=['tag` as `' + set_table_name] ).annotate(_default_tag=models.Value(tag_set, output_field=models.IntegerField()))

    def build_filters(self, **kwargs):
        filter = {}

        if 'assembly' in kwargs:
            if type(kwargs['assembly']) == Assembly:
                assembly = kwargs['assembly']
            else:
                assembly = Assembly.fetch_by_accession(kwargs['assembly'])

            filter.update({'assembly_id': assembly})

        elif 'assembly_id' in kwargs:
            filter.update({'assembly_id': kwargs['assembly_id']})

        if filter:
            return self.filter(**filter)
        else:
            return self

    def by_stable_id(self, stable_id, **kwargs):
#        print "trying {}".format(stable_id)
        try:
            split_stable_id = stable_id.split('.')
            # Will this be a problem for null version stable_ids?
            if len(split_stable_id) == 2:
                feature = self.build_filters(**kwargs).release(kwargs.get('release', None), **kwargs).filter(stable_id=split_stable_id[0], stable_id_version=split_stable_id[1])
                
                if feature:
                    return feature
        except Exception as e:
            print e        

        return self.build_filters(**kwargs).release(kwargs.get('release', None), **kwargs).filter(stable_id=stable_id)     

    def by_name(self, name, **kwargs):
        try:
            hgnc = Genenames.objects.get(name=name)
        except Exception as e:
            if settings.DEBUG:
                print str(e)
            return []
            
        if self.model.__name__ == 'Gene':
            return self.build_filters(**kwargs).filter(hgnc_id=hgnc.external_id).build_filters(**kwargs)
        elif self.model.__name__ == 'Transcript':
            return self.build_filters(**kwargs).filter(gene__hgnc_id=hgnc.external_id).build_filters(**kwargs)
        
        raise FeatureNotFound("Feature {} not found".format(name))

    def by_location(self, loc_region, location, **kwargs):
        return self.filter(Q(loc_region=loc_region) & 
                           ((Q(loc_start__lte=location.start) & Q(loc_end__gte=location.start)) |
                           (Q(loc_start__lte=location.end) & Q(loc_end__gte=location.end)) |
                           (Q(loc_start__gte=location.start) & Q(loc_end__lte=location.end)) )).build_filters(**kwargs)
                    
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
        #feature_id = FEATURE_LOOKUP[self.model.__name__]
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

    def by_name(self, name, **kwargs):
        return self.get_queryset().by_name(name, **kwargs)
    
    def by_location(self, loc_region, location, **kwargs):
        return self.get_queryset().by_location(loc_region, location, **kwargs)
    
    def by_stable_ids(self, stable_ids, **kwargs):
        self.stable_ids = stable_ids
        return self

class Feature(models.Model):
    objects = FeatureManager()

    @property
    def location(self):
        loc = "{}:{}:{}:{}:{}".format(str(self.assembly), self.loc_region, self.loc_start, self.loc_end, self.loc_strand)
        
        return loc
    
    @property
    def short_location(self):
        loc = "{}:{}:{}:{}".format(self.loc_region, self.loc_start, self.loc_end, self.loc_strand)

        return loc

    @property
    def length(self):
        return abs(self.loc_end - self.loc_start) + 1

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

    @property
    def sequence(self, **kwargs):
        if not self.has_sequence:
            return None

        return self.seq_checksum

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

#                    feature_obj['seq_checksum'] = str(getattr(self, field.name))
                    seq = Sequence.fetch_sequence(feature_obj[field.name], **{'feature_type': type(self).__name__})
                    feature_obj['sequence'] = seq.sequence
            elif getattr(self, field.name, None):
                '''
                We're looking if the object's value is not-None, not if it has one or not.
                '''
                feature_obj[field.name] = getattr(self, field.name)
                
        return feature_obj

    @property
    def release(self, release=None):
        return getattr(self, '_default_release', None)
    
    def default_release(self, release=None):
        if release:
            setattr(self, '_default_release', release)

        elif not hasattr(self, '_default_release'):
            release = int(self.releases.aggregate(default_release=models.Max('release__shortname'))['default_release'])
            setattr(self, '_default_release', release)

        return self

    def in_release(self, release):
        for r in self.releases.values_list('release__shortname', flat=True):
            if r == release:
                return True

        return False

    @property
    def release_tags(self):
        return [str(tag.release) for tag in self.releases.all()] or None
    
    @property
    def feature_type(self):
        return type(self).__name__.lower()
    
    def feature_location(self, start, end):
        ref_db = self.feature_type if self.feature_type != 'translation' else 'cdna'
        return FeatureLocation(start, end, strand=self.loc_strand, ref=self.stable_id, ref_db=ref_db)

    def genomic_location(self, start, end):
        return FeatureLocation(start, end, strand=self.loc_strand, ref=self.stable_id, ref_db='genomic')

    def remap2feature(self, location, **kwargs):
        return self.mapper.remap2feature(location, self, **kwargs)
    
    def remap2genomic(self, location):
        return self.mapper.remap2genomic(location)

    def genomic2cdna(self, location):
        if self.feature_type == 'transcript':
            return self.mapper.remap2feature(location,
                                             self.translation)
        elif self.feature_type == 'translation':
            return self.mapper.remap2feature(location, self)
        
        raise IncompatibleFeatureType()
    
    def genomic2cds(self, location):
        if self.feature_type == 'transcript':
            return self.mapper.remap2feature(location,
                                             self.translation,
                                             coordinates='cds')
        elif self.feature_type == 'translation':
            return self.mapper.remap2feature(location, self,
                                             coordinates='cds')
        
        raise IncompatibleFeatureType()

    @classmethod
    def pairs(self, featureset1, featureset2):
        feature_hash = {}
        feature_pairs = []
        unpaired_locations = {}

        for feature in featureset1:
            feature_hash[feature.stable_id] = {}
            feature_hash[feature.stable_id]['a'] = feature

        for feature in featureset2:
            if feature.stable_id not in feature_hash:
                feature_hash[feature.stable_id] = {}
            feature_hash[feature.stable_id]['b'] = feature

        for stable_id in feature_hash:
            a = feature_hash[stable_id]['a'] if 'a' in feature_hash[stable_id] else None
            b = feature_hash[stable_id]['b'] if 'b' in feature_hash[stable_id] else None

            if a is None:
                if b.loc_checksum not in unpaired_locations:
                    unpaired_locations[b.loc_checksum] = {}
                if 'b' not in unpaired_locations[b.loc_checksum]:
                    unpaired_locations[b.loc_checksum]['b'] = []
                unpaired_locations[b.loc_checksum]['b'].append(b)
            elif b is None:
                if a.loc_checksum not in unpaired_locations:
                    unpaired_locations[a.loc_checksum] = {}
                if 'a' not in unpaired_locations[a.loc_checksum]:
                    unpaired_locations[a.loc_checksum]['a'] = []
                unpaired_locations[a.loc_checksum]['a'].append(a)
            else:
                feature_pairs.append([stable_id,a,b])

        # Now we have to go through our unpairs locations and see if
        # any are a potential stable_id change
        for location in unpaired_locations:
            loc_features = unpaired_locations[location]

            if 'a' not in loc_features:
                for feature in loc_features['b']:
                    feature_pairs.append([feature.stable_id,None,feature])
            elif 'b' not in loc_features:
                for feature in loc_features['a']:
                    feature_pairs.append([feature.stable_id,feature,None])
            else:
                loc = loc_features['a'][0].short_location
                feature_pairs.append(['remapped', loc,
                                      [str(f) for f in loc_features['a']] if len(loc_features['a']) < 1 else str(loc_features['a'][0]),
                                      [str(f) for f in loc_features['b']] if len(loc_features['b']) < 1 else str(loc_features['b'][0])])

        return feature_pairs

    def __contains__(self, item):
        if type(item) == FeatureLocation:
#            if self.loc_strand != item.strand:
#                raise Exception("Strands of feature {} does not match given location {}".format(self, item))
            
            # start < end, always, so strand doesn't matter
            if item.start >= self.loc_start and item.end <= self.loc_end:
                return True
        elif isinstance(item, ( int, long )):
            if item >= self.loc_start and item <= self.loc_end:
                return True
        elif isinstance(item, slice):
            if item.start and item.start not in self:
                return False
            if item.stop and item.stop not in self:
                return False
                
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
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    @property
    def assembly_str(self):
        """ Fetch a definitive name for this assembly, for use in lookups like
            mod_faidx.
            ** For now we're just using assembly_name, we'll have to switch to
               something involving species name once we have more species
        """
        return self.assembly_name
    
#        return "{}_{}".format(self.genome.name, self.assembly_name)

    @classmethod
    def fetch_by_accession(cls, accession):
        """ Fetch a single assembly by accession (ie. GCA_000001405.14)
            Lookup values stored in the AssemblyAlias table/model
        """
        try:            
            assembly = AssemblyAlias.objects.get(alias=accession).assembly
            if assembly:
                return assembly
        except Exception as e:
            pass
            
        raise AssemblyNotFound("Accession " + accession + " not found")

    @classmethod
    def fetch_by_name(cls, name, assembly=None):
        full_name = SPECIES_NAMES[name.lower()]
        
        if not full_name:
            raise AssemblyNotFound("Assembly for " + name + " not found")

        if assembly:
            ''' There should be only one '''
            return cls.objects.get(genome__name=full_name, assembly_name=assembly)
            
        return cls.objects.filter(genome__name=full_name)

    @classmethod
    def fetch_by_id(cls, assembly_id):
        return cls.objects.get(pk=assembly_id)

    def __str__(self):
        return "{}".format(self.assembly_name)

    class Meta:
        managed = False
        db_table = 'assembly'
        unique_together = (('assembly_name', 'genome'),)

class AssemblyAlias(models.Model):
    assembly_alias_id = models.AutoField(primary_key=True)
    alias = models.CharField(max_length=64, blank=True, null=True)
    genome = models.ForeignKey('Genome', models.DO_NOTHING, blank=True, null=True)
    assembly = models.ForeignKey(Assembly, models.DO_NOTHING, blank=True, null=True, related_name='aliases')
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    @classmethod
    def fetch_by_accession(cls, accession):
        try:            
            split_accession = accession.rsplit('.', 1)
            if len(split_accession) == 2:
                assemblyalias = AssemblyAlias.objects.get(assembly_accession=split_accession[0], assembly_version=split_accession[1])
                return assemblyalias
            
            assemblyalias = AssemblyAlias.objects.filter(assembly_accession=accession)
            if assemblyalias:
                return assemblyalias
        except Exception as e:
            pass
            
        raise AssemblyNotFound("Accession " + accession + " not found")

    class Meta:
        managed = False
        db_table = 'assembly_alias'
        unique_together = (('genome', 'assembly'),)


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
    transcripts_m2m = models.ManyToManyField('Transcript', through='Exontranscript')

    @property
    def transcripts(self):
        release = self.release if self.release else self.default_release().release
        return self.transcripts_m2m.annotate(_default_release=models.Value(release, output_field=models.CharField()))


    def to_dict(self, **kwargs):
        feature_obj = super(Exon, self).to_dict(**kwargs)

        feature_obj['release'] = self.release if self.release else self.default_release().release
        
        return feature_obj

    @property
    def gene(self):
        # This feels a bit fragile, but there's really no other connection
        # to a gene
        for transcript in self.transcripts.all():
            return transcript.gene
    
    @property
    def mapper(self):
        return self.gene.mapper

    @classmethod
    def releasetag(cls):
        return ExonReleaseTag

    @classmethod
    def difference(cls, exon1, exon2, recursive=True):
        differences = {}

        if not exon1:
            differences['added'] = {'release': exon2.release}
            return differences

        if not exon2:
            differences['missing'] = {'release': exon1.release}
            return differences

        # Has the version changed
        if exon1.stable_id_version != exon2.stable_id_version:
            differences['version'] = {'base': exon1.stable_id_version, 'updated': exon2.stable_id_version}

        # The location has changed
        if exon1.loc_checksum != exon2.loc_checksum:
            differences['location'] = {'base': exon1.location, 'updated': exon2.location}

        # The sequence has changed
        if exon1.seq_checksum != exon2.seq_checksum:
            differences['sequence'] = {'base': exon1.seq, 'updated': exon2.seq}

        return differences

    # Does not support negative indexes for slicing
    
    def __getitem__(self, val):
        if not self.has_sequence:
            return Seq('')

        if val not in self:
            return Seq('')

        seq = self.seq
        offset = self.loc_start

        if isinstance(val, slice):                
            start = val.start - offset if val.start else 0
            end = val.stop - offset + 1 if val.stop else self.length
            return seq[start:end:val.step]
            
        elif isinstance(val, int):
            val = val - offset
            return seq[val]

        return Seq('')

#    @classmethod
#    def fetch_set(cls, transcript_id):
#        exon_transcript.objects.filter(transcript)

    class Meta:
        managed = False
        db_table = 'exon'    

class Exontranscript(models.Model):
    exon_transcript_id = models.AutoField(primary_key=True)
    transcript = models.ForeignKey('Transcript', models.DO_NOTHING, blank=True, null=True, related_name='exons_old')
    exon = models.ForeignKey(Exon, models.DO_NOTHING, blank=True, null=True, related_name='transcripts_old')
    exon_order = models.SmallIntegerField(blank=True, null=True)
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'exon_transcript'
        ordering = ('exon_order',)


class Gene(Feature):
    gene_id = models.AutoField(primary_key=True)
    stable_id = models.CharField(max_length=64)
    stable_id_version = models.IntegerField()
    assembly = models.ForeignKey(Assembly, models.DO_NOTHING, blank=True, null=True)
    loc_start = models.IntegerField(blank=True, null=True)
    loc_end = models.IntegerField(blank=True, null=True)
    loc_strand = models.IntegerField(blank=True, null=True)
    loc_region = models.CharField(max_length=42, blank=True, null=True)
    loc_checksum = ChecksumField(unique=True, max_length=20, blank=True, null=True)
    hgnc = HGNCField('Genenames', models.DO_NOTHING, to_field='external_id', blank=True, null=True)
    gene_checksum = ChecksumField(unique=True, max_length=20, blank=True, null=True)
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)
    transcripts_m2m = models.ManyToManyField('Transcript', through='Transcriptgene')

    @property
    def transcripts(self):
        release = self.release if self.release else self.default_release().release
        return self.transcripts_m2m.annotate(_default_release=models.Value(release, output_field=models.CharField()))

    def to_dict(self, **kwargs):
        feature_obj = super(Gene, self).to_dict(**kwargs)

        seq = self.seq if not kwargs.get('skip_sequence', False) else None
        if seq:
            feature_obj['sequence'] = seq

        if self.release:
            feature_obj['release'] = self.release
            return self.expand_dict(feature_obj, self.release, **kwargs)
        else:
            feature_objs = []
            for release in self.release_tags:
                release_obj = copy.deepcopy(feature_obj)
                release_obj['release'] = release
                feature_objs.append(self.expand_dict(release_obj, release, **kwargs)) 

            return feature_objs

    def expand_dict(self, feature_obj, release, **kwargs):
#        print "release: {}".format(release)
        
        if kwargs.get('expand', False):
            transcript_ary = []
            for transcript in self.transcripts.all():
#                print exon_transcript.exon.exon_id
                transcript_obj = transcript.to_dict(**kwargs)
#                transcript_obj = transcript_gene.transcript.default_release(release).to_dict(**kwargs)
                transcript_ary.append(transcript_obj)
                        
            if transcript_ary:
                feature_obj['transcript'] = transcript_ary
    
        return feature_obj

    @property
    def has_sequence(self):
        return SeqFetcher.has_sequence(self.assembly.assembly_str, self.loc_region, self.loc_end)

    @property
    def seq(self, **kwargs):
        seq = SeqFetcher.fetch(self.assembly.assembly_str, self.loc_region, self.loc_start, self.loc_end)

        return Seq(seq) if seq else Seq('')

    @property
    def tag(self):
        tag = getattr(self, '_default_tag', None)

        return tag

    @property
    def sequence(self, **kwargs):
        return None

    @classmethod
    def difference(cls, gene1, gene2, recursive=True):
        differences = {}

        if not gene1:
            differences['added'] = {'release': gene2.release}
            return differences

        if not gene2:
            differences['missing'] = {'release': gene1.release}
            return differences

        # Has the version changed
        if gene1.stable_id_version != gene2.stable_id_version:
            differences['version'] = {'base': gene1.stable_id_version, 'updated': gene2.stable_id_version}

        # The location has changed
        if gene1.loc_checksum != gene2.loc_checksum:
            differences['location'] = {'base': gene1.location, 'updated': gene2.location}

        if not recursive:
            return differences

        # On to the transcripts
        gene1_transcripts = gene1.transcripts.order_by('stable_id')
        gene2_transcripts = gene2.transcripts.order_by('stable_id')

        transcript_pairs = cls.pairs(gene1_transcripts, gene2_transcripts)
        transcript_differences = []

        for pair in transcript_pairs:
#            print "stable_id {}, set a: {}, set b: {}".format(pair[0],
#                                                              pair[1],
#                                                              pair[2])
            if len(pair) > 3:
                transcript_differences.append({pair[0]: pair[1:]})
            else:
                transcript_difference = Transcript.difference(pair[1],
                                                              pair[2])
                if transcript_difference:
                    transcript_difference['stable_id'] = pair[0]
                    transcript_differences.append(transcript_difference)

        if transcript_differences:
            differences['transcript_differences'] = transcript_differences

        return differences

    @property
    def mapper(self):
        if not hasattr(self, '_mapper'):
            self._mapper = Mapper(self)
            
        return self._mapper

    @classmethod
    def releasetag(cls):
        return GeneReleaseTag

    class Meta:
        managed = False
        db_table = 'gene'

class Genenames(models.Model):
    gene_names_id = models.AutoField(primary_key=True)
    external_id = models.IntegerField(max_length=32, blank=True, null=True)
    name = models.CharField(max_length=32, blank=True, null=True)
    source = models.CharField(max_length=32, blank=True, null=True)
    primary_id = models.IntegerField(blank=True, null=True)
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    def __str__(self):
        return "{}".format(self.name)

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
#        pprint.pprint(kwargs)
        filter = {}
        if 'assembly' in kwargs:
            assembly = Assembly.fetch_by_accession(kwargs['assembly'])
            filter['assembly'] = assembly
        elif 'assembly_id' in kwargs:
            filter['assembly'] = Assembly.fetch_by_id(kwargs['assembly_id'])
        
        if 'tag' in kwargs:
            filter['shortname'] = kwargs['tag']
            
        if 'description' in kwargs:
            filter['description__contains'] = kwargs['description']
   
        try:
            return Releaseset.objects.get(**filter)
        
        except Exception as e:
            raise ReleaseNotFound("Unique release not found, trying specifying one or more of assembly, tag or description")
        
    def __str__(self):
        return self.shortname

    class Meta:
        managed = False
        db_table = 'release_set'
        unique_together = (('shortname', 'assembly'),)


class Sequence(models.Model):
    seq_checksum = ChecksumField(primary_key=True, max_length=20)
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

    @classmethod
    def by_location(cls, assembly, region, location):
#        if SeqFetcher.has_sequence(assembly.assembly_str, self.loc_region, self.loc_end):
            pass

#    @property
#    def seq(self, **kwargs):
#        seq = SeqFetcher.fetch(self.assembly.assembly_str, self.loc_region, self.loc_start, self.loc_end)

#        return Seq(seq) if seq else Seq('')
    
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
    session = models.ForeignKey(Session, models.DO_NOTHING, blank=True, null=True)
    genes_m2m = models.ManyToManyField('Gene', through='Transcriptgene')
    exons_m2m = models.ManyToManyField('Exon', through='Exontranscript')
    translations_m2m = models.ManyToManyField('Translation', through='Translationtranscript')

    @property
    def genes(self):
        release = self.release if self.release else self.default_release().release
        return self.genes_m2m.annotate(_default_release=models.Value(release, output_field=models.CharField()))

    @property
    def exons(self):
        release = self.release if self.release else self.default_release().release
        return self.exons_m2m.annotate(_default_release=models.Value(release, output_field=models.CharField()))

    # There should only be one translation, if there's more than one
    # we're in trouble, we're just returning the first
    @property
    def translation(self):
        release = self.release if self.release else self.default_release().release
        for translation in self.translations_m2m.annotate(_default_release=models.Value(release, output_field=models.CharField())).all():
            return translation
        
        raise FeatureNonCoding()

    @property
    def is_coding(self):
        if self.translations_m2m.exists():
            return True
        
        return False

    @property
    def gene(self):
        release = self.release if self.release else self.default_release().release

        for gene in self.genes.all():
            if gene.in_release(release):
                return gene.default_release(release)
            
        return None

    @property
    def cdna_coding_start(self):
        if not hasattr(self, '_cdna_coding_start'):
            self._cdna_coding_start = self.mapper.cdna_coding_start()
            
        return self._cdna_coding_start

    @property
    def cdna_coding_end(self):
        if not hasattr(self, '_cdna_coding_end'):
            self._cdna_coding_end = self.mapper.cdna_coding_end()
            
        return self._cdna_coding_end

    @property
    def mapper(self):
        return self.gene.mapper.transcript(self)
    
    def genomic2cdna(self, location):
        return self.mapper.genomic2cdna(location)

    def genomic2cds(self, location):
        return self.mapper.genomic2cds(location)
    
    def genomic2pep(self, location):
        return self.mapper.genomic2pep(location)

    def to_dict(self, **kwargs):
        feature_obj = super(Transcript, self).to_dict(**kwargs)

        if self.release:
            feature_obj['release'] = self.release
            return self.expand_dict(feature_obj, self.release, **kwargs)
        else:
            feature_objs = []
            for release in self.release_tags:
                release_obj = copy.copy(feature_obj)
                release_obj['release'] = release
                feature_objs.append(self.expand_dict(release_obj, release, **kwargs)) 

        return feature_objs

    def expand_dict(self, feature_obj, release, **kwargs):
#        print "release: {}".format(release)
        
        if kwargs.get('expand', False):
            if self.is_coding:
                feature_obj['translation'] = self.translation.to_dict(**kwargs)
#            children = Translation.objects.filter(transcript_id=self.transcript_id).annotate(_default_release=models.Value(release, output_field=models.CharField())).to_dict(**kwargs)
#            if children:
#                feature_obj['translation'] = children

            exons_ary = []
            for exon in self.exons.all():
                exon_obj = exon.to_dict(**kwargs)
                exons_ary.append(exon_obj)

            if exons_ary:
                feature_obj['exon'] = exons_ary
    
        return feature_obj

    def subseq(self, start, end, introns=False):
        return self.mapper.transcript_subseq(start, end, introns)

    @classmethod
    def releasetag(cls):
        return TranscriptReleaseTag

    @classmethod
    def difference(cls, transcript1, transcript2, recursive=True):
        differences = {}

        if not transcript1:
            differences['added'] = {'release': transcript2.release}
            return differences

        if not transcript2:
            differences['missing'] = {'release': transcript1.release}
            return differences

        # Has the version changed
        if transcript1.stable_id_version != transcript2.stable_id_version:
            differences['version'] = {'base': transcript1.stable_id_version, 'updated': transcript2.stable_id_version}

        # The location has changed
        if transcript1.loc_checksum != transcript2.loc_checksum:
            differences['location'] = {'base': transcript1.location, 'updated': transcript2.location}

        # The sequence has changed
        if transcript1.seq_checksum != transcript2.seq_checksum:
            differences['sequence'] = {'base': transcript1.seq, 'updated': transcript2.seq}

        if not recursive:
            return differences

        # On to the exons
        if transcript1.exon_set_checksum != transcript2.exon_set_checksum:
            transcript1_exons = transcript1.exons.order_by('stable_id')
            transcript2_exons = transcript2.exons.order_by('stable_id')

            exon_pairs = cls.pairs(transcript1_exons, transcript2_exons)
            exon_differences = []

            for pair in exon_pairs:
#                print "stable_id {}, set a: {}, set b: {}".format(pair[0],
#                                                                  pair[1],
#                                                                  pair[2])
                if len(pair) > 3:
                    exon_differences.append({pair[0]: pair[1:]})
                else:
                    exon_difference = Exon.difference(pair[1], pair[2])
                    if exon_difference:
                        exon_difference['stable_id'] = pair[0]
                        exon_differences.append(exon_difference)

            if exon_differences:
                differences['exon_differences'] = exon_differences

        # One of them is coding, but are both?
        # Annoyingly complex logic, but needed because a transcript
        # could go from coding to non-coding, or vice versa
        if transcript1.is_coding or transcript2.is_coding:
            if transcript1.is_coding == transcript2.is_coding:
                # Both are coding, this is the good situation

                translation_differences = Translation.difference(transcript1.translation,
                                                                transcript2.translation)

                if translation_differences:
                    translation_differences['stable_id'] = transcript1.translation.stable_id if transcript1.translation else transcript2.translation.stable_id
                    differences['translation_differences'] = translation_differences

            elif transcript1.is_coding:
                translation = transcript1.translation
                differences['translation_lost'] = {'stable_id': translation.stable_id, 'release': translation.release}

            elif transcript2.is_coding:
                translation = transcript2.translation
                differences['translation_gained'] = {'stable_id': translation.stable_id, 'release': translation.release}

        return differences

    # Does not support negative indexes for slicing
    
    def __getitem__(self, val):
        if val not in self:
            return Seq('')

        if isinstance(val, slice):                
            start = val.start if val.start else self.loc_start
            end = val.stop if val.stop else self.loc_end
#            return seq[start:end:val.step]
            
        elif isinstance(val, int):
            start = val
            end = val
#            return seq[val]

        return self.subseq(start, end)

    class Meta:
        managed = False
        db_table = 'transcript'

class Transcriptgene(models.Model):
    gene_transcript_id = models.AutoField(primary_key=True)
    gene = models.ForeignKey(Gene, models.DO_NOTHING, blank=True, null=True, related_name='transcripts_old')
    transcript = models.ForeignKey('Transcript', models.DO_NOTHING, blank=True, null=True, related_name='genes_old')
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'transcript_gene'

class Translation(Feature):
    translation_id = models.AutoField(primary_key=True)
    stable_id = models.CharField(max_length=64)
    stable_id_version = models.IntegerField()
    assembly = models.ForeignKey(Assembly, models.DO_NOTHING, blank=True, null=True)
    loc_start = models.IntegerField(blank=True, null=True)
    loc_end = models.IntegerField(blank=True, null=True)
    loc_strand = models.IntegerField(blank=True, null=True)
    loc_region = models.CharField(max_length=42, blank=True, null=True)
    loc_checksum = ChecksumField(max_length=20, blank=True, null=True)
    translation_checksum = ChecksumField(unique=True, max_length=20, blank=True, null=True)
    seq_checksum = models.ForeignKey(Sequence, models.DO_NOTHING, db_column='seq_checksum', blank=True, null=True)
    session = models.ForeignKey(Session, models.DO_NOTHING, blank=True, null=True)
    transcripts_m2m = models.ManyToManyField('Transcript', through='Translationtranscript')

    @property
    def transcripts(self):
        release = self.release if self.release else self.default_release().release
        return self.transcripts_m2m.annotate(_default_release=models.Value(release, output_field=models.CharField()))

    def to_dict(self, **kwargs):
        feature_obj = super(Translation, self).to_dict(**kwargs)

        feature_obj['release'] = self.release if self.release else self.default_release().release
        
        return feature_obj

    @property
    def mapper(self):
        return self.transcript.mapper

    def genomic2cdna(self, location):
        return self.mapper.genomic2cdna(location)

    def genomic2cds(self, location):
        return self.mapper.genomic2cds(location)

    def genomic2pep(self, location):
        return self.mapper.genomic2pep(location)

    @classmethod
    def releasetag(cls):
        return TranslationReleaseTag

    @classmethod
    def difference(cls, translation1, translation2):
        differences = {}

        # Has the version changed
        if translation1.stable_id_version != translation2.stable_id_version:
            differences['version'] = {'base': translation1.stable_id_version, 'updated': translation2.stable_id_version}

        # The location has changed
        if translation1.loc_checksum != translation2.loc_checksum:
            differences['location'] = {'base': translation1.location, 'updated': translation2.location}

        # The sequence has changed
        if translation1.seq_checksum != translation2.seq_checksum:
            differences['sequence'] = {'base': translation1.seq, 'updated': translation2.seq}

        return differences

    class Meta:
        managed = False
        db_table = 'translation'

class Translationtranscript(models.Model):
    transcript_translation_id = models.AutoField(primary_key=True)
    transcript = models.ForeignKey(Transcript, models.DO_NOTHING, blank=True, null=True)
    translation = models.ForeignKey(Translation, models.DO_NOTHING, blank=True, null=True)
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'translation_transcript'

class Releasetag(models.Model):
#    feature_id = models.PositiveIntegerField(primary_key=True)
#    feature_type = models
#    feature_id = models.IntegerField()
    feature_type = models.IntegerField()
#    release = models.ForeignKey(Releaseset, models.DO_NOTHING, related_name='features')
    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
#        db_table = 'release_tag'
        abstract = True
#        unique_together = (('feature_id', 'feature_type', 'release'),)

class GeneReleaseTagManager(models.Manager):
    def get_queryset(self):
        return super(GeneReleaseTagManager, self).get_queryset().filter(feature_type=FEATURE_LOOKUP['Gene'])
    
class GeneReleaseTag(Releasetag):
    objects = GeneReleaseTagManager()
    feature = models.ForeignKey('Gene', to_field='gene_id', primary_key=True, related_name='releases')
#    feature_type = models.IntegerField()
    release = models.ForeignKey(Releaseset, models.DO_NOTHING, related_name='gene_features')
#    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)
    
    class Meta:
        managed = False
        db_table = 'release_tag'
        unique_together = (('feature', 'feature_type', 'release'),)
#        proxy = True

class TranscriptReleaseTagManager(models.Manager):
    def get_queryset(self):
        return super(TranscriptReleaseTagManager, self).get_queryset().filter(feature_type=FEATURE_LOOKUP['Transcript'])
    
class TranscriptReleaseTag(Releasetag):
    objects = TranscriptReleaseTagManager()
    feature = models.ForeignKey('Transcript', to_field='transcript_id', primary_key=True, related_name='releases')
#    feature_type = models.IntegerField()
    release = models.ForeignKey(Releaseset, models.DO_NOTHING, related_name='transcript_features')
#    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)
    
    class Meta:
        managed = False
        db_table = 'release_tag'
        unique_together = (('feature', 'feature_type', 'release'),)
#        proxy = True

class ExonReleaseTagManager(models.Manager):
    def get_queryset(self):
        return super(ExonReleaseTagManager, self).get_queryset().filter(feature_type=FEATURE_LOOKUP['Exon'])
    
class ExonReleaseTag(Releasetag):
    objects = ExonReleaseTagManager()
    feature = models.ForeignKey('Exon', to_field='exon_id', primary_key=True, related_name='releases')
#    feature_type = models.IntegerField()
    release = models.ForeignKey(Releaseset, models.DO_NOTHING, related_name='exon_features')
#    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)
    
    class Meta:
        managed = False
        db_table = 'release_tag'
        unique_together = (('feature', 'feature_type', 'release'),)
#        proxy = True

class TranslationReleaseTagManager(models.Manager):
    def get_queryset(self):
        return super(TranslationReleaseTagManager, self).get_queryset().filter(feature_type=FEATURE_LOOKUP['Translation'])
    
class TranslationReleaseTag(Releasetag):
    objects = TranslationReleaseTagManager()
    feature = models.ForeignKey('Translation', to_field='translation_id', primary_key=True, related_name='releases')
#    feature_type = models.IntegerField()
    release = models.ForeignKey(Releaseset, models.DO_NOTHING, related_name='translation_features')
#    session = models.ForeignKey('Session', models.DO_NOTHING, blank=True, null=True)
   
    class Meta:
        managed = False
        db_table = 'release_tag'
        unique_together = (('feature', 'feature_type', 'release'),)
#        proxy = True

class ReleaseDifference(models.Model):
    stable_id = models.CharField(max_length=20, primary_key=True)
    gene_count = models.IntegerField
    gene_set = GeneSetField(max_length=200)
    release1 = models.IntegerField
    assembly1 = models.IntegerField
    release2 = models.IntegerField
    assembly2 = models.IntegerField

    @property
    def is_different(self):
        return len(self.gene_set) < 2 or self.gene_set[0][3] != self.gene_set[1][3]

    def difference(self):
        if not self.is_different:
            return {}
        
        differences = {'stable_id': self.stable_id}

        gene1 = Gene.objects.by_stable_id(self.stable_id, assembly_id=self.assembly1, release=self.release1).first()
        gene2 = Gene.objects.by_stable_id(self.stable_id, assembly_id=self.assembly2, release=self.release2).first()

        differences.update(Gene.difference(gene1, gene2))


        return differences

    @classmethod
    def fetch_differences(cls, release1, release2, assembly1, assembly2=None):
        if not assembly2:
            assembly2 = assembly1

        for gene_set in ReleaseDifference.compare(release1, release2, assembly1, assembly2):
            if not gene_set.is_different:
                continue

            yield gene_set.difference()

    @classmethod
    def compare(cls, release1, release2, assembly1, assembly2=None):
        if not assembly2:
            assembly2 = assembly1

        query = ("select stable_id, "
                 "GROUP_CONCAT(CONCAT(stable_id_version, ',', rs.shortname, ',', rs.assembly_id, ',', gene_checksum) SEPARATOR '#!#!#') as gene_set, count(stable_id) as gene_count from gene left join release_tag r1 on gene.gene_id = r1.feature_id and r1.feature_type = 1 join release_set rs on rs.release_id = r1.release_id and ((rs.shortname = %s and rs.assembly_id = %s) or (rs.shortname = %s and rs.assembly_id = %s)) group by stable_id order by gene_count desc")
        
        return cls.objects.raw("SELECT stable_id, GROUP_CONCAT(CONCAT(stable_id_version, ',', rs.shortname, ',', rs.assembly_id, ',', gene_checksum) ORDER BY FIELD(rs.shortname, %s, %s) SEPARATOR '#!#!#') as gene_set, count(stable_id) as gene_count, %s as release1, %s as assembly1, %s as release2, %s as assembly2 from gene left join release_tag r1 on gene.gene_id = r1.feature_id and r1.feature_type = 1 join release_set rs on rs.release_id = r1.release_id and ((rs.shortname = %s and rs.assembly_id = %s) or (rs.shortname = %s and rs.assembly_id = %s)) group by stable_id order by gene_count desc", [release1, release2, release1, assembly1, release2, assembly2, release1, assembly1, release2, assembly2])
    
    class Meta:
        managed = False
        db_table = 'null'
