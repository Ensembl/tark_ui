from Bio.SeqFeature import FeatureLocation
import pprint

class Mapper():
    def __init__(self, gene):
        self._gene = gene
        print "gene start:end {}:{}".format(gene.loc_start, gene.loc_end)

    
    def transcript(self, transcript=None):
        if transcript:
            self._transcript = transcript
            return self
        
        if not hasattr(self, '_transcript'):
            raise MissingTranscriptError()
        
        return self._transcript

    def is_exonic(self, location):
        if location.ref_db and location.ref_db != 'genomic':
            location = self.remap2genomic(location)

        return True if self.fetch_exon(location) else False
       
    def fetch_exon(self, location):
        if location.ref_db and location.ref_db != 'genomic':
            location = self.remap2genomic(location)
            
        pprint.pprint(location)
        for exon in self.exons:
            if location in exon:
                return exon
            
        return None

    def find_exon(self, stable_id):
        for exon in self.exons:
            if exon.stable_id == stable_id:
                return exon
            
        return None

    def find_translation(self, stable_id):
        for translation in self.translations:
            if translation.stable_id == stable_id:
                return translation
            
        return None

    @property
    def exons(self):
        if not hasattr(self, '_transcript'):
            raise MissingTranscriptError()
        
        if not hasattr(self, '_exons'):
            self._exons = []
            for exon_transcript in self.transcript().exons.all().order_by('exon_order'):
                exon = exon_transcript.exon
                self.exons.append(exon)

        return self._exons
    
    @property
    def translations(self):
        if not hasattr(self, '_transcript'):
            raise MissingTranscriptError()

        if not hasattr(self, '_translations'):
            self._translations = []
            if not self.transcript().is_coding():
                raise FeatureNonCoding()
            
            for translation in self.transcript().translations.all():
                self._translations.append(translation)
                
        return self._translations
        

    def remap2genomic(self, location):
        # If no type is given we assume it's already genomic
        if not location.ref_db:
            return FeatureLocation(location.start, location.end, strand=location.strand, ref_db='genomic')

        if location.ref_db == 'gene':
            feature = self._gene
        elif location.ref_db == 'transcript':
            feature = self.transcript()
        elif location.ref_db == 'exon' and location.ref:
            feature = self.find_exon(location.ref)
            if not feature:
                raise FeatureNotFound()

        elif location.ref_db == 'translation' and location.ref:
            feature = self.find_translation(location.ref)
            if not feature:
                raise FeatureNonCoding()
        else:
            raise UnknownCoordinateSystem()
            
        if location.strand == 1:
            start = location.start + feature.loc_start
            end = location.end + feature.loc_start
        else:
            start = feature.loc_end - location.start
            end = feature.loc_end - location.end
            
        return FeatureLocation(start, end, strand=location.strand, ref='genomic')
            
    def remap2feature(self, location, feature):
        if location.strand != feature.loc_strand:
            raise WrongStrandError()
        
        if location.ref_db != 'genomic':
            location = self.remap2genomic(location)
            
        print "feature start:end {}:{}".format(feature.loc_start, feature.loc_end)
        if location.strand == 1:
            start = location.start - feature.loc_start
            end = location.end - feature.loc_start
        else:
            start = feature.loc_start - location.start
            end = feature.loc_end - location.end
            
        return FeatureLocation(start, end, strand=location.strand, ref=feature.stable_id, ref_db=feature.feature_type)

            

class TranscriptMapper():
    
    def __init__(self, transcript):
        self.transcript = transcript
        self.gene_start = transcript.gene.loc_start
        self.gene_end = transcript.gene.loc_end
        print "gene start:end {}:{}".format(self.gene_start, self.gene_end)
        print "tran start:end {}:{}".format(transcript.loc_start, transcript.loc_end)
        
        self.build_exons()

    def is_exonic(self, location):
        if location.ref and location.ref == 'feature':
            location = self.feature2genomic(location)
                    
#    def is_exonic(self, start, end, coordinate='genomic'):
#        if coordinate == 'feature':
#            location = self.feature2genomic(start, end, self.transcript.loc_strand)
            
        return True if self.fetch_exon(location) else False

    def fetch_exon(self, location):
        if location.ref and location.ref == 'feature':
            location = self.feature2genomic(location)
        pprint.pprint(location)
        for exon in self.exons:
            if location in exon:
                return exon

        return None
    
    def feature2genomic(self, location):
        if location.strand == 1:
            start = location.start + self.gene_start
            end = location.end + self.gene_start
        else:
            start = self.gene_end - location.start
            end = self.gene_end - location.end
            
        return FeatureLocation(start, end, strand=location.strand, ref='genomic')

#    def feature2genomic(self, start, end, strand=1):
#        if strand == 1:
#            start = start + self.gene_start - 1
#            end = end + self.gene_start - 1
#        else:
#            start = self.gene_end - start
#            end = self.gene_end - end
#
#        return FeatureLocation(start, end, strand=strand)
            
    def build_exons(self):
        self.exons = []
        for exon_transcript in self.transcript.exons.all().order_by('exon_order'):
            exon = exon_transcript.exon
            self.exons.append(exon)
            
    
class ExonMapper():
    def __init__(self, exon):
        self.exon = exon
        
    def genomic2feature(self, location):
        print "exon start:end {}:{}".format(self.exon.loc_start, self.exon.loc_end)
        if location.strand == 1:
            start = location.start - self.exon.loc_start
            end = location.end - self.exon.loc_start
        else:
            start = self.exon.loc_start - location.start
            end = self.exon.loc_end - location.end
            
        return FeatureLocation(start, end, strand=location.strand, ref='feature')

class MissingTranscriptError(LookupError):
    "No transcript has been associated with the mapper"

class FeatureNotFound(LookupError):
    "The feature was not found"
    
class FeatureNonCoding(LookupError):
    "This transcript doesn't have a translation"

class UnknownCoordinateSystem(ValueError):
    "The coordinate system isn't known"
    
class WrongStrandError(TypeError):
    "The strand doesn't match between the given features"

