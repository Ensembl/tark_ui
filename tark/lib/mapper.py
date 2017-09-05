from Bio.SeqFeature import FeatureLocation
from Bio.Seq import Seq

#from tark.exceptions import MissingTranscriptError, FeatureNonCoding, UnknownCoordinateSystem, FeatureNotFound,\
from ..exceptions import MissingTranscriptError, FeatureNonCoding, UnknownCoordinateSystem, FeatureNotFound,\
WrongStrandError, IncompatibleFeatureType, LocationNotInFeature
import pprint

class Mapper():
    def __init__(self, feature):
        if feature.feature_type == 'gene':
            gene = feature
        if feature.feature_type == 'transcript':
            gene = feature.gene
            self.transcript(feature)
            
        self._gene = gene
        print "gene start:end {}:{}".format(gene.loc_start, gene.loc_end)

    
    def transcript(self, transcript=None):
        if transcript:
            self.flush_mapper()
            self._transcript = transcript
            return self
        
        if not hasattr(self, '_transcript'):
            raise MissingTranscriptError()
        
        return self._transcript

    def is_exonic(self, location):
        if location.ref_db and location.ref_db != 'genomic':
            location = self.remap2genomic(location)

        return True if self.fetch_exon(location) else False

    def is_cds(self, location):
        if not self.is_exonic(location):
            return False
        
        if location in self.transcript().translation:
            return True
       
        return False
       
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

    def flush_mapper(self):
        self.flush_exons()
        self.flush_translations()
        
    @property
    def exons(self):
        if not hasattr(self, '_transcript'):
            raise MissingTranscriptError()
        
        if not hasattr(self, '_exons'):
            self._exons = list(self.transcript().exons)

        return self._exons

    def flush_exons(self):
        if hasattr(self, '_exons'):
            del self._exons
    
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

    def flush_translations(self):
        if hasattr(self, '_translations'):
            del self._translations

    def remap2genomic(self, location):
        # If no type is given we assume it's already genomic
        if not location.ref_db or location.ref_db == 'genomic':
            return FeatureLocation(location.start, location.end, strand=location.strand, ref_db='genomic')

        if location.ref_db == 'gene':
            feature = self._gene
        elif location.ref_db == 'transcript':
            feature = self.transcript()
        elif location.ref_db == 'exon' and location.ref:
            feature = self.find_exon(location.ref)
            if not feature:
                raise FeatureNotFound()

        elif (location.ref_db == 'cdna' or location.ref_db == 'cds') and location.ref:
            feature = self.find_translation(location.ref)
            if not feature:
                raise FeatureNonCoding()

            if location.ref_db == 'cds':
                location = FeatureLocation( location.start * 3, 
                                            location.end * 3, 
                                            strand=location.strand,
                                            ref=location.ref,
                                            ref_db=location.ref_db )
        else:
            raise UnknownCoordinateSystem()
            
        print feature
        print feature.location
        if location.strand == 1:
            start = location.start + feature.loc_start
            end = location.end + feature.loc_start
        else:
            start = feature.loc_end - location.start
            end = feature.loc_end - location.end
            
        return FeatureLocation(start, end, strand=location.strand, ref=feature.stable_id, ref_db='genomic')

    @property
    def cdna_coding_start(self):
        translation = self.transcript().translation
        cdna_offset = 0
        
        for exon in self.exons:
            if translation.loc_start in exon:
                cdna_offset = translation.loc_start - exon.loc_start + cdna_offset + 1
                return cdna_offset
            
            cdna_offset = cdna_offset + exon.length
            
        return None

    @property
    def cdna_coding_end(self):
        translation = self.transcript().translation
        cdna_offset = 0
        
        for exon in self.exons:
            if translation.loc_end in exon:
                cdna_offset =  translation.loc_end - exon.loc_start + cdna_offset + 1
                return cdna_offset
            
            cdna_offset = cdna_offset + exon.length
            
        return None

    def genomic2cdna(self, location):
        feature = self.transcript()
        cdna_offset = 0
        prior_exon_end = 0
        coordinates = []

        if location.strand != feature.loc_strand:
            raise WrongStrandError()
        
        if location.ref_db != 'genomic':
            location = self.remap2genomic(location)

        # It's all gap
        if location not in feature:
            coordinates.append(FeatureLocation(location.start, location.end,
                                               strand=location.strand,
                                               ref_db='gap'))
            return coordinates
        
        for exon in self.exons:
            # We haven't hit the start exon yet, move on
            if not coordinates:
                if location.start > exon.loc_end:
                    cdna_offset = cdna_offset + exon.length
                    continue
                elif location.start < exon.loc_start:
                    # We started in a gap
                    coordinates.append(FeatureLocation(location.start,
                                                       exon.loc_start - 1,
                                                       strand=location.strand,
                                                       ref_db='gap'))                    

            if location.end < exon.loc_start:
                # If the region ends in the gap, we need to push
                # one more gap on the coordinates
                if prior_exon_end and location.end > prior_exon_end:
                    coordinates.append(FeatureLocation(prior_exon_end + 1,
                                                       location.end,
                                                       strand=location.strand,
                                                       ref_db='gap'))

                # And we're done, end
                break

            # We're in the desired location, add any gaps we're
            # encountered until now
            if prior_exon_end:
                coordinates.append(FeatureLocation(prior_exon_end + 1,
                                                   exon.loc_start - 1,
                                                   strand=location.strand,
                                                   ref_db='gap') )
            prior_exon_end = exon.loc_end
            
            if location.start in exon:
                feature_start = location.start - exon.loc_start + cdna_offset + 1
            else:
                feature_start = cdna_offset + 1
                
            if location.end in exon:
                feature_end = location.end - exon.loc_start + cdna_offset + 1
            else:
                feature_end = cdna_offset + exon.length
                
            coordinates.append(FeatureLocation(feature_start, feature_end,
                                               strand=location.strand,
                                               ref_db='cdna', ref=exon.stable_id))
            
            cdna_offset = cdna_offset + exon.length
            
        return coordinates
            
 
    def genomic2cds(self, location):
        if location.start > location.end + 1:
            raise WrongStrandError()
        
        out = []
        try:
            translation = self.transcript().translation
        except FeatureNonCoding:
            # Pseudogene, return a gap
            out.append(FeatureLocation(location.start, location.end,
                                       strand=location.strand,
                                       ref_db='gap'))
            return out
            
        cdna_cstart = self.cdna_coding_start
        cdna_cend = self.cdna_coding_end
        
        coords = self.genomic2cdna(location)
        print coords
        
        for coord in coords:
            if coord.ref_db == 'gap':
                out.append(coord)
                
            else:
                start = coord.start
                end = coord.end
                
                if (coord.strand == -1) or (end < cdna_cstart) or (start > cdna_cend):
                    # is all gap - does not map to peptide
                    # WRONG, copied from perl api, start/end are in cdna coords not genomic!
                    out.append(FeatureLocation(start, end,
                                               strand=location.strand,
                                               ref_db='gap'))
                else:
                    # we know area is at least partially overlapping CDS
                    cds_start = start - cdna_cstart + 1
                    cds_end = end - cdna_cstart + 1
                    print "start: {}, cdna_cstart: {}".format(start, cdna_cstart)
                    print "cds_start: {} cds_end: {}".format(cds_start, cds_end)
                    
                    if start < cdna_cstart:
                        # start of coordinates are in the 4prime UTR
                        out.append(FeatureLocation(start, cdna_cstart-1,
                                                   strand=location.strand,
                                                   ref_db='gap'))
                        
                        # start is now relative to start of CDS
                        cds_start = 1
                        
                    end_gap = None
                    if end > cdna_cend:
                        # end of cooordinates are in the 3prime UTR
                        end_gap = FeatureLocation(cdna_cend+1, end,
                                                  strand=location.strand,
                                                  ref_db='gap')
                        # adjust end to realtive to CDS start
                        cds_end = cdna_cend - cdna_cstart + 1
                        
                    # start and end are now entirely in CDS and relative to CDS start
                    print "cds_start: {} cds_end: {}".format(cds_start, cds_end)
                    adjusted_coord = FeatureLocation(cds_start, cds_end,
                                                     strand=location.strand,
                                                     ref_db='cds')
                    out.append(adjusted_coord)
                    
                    if end_gap:
                        out.append(end_gap)
                        
        return out

    def genomic2pep(self, location):
        if location.start > location.end + 1:
            raise WrongStrandError()

        coords = self.genomic2cds(location)
        out = []
        
        for coord in coords:
            if coord.ref_db == 'gap':
                out.append(coord)
                continue
              
            shift = coord.start % 3
            adjusted_coord = FeatureLocation(int((coord.start + shift + 2) / 3.0),
                                             int((coord.end   + shift) / 3.0),
                                             strand=location.strand,
                                             ref_db='pep')
            out.append(adjusted_coord)
            
        return out

    def transcript_subseq(self, start, end, introns=False):
        transcript = self.transcript()
        location = transcript.genomic_location(start, end)
        if location not in transcript:
            return Seq('')

        final_seq = Seq('')
        seq = transcript.seq
        coords = self.genomic2cdna(location)
        
        for coord in coords:
            if coord.ref_db == 'gap' and introns:
                final_seq += Seq('N' * (len(coord)+1))
            elif coord.ref_db == 'cdna':
                final_seq += seq[coord.start-1:coord.end]

        return final_seq
        
    def remap2feature(self, location, feature, coordinates='cdna'):
        if location.strand != feature.loc_strand:
            raise WrongStrandError()
        
        if location.ref_db != 'genomic':
            location = self.remap2genomic(location)

        if coordinates == 'cds':
            if feature.feature_type == 'transcript':
                feature = feature.translation
            elif feature.feature_type == 'translation':
                pass
            else:
                raise IncompatibleFeatureType()
            
            feature_location = self._remap2feature(location, feature)
            
            return FeatureLocation( int( feature_location.start / 3.0 ), 
                                    int( feature_location.end / 3.0 ), 
                                    strand=feature_location.strand,
                                    ref=feature_location.ref,
                                    ref_db='cds' )

        else:
            return self._remap2feature(location, feature)

    def _remap2feature(self, location, feature):

        if not location in feature:
            print "TYPE: {}".format(feature.feature_type)
            print "FEATURE: {}".format(feature.location)
            print "LOCATION: {}".format(location)
            raise LocationNotInFeature()
          
        print "feature start:end {}:{}".format(feature.loc_start, feature.loc_end)
        if location.strand == 1:
            start = location.start - feature.loc_start
            end = location.end - feature.loc_start
        else:
            start = feature.loc_start - location.start
            end = feature.loc_end - location.end

        ref_db = feature.feature_type if feature.feature_type != 'translation' else 'cdna'
        return FeatureLocation(start, end, strand=location.strand, ref=feature.stable_id, ref_db=ref_db)
            

    def __str__(self):
        return "Mapper for gene {}, location: {}".format(self._gene, 
                                                         self._gene.location)

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
        for exon_transcript in self.transcript.exons.order_by('exon_order'):
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

