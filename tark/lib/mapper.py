from Bio.SeqFeature import FeatureLocation
from tark.models import Feature, Transcript
import pprint

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
