import hgvs.parser
from hgvs.edit import NARefAlt
from Bio.Seq import Seq
from tark.models import Transcript, Assembly, VariationFeature
from tark.lib.mapper import TranscriptMapper, ExonMapper, Mapper
from Bio.SeqFeature import FeatureLocation
from tark.lib.iterables import iterlist
import pprint

def fetch_by_hgvs(hgvs_str, species, **kwargs):
    hp = hgvs.parser.Parser()
    variant = hp.parse_hgvs_variant(hgvs_str)
    pprint.pprint(variant)
    
    assemblies = Assembly.fetch_by_name(species)    
    # Check is assemblies came back
    
    if variant.type == 'g':
        matched_features = _fetch_by_hgvs_g(variant, assemblies, **kwargs)

    return matched_features


def _fetch_by_hgvs_g(variant, assemblies, **kwargs):
    print "name: {}".format(variant.ac)
    results = []
    for assembly in assemblies.all():
        
        transcripts = Transcript.objects.fetch_by_name(variant.ac, assembly=assembly, **kwargs)
        print transcripts
        matched_features = []
    
        for transcript in transcripts: 
            mapper = TranscriptMapper(transcript)
        
            location = FeatureLocation(variant.posedit.pos.start.base, variant.posedit.pos.end.base, transcript.loc_strand, ref='feature')
#            if not mapper.is_exonic(variant.posedit.pos.start.base, variant.posedit.pos.end.base, 'feature'):
            if not mapper.is_exonic(location):
#                matched_features.append(transcript)
#                matched_features.append(VariationFeature(feature=transcript))
                continue
            
            # It is exonic, so we need to convert the coordinates of the hgvs to feature based
            #location = FeatureLocation(variant.posedit.pos.start.base, variant.posedit.pos.end.base, transcript.loc_strand)
            exon = mapper.fetch_exon(location)
            exon_mapper = ExonMapper(exon)
            exon_location = exon_mapper.genomic2feature(mapper.feature2genomic(location))
            alt_seq = _apply_hgvs_g(variant, exon, exon_location)
            if alt_seq:
                matched_features.append({'transcript': transcript, 'exon': exon, 'alt_seq': alt_seq})
            
        if matched_features:
            print "type: {}".format(type(matched_features))
            pprint.pprint(matched_features)
            results.append({'assembly': assembly, 'transcripts': matched_features}) 
        
#    pprint.pprint(results)
    return results

def _apply_hgvs_g(variant, feature, location):

    print "Applying"
    print feature
    print location
    edit = variant.posedit.edit
#    if not location in feature:
#        raise Exception("Location {} isn't in feature {}".format(location, feature))
    
    seq = feature.seq_checksum.sequence
    print type(seq)
    print len(seq)
    
    # First figure out the type of edit to apply
    if type(edit) == NARefAlt:
        print seq
        print location
        subseq = seq[location.start:location.end+1]
        if edit.ref and edit.ref != subseq:
            return None
#            raise Exception("Reference {} does not match hgvs edit {}".format(subseq, edit.ref))
       
        # It matches! (or there's no reference given in the hgvs)
        alt_seq = seq[:location.start-1] + Seq(edit.alt) + seq[location.end:]
        if len(alt_seq) != (len(seq) + variant.posedit.length_change()):
            raise Exception("New sequence not the length it should be, expecting {}, got {}".format( (len(seq) + edit.length_change()), len(alt_seq)))
        
        return alt_seq
    
    