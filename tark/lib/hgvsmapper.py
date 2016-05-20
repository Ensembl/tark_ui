import hgvs.parser
from hgvs.edit import NARefAlt, Dup, Inv
from Bio.Seq import Seq
from ..models import Gene, Transcript, Assembly, VariationFeature
from ..lib.mapper import TranscriptMapper, ExonMapper
from Bio.SeqFeature import FeatureLocation
import pprint

def fetch_by_hgvs(hgvs_str, species, **kwargs):
    hp = hgvs.parser.Parser()
    variant = hp.parse_hgvs_variant(hgvs_str)
    pprint.pprint(variant)
    
    assemblies = Assembly.fetch_by_name(species)    
    # Check if assemblies came back
    
    if variant.type == 'g':
        matched_features = _fetch_by_hgvs_g(variant, assemblies, **kwargs)

    return matched_features

def _fetch_by_hgvs_g(variant, assemblies, **kwargs):
    results = []
    
    for assembly in assemblies.all():
        location = None
        matched_features = []

        if variant.ac.startswith('chr'):
            # We have to do this hack for searching by genomic coordinates because
            # the hgvs library doesn't support numbers as the "location" in an hgvs string
            region = variant.ac.strip('chr')
            location = FeatureLocation(variant.posedit.pos.start.base, variant.posedit.pos.end.base, strand=1, ref='genomic')
            genes = Gene.objects.fetch_by_location(region, location, assembly=assembly, **kwargs)
        else:
            genes = Gene.objects.fetch_by_name(variant.ac, assembly=assembly, **kwargs)
        if not genes:
            print "we're in real trouble, we found nothing"
#            return results
            
        for gene in genes:
            matched_transcripts = []

#         transcripts = Transcript.objects.fetch_by_name(variant.ac, assembly=assembly, **kwargs)
#         if not transcripts:
#             # NOT IMPLEMENTED YET
#             transcripts = Transcript.objects.fetch_by_location(variant.loc, assembly=assembly, **kwargs)
#             location = FeatureLocation(variant.posedit.pos.start.base, variant.posedit.pos.end.base, ref='genomic')
#             if not transcripts:
#                 print "we're in real trouble, we found nothing"
#                 return results

        # We're going to do everything in genomic coordinates            
            for transcript in gene.transcripts.all():
                gene.mapper.transcript(transcript)
                if not location:
                    location = gene.feature_location(variant.posedit.pos.start.base, variant.posedit.pos.end.base)
                    location = transcript.mapper.remap2genomic(location)
                    print location
                # If there's no location it means we've been given a feature based location, remap!
                    
                    
                # We now have a genomic based location, let's apply to edit
    
                # If the location isn't exonic, we can't do anything to it
                # We only support locations that don't span exons for g right now
                # since we don't have the genomic sequence available, this may change
                if not transcript.mapper.is_exonic(location):
                    print "{} at {} is not exonic".format(transcript, location)
                    matched_transcripts.append({'transcript': transcript, 'info': 'Location does not fall in an exon'})
                    continue
    
                alt_transcript_seq = _apply_hgvs_g(variant, transcript, location)
                
                if alt_transcript_seq:
                    matched_transcripts.append({'transcript': transcript, 'alt_seq': alt_transcript_seq})
                else:
                    matched_transcripts.append({'transcript': transcript, 'info': "Transcript was not changed"})
                    
            if matched_transcripts:
                matched_features.append({'gene': gene, 'transcripts': matched_transcripts})
            
#        if matched_features:
        results.append({'assembly': assembly, 'genes': matched_features}) 
            
    return results
            

def _fetch_by_hgvs_g_old(variant, assemblies, **kwargs):
    print "name: {}".format(variant.ac)
    results = []
    for assembly in assemblies.all():
        
        transcripts = Transcript.objects.fetch_by_name(variant.ac, assembly=assembly, **kwargs)
        print transcripts
        matched_features = []
    
        for transcript in transcripts: 
            mapper = TranscriptMapper(transcript)
        
            location = FeatureLocation(variant.posedit.pos.start.base, variant.posedit.pos.end.base, strand=transcript.loc_strand, ref='feature')
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
    print "Applying!"
    
    edit = variant.posedit.edit

    if type(edit) == NARefAlt:
        subseq = feature[location.start:location.end]
        print "Found subseq: {}, type: {}".format(subseq, type(subseq))
        if edit.ref and edit.ref != subseq:
            return None
        
        # To handle deletions where there is no alternative sequence
        alt = Seq(edit.alt) if edit.alt else Seq('')
        
        # To handle insertions where we're keeping the existing sequences
        # and sliding the alt between the two base pair (and there should
        # only be two)
        if edit.type == 'ins':
            if len(subseq) != 2:
                raise Exception("Length of subseq for insertions must be exactly 2")
            alt = subseq[0:1] + alt + subseq[1:2]
            
        alt_seq = feature[:location.start-1] + alt + feature[location.end+1:]
        if len(alt_seq) != (len(feature.sequence.seq) + variant.posedit.length_change()):
            raise Exception("New sequence not the length it should be, expecting {}, got {}".format( (len(feature.sequence.seq) + variant.posedit.length_change()), len(alt_seq)))
    elif type(edit) == Dup:
        subseq = feature[location.start:location.end]
        print "Found subseq: {}, type: {}".format(subseq, type(subseq))
        if edit.ref and edit.ref != subseq:
            return None

        alt_seq = feature[:location.start-1] + subseq + subseq + feature[location.end+1:]
        if len(alt_seq) != (len(feature.sequence.seq) + variant.posedit.length_change()):
            raise Exception("New sequence not the length it should be, expecting {}, got {}".format( (len(feature.sequence.seq) + variant.posedit.length_change()), len(alt_seq)))
        
    elif type(edit) == Inv:
        subseq = feature[location.start:location.end]
        print "Found subseq: {}, type: {}".format(subseq, type(subseq))
        if edit.ref and edit.ref != subseq:
            return None

        alt = subseq.complement()

        alt_seq = feature[:location.start-1] + alt + feature[location.end+1:]
        if len(alt_seq) != (len(feature.sequence.seq) + variant.posedit.length_change()):
            raise Exception("New sequence not the length it should be, expecting {}, got {}".format( (len(feature.sequence.seq) + variant.posedit.length_change()), len(alt_seq)))
        
    else:
        print type(edit)
        print type(variant.posedit)
        print dir(edit)
        pprint.pprint(edit)
        
    return alt_seq

def _apply_hgvs_g_old(variant, feature, location):

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
    
    