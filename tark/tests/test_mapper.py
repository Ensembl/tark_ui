from django.test import TestCase
from ..lib.mapper import Mapper, MissingTranscriptError, LocationNotInFeature
from ..models import Gene, Transcript, Exon
from django.core.management import call_command
import json
import pprint

class MapperTestCase(TestCase):
    test_fixtures = ['session.json', 'genome.json', 'assembly.json', 'sequences.json', 'features.json', 'gene_names.json']
    
    def setUp(self):
        for fixture in self.test_fixtures:
            full_fixture = 'tark/tests/test-data/' + fixture
#            print "Loading fixture {}".format(full_fixture)
            call_command('loaddata', full_fixture, verbosity=1)
            
    def testMapperGene(self):
        print "Running mapper(gene) test"

        gene = Gene.objects.get(pk=36751)
        mapper = Mapper(gene)
        print mapper
        
#        transcript = mapper.transcript()
        with self.assertRaises(LookupError, msg="MissingTranscriptError wasn't thrown when accessing transcript"):
            mapper.transcript()
        
    def testMapperTranscript(self):
        print "Running mapper(transcript) test"

        transcript = Transcript.objects.get(pk=105417)
        mapper = Mapper(transcript)
        print mapper

        transcript = mapper.transcript()
        print "Fetched transcript {} from mapper: {}".format(transcript, transcript.location)
        
        
    def testMapperfromGene(self):
        print "Running fetching mapper from a gene test"

        gene = Gene.objects.get(pk=36751)
        mapper = gene.mapper
        print mapper

        location = gene.genomic_location(43758944, 43758947)
        feature_location = gene.remap2feature(location)
        print "Location remapped to feature {} is {}".format(gene, feature_location)
        self.assertEqual(feature_location.start, 0, "Remapped feature start is 0")
        self.assertEqual(feature_location.end, 3, "Remapped feature end is 3")

    def testMapperfromTranscript(self):
        print "Running fetching mapper from a transcript test"

        transcript = Transcript.objects.get(pk=105417)
        mapper = transcript.mapper
        print mapper
        
        transcript = mapper.transcript()
        print "Fetched transcript {} from mapper: {}".format(transcript, transcript.location)
        
        location = transcript.genomic_location(43758944, 43758947)
        print "Is {} exonic: {}".format( location, mapper.is_exonic(location) )
        
        print "Is {} coding: {}".format( transcript, transcript.is_coding )

        feature_location = transcript.remap2feature(location)
        print "Location remapped to feature {} is {}".format(transcript, feature_location)
        self.assertEqual(feature_location.start, 0, "Remapped feature start is 0")
        self.assertEqual(feature_location.end, 3, "Remapped feature end is 3")
        
        feature_location2 = transcript.feature_location(0, 3)
        print "Created feature_location for {} at {}".format(transcript, feature_location2)
        self.assertEqual(str(feature_location), str(feature_location2), "feature_locations are not equal, 0-3")

        genomic_location = transcript.remap2genomic(feature_location)
        print "Remapped to genomic for {} give {}".format(transcript, genomic_location)
        self.assertEqual(str(genomic_location), str(location), "genomic locations don't match")
        
        with self.assertRaises(LocationNotInFeature, msg="Location {} isn't in the coding region of {}".format(location, transcript)):
            cds_location = transcript.remap2feature(location, coordinates='cds')
            
        genomic_cds_location = transcript.genomic_location(43772280, 43772300)
        print "Created genomic location inside coding region of {} at {}".format(transcript, genomic_cds_location)
        cds_location = transcript.remap2feature(genomic_cds_location, coordinates='cds')
        print "Mapped to cds location of {} at {}".format(transcript, cds_location)
        cdna_location = transcript.genomic2cdna(genomic_cds_location)
        print "Mapped to cdna location of {} at {}".format(transcript, cdna_location)
        transcript_location = transcript.remap2feature(genomic_cds_location)
        print "Mapped to transcript feature location of {} at {}".format(transcript, transcript_location)
        
        genomic_cdna_location = transcript.genomic_location(43772280, 43772300)
        print "\nTesting genomic2cdna of {} at {}".format(transcript, genomic_cdna_location)
        cdna_location = transcript.genomic2cdna(genomic_cdna_location)
        print "Mapped to transcript feature location of {} at {}\n".format(transcript, cdna_location)
        self.assertEqual(cdna_location[0].start, 344, "genomic2cdna within an exon start should be 344, got {}".format(cdna_location[0].start))
        self.assertEqual(cdna_location[0].end, 364, "genomic2cdna within an exon end should be 364, got {}".format(cdna_location[0].end))

        genomic_cdna_location = transcript.genomic_location(43772280, 43772920)
        print "\nTesting genomic2cdna of {} at {}".format(transcript, genomic_cdna_location)
        cdna_location = transcript.genomic2cdna(genomic_cdna_location)
        print "Mapped to transcript feature location of {} at {}\n".format(transcript, cdna_location)
        self.assertEqual(cdna_location[0].start, 344, "genomic2cdna across introns start should be 344, got {}".format(cdna_location[0].start))
        self.assertEqual(cdna_location[0].end, 426, "genomic2cdna across introns end should be 426, got {}".format(cdna_location[0].end))
        self.assertEqual(cdna_location[1].start, 43772363, "genomic2cdna intron gap should be 43772363, got {}".format(cdna_location[1].start))
        self.assertEqual(cdna_location[2].end, 435, "genomic2cdna second exon end should be 435, got {}".format(cdna_location[2].end))

        genomic_cds_location = transcript.genomic_location(43772280, 43772300)
        print "\nTesting genomic2cds of {} at {}".format(transcript, genomic_cds_location)
        cds_location = transcript.genomic2cds(genomic_cds_location)
        print "Mapped to transcript feature location of {} at {}\n".format(transcript, cds_location)
        self.assertEqual(cds_location[0].start, 36, "genomic2cds within an exon start should be 36, got {}".format(cds_location[0].start))
        self.assertEqual(cds_location[0].end, 56, "genomic2cds within an exon end should be 56, got {}".format(cds_location[0].end))

        genomic_cds_location = transcript.genomic_location(43772280, 43772920)
        print "\nTesting genomic2cds of {} at {}".format(transcript, genomic_cds_location)
        cds_location = transcript.genomic2cds(genomic_cds_location)
        print "Mapped to transcript feature location of {} at {}\n".format(transcript, cds_location)
        self.assertEqual(cds_location[0].start, 36, "genomic2cds across introns start should be 36, got {}".format(cds_location[0].start))
        self.assertEqual(cds_location[0].end, 118, "genomic2cdna across introns end should be 118, got {}".format(cds_location[0].end))
        self.assertEqual(cds_location[1].start, 43772363, "genomic2cdna intron gap should be 43772363, got {}".format(cds_location[1].start))
        self.assertEqual(cds_location[2].end, 127, "genomic2cdna second exon end should be 127, got {}".format(cds_location[2].end))

        genomic_pep_location = transcript.genomic_location(43772280, 43772921)
        print "\nTesting genomic2pep of {} at {}".format(transcript, genomic_pep_location)
        pep_location = transcript.genomic2pep(genomic_pep_location)
        print "Mapped to transcript feature location of {} at {}\n".format(transcript, pep_location)
        self.assertEqual(pep_location[0].start, 12, "genomic2pep across introns start should be 12, got {}".format(pep_location[0].start))
        self.assertEqual(pep_location[0].end, 39, "genomic2cdna across introns end should be 39, got {}".format(pep_location[0].end))
        self.assertEqual(pep_location[1].start, 43772363, "genomic2cdna intron gap should be 43772363, got {}".format(pep_location[1].start))
        self.assertEqual(pep_location[2].end, 43, "genomic2cdna second exon end should be 43, got {}".format(pep_location[2].end))

    def testMapperfromExon(self):
        print "Running fetching mapper from an exon test"

        exon = Exon.objects.get(pk=583852)

        mapper = exon.mapper
        print mapper
