from django.test import TestCase
from ..lib.mapper import Mapper, MissingTranscriptError
from ..models import Exon, Transcript
from django.core.management import call_command
import json
import pprint

from . import test_fixtures_set

class MapperTestCase(TestCase):
    test_fixtures = test_fixtures_set
    
    def setUp(self):
        for fixture in self.test_fixtures:
            full_fixture = 'tark/tests/test-data/' + fixture
            print "Loading fixture {}".format(full_fixture)
            call_command('loaddata', full_fixture, verbosity=1)

    def testExon(self):
        print "Running Exon model test"

        exon = Exon.objects.get(pk=583852)
       
        print "Fetched exon {} at location {}".format(exon, exon.location)
        
        gene = exon.gene
        print "Exon {} belongs to gene {}".format(exon, gene)
        
        print "Exon feature has a sequence associated: {}".format(exon.has_sequence)
        print "Sequence: {}".format(exon.sequence.seq)
        print "Sub-sequence: {}".format(exon[43758945:43758949])
        print "Sub-sequence: {}".format(exon[:43758949])
        print "Sub-sequence: {}".format(exon[43759016:])
        print "Sub-sequence: {}".format(exon[43758945:43758949:2])
        print "Sub-sequence: {}".format(exon[43758947])
        print "Type of sub-sequence: {}".format (type(exon[43758945:43758949]) )

    def testTranscript(self):
        print "Running Transcript model test"
        
        transcript = Transcript.objects.get(pk=105417)
        
        print "Fetched transcript {} at location {}".format(transcript, transcript.location)
        print "Sub-sequence: {}".format(transcript[43759016:43768103])
        print "Sub-sequence: {}".format(transcript[43768101:43768103])
        print "Sub-sequence: {}".format(transcript[43768102:43768103])
        print "Sub-sequence: {}".format(transcript[43768103:43768103])
        print "Sub-sequence: {}".format(transcript[43768103])
        print "Sub-sequence: {}".format(transcript[43789542:])
        
        