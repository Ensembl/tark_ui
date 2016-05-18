from django.test import TestCase
from tark.lib.mapper import Mapper, MissingTranscriptError
from tark.models import Exon
from django.core.management import call_command
import json
import pprint

class MapperTestCase(TestCase):
    test_fixtures = ['session.json', 'genome.json', 'assembly.json', 'sequences.json', 'features.json', 'gene_names.json']
    
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