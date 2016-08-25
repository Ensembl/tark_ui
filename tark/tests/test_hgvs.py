from django.test import TestCase, RequestFactory, Client
from django.core.management import call_command
from django.http import HttpResponse, StreamingHttpResponse
import json
import os
import pprint

'''

Test the /hgvs/ endpoint

'''

REF_FILE_PATH = os.path.join(os.path.dirname(__file__), 'references')
WRITE_REFS = False

class HGVSTestCase(TestCase):
    test_fixtures = ['session.json', 'genome.json', 'assembly.json', 'sequences.json', 'features.json', 'gene_names.json']
    
    def setUp(self):
        for fixture in self.test_fixtures:
            full_fixture = 'tark/tests/test-data/' + fixture
            print "Loading fixture {}".format(full_fixture)
            call_command('loaddata', full_fixture, verbosity=1)

        self.factory = RequestFactory()

    def testHGVSGene(self):
        print "Running HGVS g test"

        c = Client()

        self.HGVSLookup('human', 'IRAK4:g.2C>T', WRITE_REFS)

    def HGVSLookup(self, species, hgvs_string, write_ref = False):
        print "\tHGVSLookup ({})".format(hgvs_string)
        hgvs_test_name = self.cleanHGVS(hgvs_string)

        c = Client()

        url = '/tark/lookup/hgvs/{}/{}/'.format(species, hgvs_string)

        if write_ref:
            refs = {}
        else:
            json_file = open( os.path.join( REF_FILE_PATH, "hgvs_{}_{}_reference.json".format(species, hgvs_test_name) ) )
            refs = json.load(json_file)

        response = c.get(url)
        content = json.loads(''.join( self.fetchContent(response) ) )

        if write_ref:
            outfile_name = os.path.join( REF_FILE_PATH, "hgvs_{}_{}_reference.json".format(species, hgvs_test_name) )
            print "\tWriting reference outfile {}".format(outfile_name)
            with open(outfile_name, 'w') as outfile:
                json.dump(content, outfile, sort_keys=True, indent=4)
        else:
            self.assertEqual(refs, content, "HGVS {} {} doesn't match content returned from the server".format(species, hgvs_string))


    def cleanHGVS(self, hgvs_string):

        hgvs_string = hgvs_string.replace('>', '_');
        hgvs_string = hgvs_string.replace(':', '_');
        
        return hgvs_string

    def fetchContent(self, response):
        if type(response) == StreamingHttpResponse:
            return response.streaming_content
        else:
            return response.content
