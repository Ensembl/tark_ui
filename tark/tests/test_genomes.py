from django.test import TestCase, RequestFactory, Client
from django.core.management import call_command
import json

'''

Test the /genome/ endpoint, fetching all available genomes

'''

genomes = [{'taxonomy': 9606, 'name': 'homo_sapiens'}]

class GenomesTestCase(TestCase):
    test_fixtures = ['session.json', 'genome.json', 'assembly.json', 'sequences.json', 'features.json', 'gene_names.json']

    def setUp(self):
        for fixture in self.test_fixtures:
            full_fixture = 'tark/tests/test-data/' + fixture
            print "Loading fixture {}".format(full_fixture)
            call_command('loaddata', full_fixture, verbosity=1)

        self.factory = RequestFactory()

    def testFetchGenomes(self):
        c = Client()

        response = c.get('/tark/genome/')
        content = json.loads(''.join(response.streaming_content))
        self.assertEquals(genomes, content, "Retrieved available genomes don't match")
