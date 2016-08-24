from django.test import TestCase, RequestFactory, Client
from django.core.management import call_command
import json

'''

Test the /assembly/{species}/ endpoint, fetching all available assemblies
for a species

'''

assemblies = [{"version": 14, "accession": "GCA_000001405", "genome": "homo_sapiens (9606)", "name": "GRCh37.p13"},{"version": 20, "accession": "GCA_000001405", "genome": "homo_sapiens (9606)", "name": "GRCh38.p5"}]

class AssembliesTestCase(TestCase):
    test_fixtures = ['session.json', 'genome.json', 'assembly.json', 'sequences.json', 'features.json', 'gene_names.json']

    def setUp(self):
        for fixture in self.test_fixtures:
            full_fixture = 'tark/tests/test-data/' + fixture
            print "Loading fixture {}".format(full_fixture)
            call_command('loaddata', full_fixture, verbosity=1)

        self.factory = RequestFactory()

    def testFetchAssemblies(self):
        c = Client()

        response = c.get('/tark/assembly/human/')
        content = json.loads(''.join(response.streaming_content))
        self.assertEquals(assemblies, content, "Retrieved available human assemblies don't match")
