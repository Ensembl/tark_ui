from django.test import TestCase
from tark.models import Feature, Assembly, Sequence
from django.core.management import call_command

# Create your tests here.
class LookupTestCase(TestCase):
#    test_fixtures = ['session.json', 'genome.json', 'assembly.json', 'sequences.json']
    test_fixtures = ['session.json', 'genome.json', 'assembly.json', 'sequences.json', 'features.json', 'gene_names.json']
    
    def setUp(self):
        for fixture in self.test_fixtures:
            full_fixture = 'tark/tests/test-data/' + fixture
            print "Loading fixture {}".format(full_fixture)
            call_command('loaddata', full_fixture, verbosity=1)
    
    def testFeatureLookup(self):
        print "Running test"
        for a in Assembly.objects.all():
            print str(a)
            
        for s in Sequence.objects.all():
            print s
            
    