from django.test import TestCase, RequestFactory, Client
from django.http import HttpResponse, StreamingHttpResponse
from django.core.management import call_command
import json
import os
import pprint

'''

Test the /genome/ endpoint, fetching all available genomes

'''

REF_FILE_PATH = os.path.dirname(__file__)

class GenomesTestCase(TestCase):
    test_fixtures = ['session.json', 'genome.json', 'assembly.json', 'sequences.json', 'features.json', 'gene_names.json']

    def setUp(self):
        for fixture in self.test_fixtures:
            full_fixture = 'tark/tests/test-data/' + fixture
            print "Loading fixture {}".format(full_fixture)
            call_command('loaddata', full_fixture, verbosity=1)

        self.factory = RequestFactory()

    def fetchContent(self, response):
        if type(response) == StreamingHttpResponse:
            return response.streaming_content
        else:
            return response.content

    def fetchRecord(self, structure, key, value):
        for rec in structure:
            if rec.get(key, None) == value:
                return [rec]

        return []

    def testLookupGeneJson(self):
        print "Testing Lookup Gene (JSON)"
        c = Client()
        
        base_url = '/tark/lookup/gene/?id={}'

        json_file = open(os.path.join(REF_FILE_PATH, 'gene_json_reference.json'))
        ref_gene = json.load(json_file)
        genes = ['ENSG00000198001', 'ENSG00000134070']
        assemblies = {'GCA_000001405.20':'GRCh38.p5:GCA_000001405.20', 'GCA_000001405.14':'GRCh37.p13:GCA_000001405.14'}

        for gene in genes:
            url = base_url.format(gene)

            response = c.get(url)
            content = json.loads(''.join(self.fetchContent(response)))
            self.assertEquals(ref_gene[gene]['ref_gene'], content, "Retrieved genese for {} don't match".format(gene))

            response = c.get(url + '&expand=1')
            content = json.loads(''.join(self.fetchContent(response)))
            self.assertEquals(ref_gene[gene]['ref_gene_expanded'], content, "Retrieved genese for {} don't match".format(gene))

            for assembly in assemblies:
                response = c.get(url + '&assembly={}'.format(assembly))
                record = self.fetchRecord(ref_gene[gene]['ref_gene'], u'assembly', assemblies[assembly]);
                content = json.loads(''.join(self.fetchContent(response)))
                self.assertEquals(record, content, "Retrieved genese for {}, assembly {} don't match".format(gene, assembly))
            
    def testLookupTranscriptJson(self):
        print "Testing Lookup Transcript (JSON)"
        c = Client()
        
        base_url = '/tark/lookup/transcript/?id={}'

        json_file = open(os.path.join(REF_FILE_PATH, 'transcript_json_reference.json'))
        ref_transcript = json.load(json_file)
        transcripts = ['ENST00000448290', 'ENST00000256458']

        for transcript in transcripts:
            url = base_url.format(transcript)

            response = c.get(url)
            content = json.loads(''.join(self.fetchContent(response)))
            self.assertEquals(ref_transcript[transcript]['ref_transcript'], content, "Retrieved transcript for {} don't match".format(transcript))

            response = c.get(url + '&expand=1')
            content = json.loads(''.join(self.fetchContent(response)))
            self.assertEquals(ref_transcript[transcript]['ref_transcript_expanded'], content, "Retrieved transcript for {} (expanded) don't match".format(transcript))

            response = c.get(url + '&skip_sequence=1')
            content = json.loads(''.join(self.fetchContent(response)))
            self.assertEquals(ref_transcript[transcript]['ref_transcript_skip_seq'], content, "Retrieved transcript for {} (skip_sequence) don't match".format(transcript))


    def testLookupTranscriptFasta(self):
        print "Testing Lookup Transcript (Fasta)"
        json_data = {}
        c = Client()
        
        base_url = '/tark/lookup/transcript/?id={}'

        json_file = open(os.path.join(REF_FILE_PATH, 'transcript_fasta_reference.json'))
        ref_transcript = json.load(json_file)
        transcripts = ['ENST00000448290', 'ENST00000256458']

        for transcript in transcripts:
            json_data[transcript] = {}
            url = base_url.format(transcript)

            response = c.get(url, CONTENT_TYPE='text/x-fasta')
            content = self.fetchContent(response)
            json_data[transcript] = [x for x in content]
            self.assertEquals(ref_transcript[transcript], json_data[transcript], "Retrieved transcript for {} (skip_sequence) don't match".format(transcript))

#        with open(os.path.join(REF_FILE_PATH, 'transcript_fasta_reference.json'), 'w') as outfile:
#            pprint.pprint(json_data)
#            json.dump(json_data, outfile, sort_keys=True, indent=4)

    def testLookupTranscriptPostFasta(self):
        print "Testing Lookup Transcript Post (Fasta)"
        json_data = {}
        c = Client()

        url = '/tark/lookup/transcript/'

        json_file = open(os.path.join(REF_FILE_PATH, 'transcript_post_json_reference.json'))
        ref_transcript = json.load(json_file)
        transcript_data = {'id': ['ENST00000448290', 'ENST00000256458']}

        response = c.post(url, json.dumps(transcript_data),
                          content_type='application/json')
        content = json.loads(''.join(self.fetchContent(response)))
        self.assertEquals(ref_transcript, content, "Retrieved transcripts via POST don't match")

#        with open(os.path.join(REF_FILE_PATH, 'transcript_post_json_reference.json'), 'w') as outfile:
#            json.dump(content, outfile, sort_keys=True, indent=4)
