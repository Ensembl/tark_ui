from django.test import TestCase, RequestFactory, Client
from django.http import HttpResponse, StreamingHttpResponse
from django.core.management import call_command
import json
import os
import pprint

'''

Test the /lookup/ endpoint

'''

REF_FILE_PATH = os.path.join(os.path.dirname(__file__), 'references')
WRITE_REFS = False

class GenomesTestCase(TestCase):
    test_fixtures = ['session.json', 'genome.json', 'assembly.json', 'sequences.json', 'features.json', 'gene_names.json']

    def setUp(self):
        for fixture in self.test_fixtures:
            full_fixture = 'tark/tests/test-data/' + fixture
            print "Loading fixture {}".format(full_fixture)
            call_command('loaddata', full_fixture, verbosity=1)

        self.factory = RequestFactory()

    def testLookupGene(self):
        print "Testing Lookup Gene"

        genes = ['ENSG00000198001', 'ENSG00000134070']
        assemblies = {'GCA_000001405.20':'GRCh38.p5:GCA_000001405.20', 'GCA_000001405.14':'GRCh37.p13:GCA_000001405.14'}

        self.doJsonLookupTest('gene', genes, assemblies, WRITE_REFS)

        self.doJsonPostTest('gene', genes, WRITE_REFS)

    def testLookupTranscript(self):
        print "Testing Lookup Transcript"

        transcripts = ['ENST00000448290', 'ENST00000256458', 'ENST12345']

        self.doJsonLookupTest('transcript', transcripts, None, WRITE_REFS)

        self.doFastaLookupTest('transcript', transcripts, WRITE_REFS)

        self.doJsonPostTest('transcript', transcripts, WRITE_REFS)

    def testLookupExon(self):
        print "Testing Lookup Exon"

        exons = ['ENSE00002345391', 'ENSE00003542427', 'ENSE12345']

        self.doJsonLookupTest('exon', exons, None, WRITE_REFS)

        self.doFastaLookupTest('exon', exons, WRITE_REFS)

        self.doJsonPostTest('exon', exons, WRITE_REFS)

    def testLookupTranslation(self):
        print "Testing Lookup Translation"

        translations = ['ENSP00000256458', 'ENSP00000390651']

        self.doJsonLookupTest('translation', translations, None, WRITE_REFS)

        self.doFastaLookupTest('translation', translations, WRITE_REFS)

        self.doJsonPostTest('translation', translations, WRITE_REFS)


    def doFastaLookupTest(self, feature_type, features, write_ref = False):
        print "\tdoFastaLookupTest ({})".format(feature_type)
        c = Client()
        
        base_url = "/tark/lookup/{}/?id={}"

        if write_ref:
            refs = {}
        else:
            json_file = open( os.path.join( REF_FILE_PATH, "{}_fasta_reference.json".format(feature_type) ) )
            refs = json.load(json_file)

        for feature in features:
            if write_ref:
                refs[feature] = {}

            url = base_url.format(feature_type, feature)

            response = c.get(url, CONTENT_TYPE='text/x-fasta')
            content = ''.join(self.fetchContent(response))
            if write_ref:
                refs[feature] = content
            else:
                self.assertEquals(refs[feature], content, "Retrieved {} for {} (fasta) don't match".format(feature_type, feature))

        if write_ref:
            outfile_name = os.path.join( REF_FILE_PATH, "{}_fasta_reference.json".format(feature_type) )
            print "\tWriting reference outfile {}".format(outfile_name)
            with open(outfile_name, 'w') as outfile:
                json.dump(refs, outfile, sort_keys=True, indent=4)


    def doJsonLookupTest(self, feature_type, features, assemblies, write_ref = False):
        print "\tdoJsonLookupTest ({})".format(feature_type)
        c = Client()
        
        base_url = "/tark/lookup/{}/?id={}"

        if write_ref:
            refs = {}
        else:
            json_file = open( os.path.join( REF_FILE_PATH, "{}_json_reference.json".format(feature_type) ) )
            refs = json.load(json_file)

        for feature in features:
            if write_ref:
                refs[feature] = {}

            url = base_url.format(feature_type, feature)

            response = c.get(url)
            content = json.loads(''.join(self.fetchContent(response)))
            if write_ref:
                refs[feature]['ref'] = content
            else:
                self.assertEquals( refs[feature]['ref'], content, "Retrieved {} for {} don't match".format(feature_type, feature) )

            response = c.get(url + '&expand=1')
            content = json.loads(''.join(self.fetchContent(response)))
            if write_ref:
                refs[feature]['ref_expanded'] = content
            else:
                self.assertEquals( refs[feature]['ref_expanded'], content, "Retrieved {} for {} (expanded) don't match".format(feature_type, feature) )

            response = c.get(url + '&skip_sequence=1')
            content = json.loads(''.join(self.fetchContent(response)))
            if write_ref:
                refs[feature]['ref_skip_sequence'] = content
            else:
                self.assertEquals( refs[feature]['ref_skip_sequence'], content, "Retrieved {} for {} (skip_sequence) don't match".format(feature_type, feature) )

            response = c.get(url + '&skip_sequence=1&expand=1')
            content = json.loads(''.join(self.fetchContent(response)))
            if write_ref:
                refs[feature]['ref_skip_sequence_expanded'] = content
            else:
                self.assertEquals( refs[feature]['ref_skip_sequence_expanded'], content, "Retrieved {} for {} (skip_sequence, expanded) don't match".format(feature_type, feature) )

            for assembly in assemblies or []:
                response = c.get(url + '&assembly={}'.format(assembly))
                record = self.fetchRecord(refs[feature]['ref'], u'assembly', assemblies[assembly]);
                content = json.loads(''.join(self.fetchContent(response)))
                if not write_ref:
                    self.assertEquals(record, content, "Retrieved {} for {}, assembly {} don't match".format(feature_type, feature, assembly))

        if write_ref:
            outfile_name = os.path.join( REF_FILE_PATH, "{}_json_reference.json".format(feature_type) )
            print "\tWriting reference outfile {}".format(outfile_name)
            with open(outfile_name, 'w') as outfile:
                json.dump(refs, outfile, sort_keys=True, indent=4)

    def doJsonPostTest(self, feature_type, features, write_ref = False):
        print "\tdoJsonPostTest ({})".format(feature_type)
        c = Client()

        url = "/tark/lookup/{}/".format(feature_type)

        if write_ref:
            refs = {}
        else:
            json_file = open( os.path.join( REF_FILE_PATH, "{}_post_json_reference.json".format(feature_type) ) )
            refs = json.load(json_file)

        features_data = {'id': features}

        response = c.post(url, json.dumps(features_data),
                          content_type='application/json')
        content = json.loads(''.join(self.fetchContent(response)))

        if write_ref:
            outfile_name = os.path.join( REF_FILE_PATH, "{}_post_json_reference.json".format(feature_type) )
            print "\tWriting reference outfile {}".format(outfile_name)
            with open(outfile_name, 'w') as outfile:
                json.dump(content, outfile, sort_keys=True, indent=4)
        
        else:
            self.assertEquals( refs, content, "Retrieved {} via POST don't match".format(feature_type) )


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
