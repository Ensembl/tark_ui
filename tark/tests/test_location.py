from django.test import TestCase, RequestFactory, Client
from django.http import HttpResponse, StreamingHttpResponse
from django.core.management import call_command
import json
import os
import pprint

from . import test_fixtures_set

'''

Test the /location/ endpoint

'''

REF_FILE_PATH = os.path.join(os.path.dirname(__file__), 'references')
WRITE_REFS = False

locations = [{"start": 3758945, "end": 43758945, "region": "12"},
             {"start": 3757945, "end": 43759945, "region": "12"},
             {"start": 3758945, "end": 1, "region": "12"},
             {"start": 3758945, "end": 43758945, "region": "Z"} ]

class LocationTestCase(TestCase):
    test_fixtures = test_fixtures_set

    def setUp(self):
        for fixture in self.test_fixtures:
            full_fixture = 'tark/tests/test-data/' + fixture
            print "Loading fixture {}".format(full_fixture)
            call_command('loaddata', full_fixture, verbosity=1)

        self.factory = RequestFactory()

    def testLookupGene(self):
        print "Testing Location Gene"

        assemblies = {'GCA_000001405.22':'GRCh38', 'GCA_000001405.14':'GRCh37'}

        self.doJsonLookupTest('gene', locations, assemblies, WRITE_REFS)

        self.doJsonPostTest('gene', locations, WRITE_REFS)

    def testLookupTranscript(self):
        print "Testing Location Transcript"

        self.doJsonLookupTest('transcript', locations, None, WRITE_REFS)

        self.doFastaLookupTest('transcript', locations, WRITE_REFS)

        self.doJsonPostTest('transcript', locations, WRITE_REFS)

    def testLookupExon(self):
        print "Testing Location Exon"

        self.doJsonLookupTest('exon', locations, None, WRITE_REFS)

        self.doFastaLookupTest('exon', locations, WRITE_REFS)

        self.doJsonPostTest('exon', locations, WRITE_REFS)

    def testLookupTranslation(self):
        print "Testing Location Translation"

        self.doJsonLookupTest('translation', locations, None, WRITE_REFS)

        self.doFastaLookupTest('translation', locations, WRITE_REFS)

        self.doJsonPostTest('translation', locations, WRITE_REFS)


    def doFastaLookupTest(self, feature_type, features, write_ref = False):
        print "\tdoFastaLookupTest ({})".format(feature_type)
        c = Client()
        
        base_url = "/tark/location/human/{}/?start={}&end={}&region={}"

        if write_ref:
            refs = {}
        else:
            json_file = open( os.path.join( REF_FILE_PATH, "{}_location_fasta_reference.json".format(feature_type) ) )
            refs = json.load(json_file)

        for feature in features:
            if write_ref:
                refs[self.location_to_str(feature)] = {}

            url = base_url.format(feature_type, feature['start'], feature['end'], feature['region'])

            response = c.get(url, {}, HTTP_ACCEPT='text/x-fasta')
            content = ''.join(self.fetchContent(response))
            if write_ref:
                refs[self.location_to_str(feature)] = content
            else:
                self.assertEquals(refs[self.location_to_str(feature)], content, "Retrieved {} for {}:{}-{} (fasta) don't match".format(feature_type, feature['region'], feature['start'], feature['end']))

        if write_ref:
            self.writeRef( "{}_location_fasta_reference.json".format(feature_type), refs )


    def doJsonLookupTest(self, feature_type, features, assemblies, write_ref = False):
        print "\tdoJsonLookupTest ({})".format(feature_type)
        c = Client()
        
        base_url = "/tark/location/human/{}/?start={}&end={}&region={}"

        if write_ref:
            refs = {}
        else:
            json_file = open( os.path.join( REF_FILE_PATH, "{}_location_json_reference.json".format(feature_type) ) )
            refs = json.load(json_file)

        for feature in features:
            if write_ref:
                refs[self.location_to_str(feature)] = {}

            url = base_url.format(feature_type, feature['start'], feature['end'], feature['region'])

            response = c.get(url)
            content = json.loads(''.join(self.fetchContent(response)))
            if write_ref:
                refs[self.location_to_str(feature)]['ref'] = content
            else:
                self.assertEquals( refs[self.location_to_str(feature)]['ref'], content, "Retrieved {} for {}:{}-{} don't match".format(feature_type, feature['region'], feature['start'], feature['end']) )

            response = c.get(url + '&expand=1')
            content = json.loads(''.join(self.fetchContent(response)))
            if write_ref:
                refs[self.location_to_str(feature)]['ref_expanded'] = content
            else:
                self.assertEquals( refs[self.location_to_str(feature)]['ref_expanded'], content, "Retrieved {} for {}:{}-{} (expanded) don't match".format(feature_type, feature['region'], feature['start'], feature['end']) )

            response = c.get(url + '&skip_sequence=1')
            content = json.loads(''.join(self.fetchContent(response)))
            if write_ref:
                refs[self.location_to_str(feature)]['ref_skip_sequence'] = content
            else:
                self.assertEquals( refs[self.location_to_str(feature)]['ref_skip_sequence'], content, "Retrieved {} for {}:{}-{} (skip_sequence) don't match".format(feature_type, feature['region'], feature['start'], feature['end']) )

            response = c.get(url + '&skip_sequence=1&expand=1')
            content = json.loads(''.join(self.fetchContent(response)))
            if write_ref:
                refs[self.location_to_str(feature)]['ref_skip_sequence_expanded'] = content
            else:
                self.assertEquals( refs[self.location_to_str(feature)]['ref_skip_sequence_expanded'], content, "Retrieved {} for {}:{}-{} (skip_sequence, expanded) don't match".format(feature_type, feature['region'], feature['start'], feature['end']) )

            for assembly in assemblies or []:
                response = c.get(url + '&assembly={}'.format(assembly))
                content = json.loads(''.join(self.fetchContent(response)))
                if not write_ref:
                    record = self.fetchRecord(refs[self.location_to_str(feature)]['ref'], u'assembly', assemblies[assembly]);
                    self.assertEquals(record, content, "Retrieved {} for {}:{}-{}, assembly {} don't match".format(feature_type, feature['region'], feature['start'], feature['end'], assembly))

        if write_ref:
            self.writeRef( "{}_location_json_reference.json".format(feature_type), refs )

    def doJsonPostTest(self, feature_type, features, write_ref = False):
        print "\tdoJsonPostTest ({})".format(feature_type)
        c = Client()

        url = "/tark/location/human/{}/".format(feature_type)

        if write_ref:
            refs = {}
        else:
            json_file = open( os.path.join( REF_FILE_PATH, "{}_location_post_json_reference.json".format(feature_type) ) )
            refs = json.load(json_file)

        features_data = {'location': features}

        response = c.post(url, json.dumps(features_data),
                          content_type='application/json')
        content = json.loads(''.join(self.fetchContent(response)))

        if write_ref:
            self.writeRef( "{}_location_post_json_reference.json".format(feature_type), content )
        else:
            self.assertEquals( refs, content, "Retrieved {} via POST don't match".format(feature_type) )


    def fetchContent(self, response):
        if response.status_code != 200:
            return '{"status_code": ' + str(response.status_code) + '}'

        if type(response) == StreamingHttpResponse:
            return response.streaming_content
        else:
            return response.content

    def fetchRecord(self, structure, key, value):
        for rec in structure:
            try:
                if rec and rec.get(key, None) == value:
                    return [rec]
            except Exception as e:
                return structure

        return []

    def location_to_str(self, location):
        return "{}_{}_{}".format(location['region'], location['start'], location['end'])

    def writeRef(self, filename, refs):
        outfile_name = os.path.join( REF_FILE_PATH, filename )
        print "\tWriting reference outfile {}".format(outfile_name)
        with open(outfile_name, 'w') as outfile:
            json.dump(refs, outfile, sort_keys=True, indent=4)
