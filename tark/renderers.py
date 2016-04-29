from django.http import HttpResponse, StreamingHttpResponse
from tark.models import Sequence
from tark.iterables import iterlist
import json
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pprint

ALLOW_CONTENT_TYPE = ['application/json', 'text/plain', 'text/x-fasta' ]

RENDER_FORMAT = {'application/json': 'json',
                 'text/plain': 'plaintext',
                 'text/x-fasta': 'fasta'}

class renderer():
    
    @classmethod
    def render(cls, featureset, **kwargs):

        pprint.pprint(kwargs)        
        content_type = kwargs.pop('content-type', 'application/json')

        if content_type not in ALLOW_CONTENT_TYPE:
            print "Unknown content-type {}".format(content_type)
            return HttpResponse(status=500)

        return getattr(cls, RENDER_FORMAT[content_type])(featureset, **kwargs)

    @classmethod
    def json(cls, featureset, **kwargs):
#        features = featureset.to_dict( **kwargs )
        iterator = BioJSONEncoder().iterencode(iterlist( featureset.dict_iterator( **kwargs ) ))
#        data = json.dumps(features, indent=4, sort_keys=False, ensure_ascii=False, cls=BioJSONEncoder)
        
        return StreamingHttpResponse(iterator, content_type="application/json")

    @classmethod
    def fasta(cls, featureset, **kwargs):
        kwargs['expand'] = False

        iterator = iterlist(featureset.seq_iterator(format='fasta', **kwargs))
        return StreamingHttpResponse(iterator, content_type="text/x-fasta")
    
    @classmethod
    def render_error(cls, error, **kwargs):
        content_type = kwargs.pop('content-type', 'application/json')
        
        if content_type == 'application/json':
            data = json.dumps({'error': error})
            return HttpResponse(data, content_type="application/json", status=404)
        else:
            return HttpResponse("Error: " + error, content_type='text/plain')

class BioJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Seq):
            return str(obj)
        json.JSONEncoder.default(self, obj)