from django.http import HttpResponse, StreamingHttpResponse
from tark.models import Sequence, FeatureQuerySet, Feature, Assembly, Genenames
from tark.lib.iterables import iterlist
import json
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import types
import pprint
import itertools

ALLOW_CONTENT_TYPE = ['application/json', 'text/plain', 'text/x-fasta' ]

RENDER_FORMAT = {'application/json': 'json',
                 'text/plain': 'plaintext',
                 'text/x-fasta': 'fasta'}

class renderer():
    
    @classmethod
    def render(cls, resultset, **kwargs):

        content_type = kwargs.pop('content-type', 'application/json')

        if content_type not in ALLOW_CONTENT_TYPE:
            print "Unknown content-type {}".format(content_type)
            return HttpResponse(status=500)

        return getattr(cls, RENDER_FORMAT[content_type])(resultset, **kwargs)

    @classmethod
    def json(cls, resultset, **kwargs):
        BioJSONEncoder.serializing_parameters(**kwargs)
        
        if type(resultset) == FeatureQuerySet and resultset:
            iterator = BioJSONEncoder().iterencode(iterlist( resultset.iterator() ))
            return StreamingHttpResponse(iterator, content_type="application/json")
        elif hasattr(resultset, '__iter__') and resultset:
            iterator = BioJSONEncoder().iterencode(iterlist( iter(resultset) ) ) 
            return StreamingHttpResponse(iterator, content_type="application/json")
        else:
            print "HERE"
            data = json.dumps(resultset, indent=4, sort_keys=False, ensure_ascii=False, cls=BioJSONEncoder)
            return HttpResponse(data, content_type="application/json")

    @classmethod
    def fasta(cls, resultset, **kwargs):
        kwargs['expand'] = False
        
        if type(resultset) == FeatureQuerySet:
            iterator = iterlist(resultset.seq_iterator(format='fasta', **kwargs))
        elif hasattr(resultset, '__iter__'):
            iterator = itertools.chain.from_iterable([r.seq_iterator(format='fasta', **kwargs) for r in resultset if type(r) == FeatureQuerySet])
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
    @classmethod
    def serializing_parameters(self, **kwargs):
        self.object_parameters = kwargs
    
    def default(self, obj):
        if hasattr(self, 'object_parameters'):
            kwargs = self.object_parameters
        else:
            kwargs = {}
            pprint.pprint(self.object_parameters)

        if isinstance(obj, Seq):
            return str(obj)
        if isinstance(obj, Feature):
            return obj.to_dict(**kwargs)
        if isinstance(obj, FeatureQuerySet):
            return obj.to_dict(**kwargs)
        if isinstance(obj, Assembly):
            return str(obj)
        if isinstance(obj, Genenames):
            return str(obj)
        json.JSONEncoder.default(self, obj)
