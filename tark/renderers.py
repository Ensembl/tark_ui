from django.http import HttpResponse
import json
from Bio.Seq import Seq
from django.core.serializers.json import DjangoJSONEncoder

ALLOW_CONTENT_TYPE = ['application/json', 'text/plain', 'text/x-fasta' ]

RENDER_FORMAT = {'application/json': 'json',
                 'text/plain': 'plaintext',
                 'text/x-fasta': 'fasta'}

class renderer():
    
    @classmethod
    def render(cls, featureset, **kwargs):
        
        if 'content-type' in kwargs:
            content_type = kwargs['content-type']
            if content_type not in ALLOW_CONTENT_TYPE:
                return HttpResponse(status=500)
        else:
            content_type = 'application/json'

        return getattr(cls, RENDER_FORMAT[content_type])(featureset, **kwargs)

    @classmethod
    def json(cls, featureset, **kwargs):
        data = json.dumps(featureset, indent=4, sort_keys=False, ensure_ascii=False, cls=BioJSONEncoder)
        
        return HttpResponse(data, content_type="application/json")


class BioJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Seq):
            return str(obj)
        json.JSONEncoder.default(self, obj)