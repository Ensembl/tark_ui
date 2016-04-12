from django.http import HttpResponse
from django.shortcuts import render
import json
from django.apps import apps
from tark.models import FEATURE_TYPES, sequence
from tark.renderers import renderer
from django.core.serializers.json import DjangoJSONEncoder

import pprint

def getsequence(request, seqtype, seqid):
    
    if seqtype not in FEATURE_TYPES:
        return HttpResponse(status=403)
    
    model = apps.get_model('tark', seqtype)
    features = model.fetch_feature(fetch_children=True, filter_pk=True, **{'stable_id': seqid } )
    
    if not features:
        return HttpResponse(status=403)
    
    return renderer.render(features, **{'content_type': 'application/json'})
    
    data = json.dumps(features, indent=4, sort_keys=False, ensure_ascii=False, cls=DjangoJSONEncoder)
    
    return HttpResponse(data, content_type="application/json")
