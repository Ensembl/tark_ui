from django.http import HttpResponse
from django.shortcuts import render
from django.apps import apps
from django.views.decorators.csrf import csrf_exempt
import json

from tark.models import FEATURE_TYPES
from tark.decorators import render, parameter_parser

import pprint

#@render(default_content_type='text/x-fasta')
@parameter_parser
@render
def getsequence(request, seqtype, seqid, **kwargs):
    
    if seqtype not in FEATURE_TYPES:
        return HttpResponse(status=403)
    
    model = apps.get_model('tark', FEATURE_TYPES[seqtype])
    feature_set = model.objects.filter(stable_id=seqid)
    
    if not feature_set:
        return HttpResponse(status=404)
    
    return feature_set

#    return renderer.render(feature_set, **render_parameters)

@parameter_parser(allow_methods='GET')
@csrf_exempt
def idlookup_GET(request, seqtype, id, **kwargs):
    
    pprint.pprint(kwargs)

    return HttpResponse("OK")

@parameter_parser(allow_methods='POST')
@csrf_exempt
def idlookup_POST(request, seqtype, **kwargs):
    
    pprint.pprint(kwargs)

    return HttpResponse("OK")
