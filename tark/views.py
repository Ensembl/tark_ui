from django.http import HttpResponse
from django.shortcuts import render
from django.apps import apps
from tark.models import FEATURE_TYPES
from tark.decorators import render

#@render(default_content_type='text/x-fasta')
@render
def getsequence(request, seqtype, seqid):
    
    if seqtype not in FEATURE_TYPES:
        return HttpResponse(status=403)
    
    model = apps.get_model('tark', FEATURE_TYPES[seqtype])
    feature_set = model.objects.filter(stable_id=seqid)
    
    if not feature_set:
        return HttpResponse(status=404)
    
    return feature_set

#    return renderer.render(feature_set, **render_parameters)

def idlookup(request, seqtype, seqid):
    pass
