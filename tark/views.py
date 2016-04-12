from django.http import HttpResponse
from django.shortcuts import render
from django.apps import apps
from tark.models import FEATURE_TYPES
from tark.renderers import renderer

def getsequence(request, seqtype, seqid):
    
    if seqtype not in FEATURE_TYPES:
        return HttpResponse(status=403)
    
    model = apps.get_model('tark', seqtype)
    feature_set = model.objects.filter(stable_id=seqid)
#    features = feature_set.to_dict()
#    features = feature_set.to_dict( **{'fetch_children': True, 'filter_pk': True} )
    
    if not feature_set:
        return HttpResponse(status=404)
    
    return renderer.render(feature_set, **{'content_type': 'application/json', 'fetch_children': True, 'filter_pk': True})
    
#    data = json.dumps(features, indent=4, sort_keys=False, ensure_ascii=False, cls=DjangoJSONEncoder)
    
#    return HttpResponse(data, content_type="application/json")
