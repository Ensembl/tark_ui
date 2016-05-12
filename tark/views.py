from django.http import HttpResponse
from django.apps import apps
from django.views.decorators.csrf import csrf_exempt
from django.db.models.sql.datastructures import Join

from tark.models import FEATURE_TYPES, Releaseset, Transcript, Releasetag, Assembly, Tagset, Genenames, Gene, Transcript
from tark.decorators import render, parameter_parser
from tark.lib.hgvsmapper import fetch_by_hgvs

import pprint

@parameter_parser
@render(default_content_type='text/x-fasta')
def getsequence(request, seqtype, seqid, **kwargs):
    
    if seqtype not in FEATURE_TYPES:
        return HttpResponse(status=403)
    
    model = apps.get_model('tark', FEATURE_TYPES[seqtype])
    feature_set = model.objects.filter(stable_id=seqid)
    
    if not feature_set:
        return HttpResponse(status=404)
    
    return feature_set

@parameter_parser(allow_methods='GET')
@render
def idlookup_GET(request, seqtype, id, **kwargs):
    
    pprint.pprint(kwargs)

    model = apps.get_model('tark', FEATURE_TYPES[seqtype])

    feature_set = model.objects.by_stable_id(id, **kwargs)
    if not feature_set:
        return HttpResponse(status=404)
    
    return feature_set


@parameter_parser(allow_methods='POST')
@render
@csrf_exempt
def idlookup_POST(request, seqtype, **kwargs):
    
    model = apps.get_model('tark', FEATURE_TYPES[seqtype])

    if 'ids' in kwargs:
        return model.objects.by_stable_ids(kwargs['ids'], **kwargs) 

@parameter_parser(allow_methods='GET')
@render
def checksum_release_type(request, seqtype, tag, **kwargs):

    release = Releaseset.fetch_set(tag=tag, **kwargs)
    release2 = Releaseset.fetch_set(tag=84, **{'assembly': 'GCA_000001405.14'})
#    release2 = Releaseset.fetch_set(tag=83, **kwargs)
    release_tags = Releasetag.objects.filter(release_id=release2, feature_type=2)
    tag = Tagset.fetch_set(tag='CARS')
    print release
#    release = Releaseset.objects.filter(shortname=tag).all()[0]

    if seqtype not in FEATURE_TYPES:
        return HttpResponse(status=403)
    
    model = apps.get_model('tark', FEATURE_TYPES[seqtype])
    feature_set = model.objects.release(release).release(release2, exclude=True)#.filter(loc_start__lt=10000).checksum()
    #return HttpResponse(feature_set)

    return feature_set
#    pprint.pprint(feature_set.all())
    return HttpResponse(feature_set.query)
    
@parameter_parser(allow_methods='GET')
@render
def name_lookup_gene(request, **kwargs):    

    if 'name' not in kwargs:
        return HttpResponse(status=403)
        
    name = kwargs['name']

    genes = Gene.objects.filter(genenames__name=name).build_filters(**kwargs)
    
    return genes

@parameter_parser(allow_methods='GET')
@render
def name_lookup_transcript(request, **kwargs):    

    if 'name' not in kwargs:
        return HttpResponse(status=403)
        
    name = kwargs['name']

    transcripts = Transcript.objects.filter(gene__genenames__name=name).build_filters(**kwargs)

    return transcripts

@parameter_parser(allow_methods='GET')
@render
def hgvs_by_name(request, species, hgvs, **kwargs):
    match_features = fetch_by_hgvs(hgvs, species, **kwargs)
    
    return match_features
    
    
