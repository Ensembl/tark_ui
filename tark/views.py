from django.http import HttpResponse
from django.apps import apps
from django.views.decorators.csrf import csrf_exempt
from django.db.models.sql.datastructures import Join
from django.db.models import Q
from Bio.SeqFeature import FeatureLocation

from tark.models import FEATURE_TYPES, Releaseset, Transcript, Releasetag, Assembly, Tagset, Genenames, Gene, Transcript,\
    TranscriptReleaseTag, Genome
from tark.decorators import render, parameter_parser
from django.shortcuts import render as render_html
from tark.lib.hgvsmapper import fetch_by_hgvs

import pprint
from itertools import groupby
import json

def index(request):
        return render_html(request, 'index.html')

@render(default_content_type='application/json')
def getgenomes(request, **kwargs):
    """Return all available genomes"""
    
    genomes_obj = []
    
    genomes = Genome.objects.all()
    for genome in genomes:
        genomes_obj.append({'name': genome.name, 'taxonomy': genome.tax_id})

    pprint.pprint(genomes_obj)
        
    return genomes_obj

@render(default_content_type='application/json')
def getassembly(request, species, **kwargs):
    """Return all assemblies available for a given species"""
    
    assemblies_obj = []
    
    try:
        assemblies = Assembly.fetch_by_name(species)
        for assembly in assemblies:
            assemblies_obj.append({'genome': str(assembly.genome),
                                   'name': assembly.assembly_name,
                                   'accession': assembly.assembly_accession,
                                   'version': assembly.assembly_version })
    except:
        pass
    
    return assemblies_obj

@parameter_parser(allow_methods='GET')
@render
def idlookup_GET(request, seqtype, id, **kwargs):

    pprint.pprint(kwargs)

    model = apps.get_model('tark', FEATURE_TYPES[seqtype])

    feature_set = model.objects.by_stable_id(id, **kwargs)
    
    return feature_set


@parameter_parser(allow_methods='POST')
@render
def idlookup_POST(request, seqtype, **kwargs):
    results = []
    print "kwargs"
    pprint.pprint(kwargs)

    model = apps.get_model('tark', FEATURE_TYPES[seqtype])

    try:
        ids = kwargs.pop('id', [])
        if not isinstance(ids, (list, tuple)):
            ids = [ids]
        for id in ids:
            feature_set = model.objects.by_stable_id(id, **kwargs)
                    
            if feature_set:
                results.append(feature_set)
            else:
                results.append([{'stable_id': id, 'message': 'Not found'}])
        
    except:
        return HttpResponse(status=500)
    

    return results

# Needs species restriction
@parameter_parser(allow_methods='GET')
@render
def checksum_release_type(request, seqtype, tag, **kwargs):

    if seqtype not in FEATURE_TYPES:
        return HttpResponse(status=403)
    
    model = apps.get_model('tark', FEATURE_TYPES[seqtype])

    releasetag = model.releasetag()
    release = Releaseset.fetch_set(tag=tag, **kwargs)
#    release2 = Releaseset.fetch_set(tag=84, **{'assembly': 'GCA_000001405.14'})
    release2 = Releaseset.fetch_set(tag=83, **kwargs)
#    release_tags = Releasetag.objects.filter(release_id=release2, feature_type=2)
    tag = Tagset.fetch_set(tag='CARS')
    print release
#    release = Releaseset.objects.filter(shortname=tag).all()[0]

    if seqtype not in FEATURE_TYPES:
        return HttpResponse(status=403)
    
    model = apps.get_model('tark', FEATURE_TYPES[seqtype])
    feature_set = model.objects.release(release).release(release2, exclude=True)#.checksum()#.filter(loc_start__lt=10000).checksum()
    #return HttpResponse(feature_set)
    print feature_set.query
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

# Need to check the arguments to ensure they're valid (ie positive numbers for start/end)
@parameter_parser(allow_methods='GET')
@render
def location_lookup(request, species, seqtype, **kwargs):
    
    pprint.pprint(kwargs)

    print "species: {}".format(species)
    assemblies = Assembly.fetch_by_name(species)    
    print "here"
    model = apps.get_model('tark', FEATURE_TYPES[seqtype])
    
    if 'start' not in kwargs or 'end' not in kwargs or 'region' not in kwargs:
        return HttpResponse(status=403)
    
    start = int(kwargs['start'])
    end = int(kwargs['end'])
    
    features = model.objects.fetch_by_location( kwargs['region'], FeatureLocation(start, end, ref_db='genomic' ), assemblies=assemblies, **kwargs )

    return features

# A testing function for finding the difference between transcript
# sets, not meant for production
@parameter_parser(allow_methods='GET')
@render
def transcript_diff(request, species, release1, release2, **kwargs):

    print "species: {}".format(species)
    print "release1: {}".format(release1)
    print "release2: {}".format(release2)
    assembly1 = Assembly.fetch_by_accession('GCA_000001405.20')    
    assembly2 = Assembly.fetch_by_accession('GCA_000001405.14')    
    print assembly1
    releaseset1 = Releaseset.objects.get(shortname=release1, assembly=assembly1)
    print "HERE!!!"
    releaseset2 = Releaseset.objects.get(shortname=release2, assembly=assembly2)
    print "HERE!!!"
    tags1 = TranscriptReleaseTag.objects.filter(release=releaseset1).values_list("feature", flat=True)
    tags2 = TranscriptReleaseTag.objects.filter(release=releaseset2).values_list("feature", flat=True)
    q1 = Q(transcript_id__in=tags1)
    q2 = Q(transcript_id__in=tags2)
    
    features = Transcript.objects.filter(q1 | q2).order_by('stable_id')
    print features.query
    
    f = open('/tmp/features.txt', 'w')
    for key, group in groupby(features, lambda x: x.stable_id):
        count = 0
        print >> f, "stable_id: {}".format(key)
        for feat in group:
            print >> f, feat.__dict__
            count += 1
        print >> f, "Count: {}".format(count)
    
    #return features

@parameter_parser(allow_methods='GET')
@render
def hgvs_by_name(request, species, hgvs, **kwargs):
    # This is a hack for now... we'll add a table of allowed regions by species later
    # and obviously this doesn't catch X, Y, MT
    (location, rest) = hgvs.split(':', 1)
    try:
        int(location)
        hgvs = 'chr' + hgvs
    except:
        pass
    
    match_features = fetch_by_hgvs(hgvs, species, **kwargs)
    
    return match_features
    
    
