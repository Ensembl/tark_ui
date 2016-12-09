from django.http import HttpResponse
from django.apps import apps
from django.views.decorators.csrf import csrf_exempt
from django.conf import settings
#from django.db.models.sql.datastructures import Join
from django.db.models import Q
from Bio.SeqFeature import FeatureLocation

from tark.models import FEATURE_TYPES, Releaseset, Transcript, Releasetag, Assembly, Tagset, Genenames, Gene, Transcript,\
    TranscriptReleaseTag, Genome, FeatureQuerySet, ReleaseDifference
from .exceptions import BadLocationCoordinates
from tark.decorators import render, parameter_parser
from django.shortcuts import render as render_html
from tark.lib.hgvsmapper import fetch_by_hgvs

import pprint
from itertools import groupby

def index(request):
        return render_html(request, 'index.html')

@render(default_content_type='application/json')
def getgenomes(request, **kwargs):
    """Return all available genomes
       url: /genomes/
    """
    
    genomes_obj = []
    
    genomes = Genome.objects.all()
    for genome in genomes:
        genomes_obj.append({'name': genome.name, 'taxonomy': genome.tax_id})
       
    return genomes_obj

@render(default_content_type='application/json')
def getassembly(request, species, **kwargs):
    """Return all assemblies available for a given species
       url: /assembly/{species}/
    """
    
    assemblies_obj = []
    
    try:
        assemblies = Assembly.fetch_by_name(species)
        for assembly in assemblies:
            aliases = []
            for alias in assembly.aliases.all():
                aliases.append(alias.alias)

            assemblies_obj.append({'genome': str(assembly.genome),
                                   'name': assembly.assembly_name,
                                   'aliases': aliases })

    except Exception as e:
        if settings.DEBUG:
            print str(e)
        pass
    
    return assemblies_obj

@parameter_parser(allow_methods='GET')
@render
def idlookup_GET(request, seqtype, id, **kwargs):
    """ Lookup a gene, transcript, exon, or translation by stable id
        url: /lookup/{feature_type}/
        method: GET
    """

    pprint.pprint(kwargs)

    model = apps.get_model('tark', FEATURE_TYPES[seqtype])

    feature_set = model.objects.by_stable_id(id, **kwargs)
    
    return feature_set


@parameter_parser(allow_methods='POST')
@render
def idlookup_POST(request, seqtype, **kwargs):
    """ Lookup one or more  genes, transcripts, exons, or translations by stable id
        url: /lookup/{feature_type}/
        method: POST
    """

    results = []

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
        
    except Exception as e:
        if settings.DEBUG:
            print str(e)
        return HttpResponse(status=500)
    

    return results

@parameter_parser(allow_methods='GET')
@render
def name_lookup_gene_GET(request, **kwargs):
    """ Lookup a gene by HGNC name or alias (ie. IRAK4)
        url: /lookup/name/gene/
        method: GET
    """

    name = kwargs.pop('name', None)
    if not name:
        return HttpResponse(status=400)

    genes = Gene.objects.by_name(name, **kwargs)
#    genes = Gene.objects.filter(genenames__name=name).build_filters(**kwargs)
    
    return genes

@parameter_parser(allow_methods='POST')
@render
def name_lookup_gene_POST(request, **kwargs):
    """ Lookup one or more genes by HGNC name or alias (ie. IRAK4)
        url: /lookup/name/gene/
        method: POST
    """

    results = []
    
    try:
        names = kwargs.pop('name', [])
        if not isinstance(names, (list, tuple)):
            names = [names]
        for name in names:
            feature_set = Gene.objects.by_name(name, **kwargs)
            
            if feature_set:
                results.append(feature_set)
            else:
                results.append([{'name': name, 'message': 'Not found'}])
    except Exception as e:
        if settings.DEBUG:
            print str(e)
        return HttpResponse(status=500)

    return results

@parameter_parser(allow_methods='GET')
@render
def name_lookup_transcript_GET(request, **kwargs):    
    """ Lookup a transcript by HGNC name or alias (ie. IRAK4)
        url: /lookup/name/transcript/
        method: GET
    """

    name = kwargs.pop('name', None)
    if not name:
        return HttpResponse(status=400)

    transcripts = Transcript.objects.by_name(name, **kwargs)

    return transcripts

@parameter_parser(allow_methods='POST')
@render
def name_lookup_transcript_POST(request, **kwargs):
    """ Lookup one or more transcripts by HGNC name or alias (ie. IRAK4)
        url: /lookup/name/transcript/
        method: POST
    """

    results = []
    
    try:
        names = kwargs.pop('name', [])
        if not isinstance(names, (list, tuple)):
            names = [names]
        for name in names:
            feature_set = Transcript.objects.by_name(name, **kwargs)
            
            if feature_set:
                results.append(feature_set)
            else:
                results.append([{'name': name, 'message': 'Not found'}])
    except Exception as e:
        if settings.DEBUG:
            print str(e)
        return HttpResponse(status=500)

    return results

# Need to check the arguments to ensure they're valid (ie positive numbers for start/end)
@parameter_parser(allow_methods='GET')
@render
def location_lookup_GET(request, species, seqtype, **kwargs):
    """ Lookup all features of a given type (gene, transcript, exon, translation) in
        a particular region.
        url: /location/{species}/{feature_type}/
        method: GET
    """
    return location_lookup(request, species, seqtype, **kwargs)
    
def location_lookup(request, species, seqtype, **kwargs):
    """ Location lookup helper functions to do the main work.
    
    """

    assemblies = Assembly.fetch_by_name(species)    
    model = apps.get_model('tark', FEATURE_TYPES[seqtype])
    
    if 'start' not in kwargs or 'end' not in kwargs or 'region' not in kwargs:
        return HttpResponse(status=400)
    
    try:
        start = int(kwargs['start'])
        end = int(kwargs['end'])

        if start < 0 or end < 0 or end < start:
            raise BadLocationCoordinates()

        features = model.objects.by_location( kwargs['region'], FeatureLocation(start, end, ref_db='genomic' ), assemblies=assemblies, **kwargs )

    except Exception as e:
        if settings.DEBUG:
            print str(e)
        return HttpResponse(status=400)

    return features

@parameter_parser(allow_methods='POST')
@render
def location_lookup_POST(request, species, seqtype, **kwargs):
    """ Lookup all features of a given type (gene, transcript, exon, translation) in
        one or more region(s).
        url: /location/{species}/{feature_type}/
        method: POST
    """

    results = []

    try:
        locations = kwargs.pop('location', [])
        if not isinstance(locations, (list, tuple)):
            locations = [locations]
        for location in locations:
            try:
                if 'start' not in location or 'end' not in location or 'region' not in location:
                    raise BadLocationCoordinates()
                
                feature_set = location_lookup(request, species, seqtype, 
                                              start=location['start'],
                                              end=location['end'],
                                              region=location['region'],
                                              **kwargs)
                
                if feature_set and isinstance(feature_set, FeatureQuerySet):
                    results.append(feature_set)
                else:
                    results.append([])

            except Exception as e:
                if settings.DEBUG:
                    print str(e)
                results.append([])
    except Exception as e:
        if settings.DEBUG:
            print str(e)
        return HttpResponse(status=500)

    return results

def sequence_GET(request, species, assembly):
    """ Fetch an arbitrary section of sequence.
    
    """

    assembly = Assembly.fetch_by_name(species, assembly=assembly)
    
    if not assembly:
        return HttpResponse(status=400)
    
    

def sequence_POST(request, species, assembly):
    pass

# Needs species restriction
@parameter_parser(allow_methods='GET')
@render
def checksum_release_type(request, seqtype, tag, **kwargs):

    if seqtype not in FEATURE_TYPES:
        return HttpResponse(status=400)
    
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
        return HttpResponse(status=400)
    
    model = apps.get_model('tark', FEATURE_TYPES[seqtype])
    feature_set = model.objects.release(release).release(release2, exclude=True)#.checksum()#.filter(loc_start__lt=10000).checksum()
    #return HttpResponse(feature_set)
    print feature_set.query
    return feature_set
#    pprint.pprint(feature_set.all())
    return HttpResponse(feature_set.query)

@parameter_parser(allow_methods='GET')
#@render
def release_diff(request, species, release1, release2, **kwargs):
    
    print "species: {}".format(species)
    print "release1: {}".format(release1)
    print "release2: {}".format(release2)
    assembly1 = Assembly.fetch_by_accession('GCA_000001405.22')    
    assembly2 = Assembly.fetch_by_accession('GCA_000001405.22')    
#    assembly2 = Assembly.fetch_by_accession('GCA_000001405.14')    
    print assembly1
    releaseset1 = Releaseset.objects.get(shortname=release1, assembly=assembly1)
    print "releaseset1: {}".format(releaseset1.shortname)
    releaseset2 = Releaseset.objects.get(shortname=release2, assembly=assembly2)
    print "releaseset2: {}".format(releaseset2.shortname)

    f = open('/tmp/diff_set.txt', 'w')

    for gene_set in ReleaseDifference.compare(release1, release2, assembly1.assembly_id, assembly2.assembly_id): 
#        print >> f, "stable_id: {}, count: {}".format(gene_set.stable_id, gene_set.gene_count)
#        if len(gene_set.gene_set) > 1 and gene_set.gene_set[0][3] != gene_set.gene_set[1][3]:
        if gene_set.is_different:
            print >> f, gene_set.stable_id
            pprint.pprint(gene_set.gene_set, f)
            print gene_set.stable_id
            pprint.pprint(gene_set.gene_set)
            print >> f, "DIFFERENCE"
            pprint.pprint(gene_set.difference(), f)
            gene1 = Gene.objects.by_stable_id("{}.{}".format(gene_set.stable_id, gene_set.gene_set[0][0]), assembly_id=gene_set.gene_set[0][2])
            pprint.pprint(gene1.to_dict(expand=True), f)
            gene2 = Gene.objects.by_stable_id("{}.{}".format(gene_set.stable_id, gene_set.gene_set[1][0]), assembly_id=gene_set.gene_set[1][2])
            pprint.pprint(gene2.to_dict(expand=True), f)

    

# A testing function for finding the difference between transcript
# sets, not meant for production
@parameter_parser(allow_methods='GET')
@render
def transcript_diff(request, species, release1, release2, **kwargs):

    print "species: {}".format(species)
    print "release1: {}".format(release1)
    print "release2: {}".format(release2)
    assembly1 = Assembly.fetch_by_accession('GCA_000001405.22')    
#    assembly2 = Assembly.fetch_by_accession('GCA_000001405.22')    
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
    
    
