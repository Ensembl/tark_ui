from django.conf.urls import url
from django.views.generic.base import RedirectView
from django.http import Http404
from django.views.decorators.csrf import csrf_exempt

from tark import views

@csrf_exempt
def method_router(request, *args, **kwargs):
    get_view = kwargs.pop('GET', None)
    post_view = kwargs.pop('POST', None)
    if request.method == 'GET' and get_view is not None:
        return get_view(request, *args, **kwargs)
    elif request.method == 'POST' and post_view is not None:
        return post_view(request, *args, **kwargs)
    raise Http404

urlpatterns = [
                url(r'^$', views.index, name='index'),
                url(r'^genome/', views.getgenomes, name='getgenomes'),
                url(r'^assembly/(?P<species>[\w_]+)/', views.getassembly, name='getgenomes'),
                url('^lookup/(?P<seqtype>[\w]+)/$', method_router, {'GET': views.idlookup_GET, 'POST': views.idlookup_POST}, name='idlookup'),
                url(r'^lookup/name/gene/$', method_router, {'GET': views.name_lookup_gene_GET, 'POST': views.name_lookup_gene_POST}, name='name_lookup_gene'),
                url(r'^lookup/name/transcript/$', method_router, {'GET': views.name_lookup_transcript_GET, 'POST': views.name_lookup_transcript_POST}, name='name_lookup_transcript'),
#                       url(r'^lookup/(?P<seqtype>[\w]+)/(?P<id>[\w\.]+)/$', views.idlookup_GET, name='idlookup_get'),
#                       url(r'^lookup/(?P<seqtype>[\w]+)/$', views.idlookup_GET, name='idlookup_get'),
#                       url(r'^lookup/(?P<seqtype>[\w]+)/$', views.idlookup_POST, name='idlookup_post'),
                url(r'^checksum/release/(?P<seqtype>[\w]+)/(?P<tag>[\w\.]+)$', views.checksum_release_type, name='checksum_release_type'),
                url(r'^location/(?P<species>[\w_]+)/(?P<seqtype>[\w]+)/$', views.location_lookup, name='location_lookup'),
                       
                url(r'^diff/(?P<species>[\w]+)/(?P<release1>[\w\.]+)/(?P<release2>[\w\.]+)/$', views.transcript_diff),
                       
                url(r'^lookup/hgvs/(?P<species>[\w]+)/(?P<hgvs>[\w\d+*\:\.\->]+)/', views.hgvs_by_name, name='hgvs_by_name'),
              ]
