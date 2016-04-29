from django.conf.urls import patterns, url
from django.views.generic.base import RedirectView

from tark import views

urlpatterns = patterns('',
                       url(r'^sequence/(?P<seqtype>[\w]+)/(?P<seqid>[\w\.]+)/$', views.getsequence, name='sequence'),
                       url(r'^lookup/(?P<seqtype>[\w]+)/(?P<id>[\w\.]+)/$', views.idlookup_GET, name='idlookup_get'),
                       url(r'^lookup/(?P<seqtype>[\w]+)/$', views.idlookup_POST, name='idlookup_post'),
                       url(r'^checksum/release/(?P<seqtype>[\w]+)/(?P<tag>[\w\.]+)$', views.checksum_release_type, name='checksum_release_type'),
                       url(r'^lookup/name/gene$', views.name_lookup_gene, name='name_lookup_gene'),
                       url(r'^lookup/name/transcript$', views.name_lookup_transcript, name='name_lookup_transcript'),
                       )
