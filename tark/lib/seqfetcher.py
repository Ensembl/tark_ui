import requests
from urlparse import urljoin
from django.conf import settings
import pprint

class SeqFetcher():

    @classmethod
    def fetch(cls, assembly, region, start, end):
        if cls.has_sequence(assembly, region, end):
            location = "{}:{}-{}".format(region, start, end)
        
            result = cls.fetch_sequence(assembly, location)
        
            # If it's not found, an empty string will be returned
            return result[location] if result[location] else None
        else:
            print "false?"
            
        return None

    
    @classmethod
    def has_sequence(cls, assembly, region, end):
        if not cls.has_archives():
            return None

        regions = cls.available_locations(assembly)
        
        return regions and region in regions and end < regions[region]
    
    @classmethod
    def fetch_sequence(cls, assembly, location):
        if not cls.has_archives():
            return None
        
        url = cls.seq_url() + "?set={}&location={}".format(assembly, location)

        r = requests.get(url, headers={ "Accept" : "application/json"})
        
        if not r.ok:
            return None
        
        decoded = r.json()
        
        return decoded
    
    @classmethod
    def available_sets(cls):
        if not cls.has_archives():
            return None
        
        url = cls.seq_url("sets") 

        r = requests.get(url, headers={ "Accept" : "application/json"})
        
        if not r.ok:
            return None
        
        decoded = r.json()
        
        return decoded
    
    @classmethod
    def available_locations(cls, assembly):
        if not cls.has_archives():
            return None
        
        url = cls.seq_url("locations/{}/".format(assembly)) 

        r = requests.get(url, headers={ "Accept" : "application/json"})
        
        if not r.ok:
            return None
        
        decoded = r.json()
        
        return decoded
    
    @classmethod
    def has_archives(cls):
        return hasattr(settings, 'SEQ_HOST_URI') and settings.SEQ_HOST_URI
    
    @classmethod
    def seq_url(cls, verb=None):
        
        if verb:
            return urljoin(settings.SEQ_HOST_URI, verb)
    
        return settings.SEQ_HOST_URI
    