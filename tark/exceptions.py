'''
Created on 15 Apr 2016

@author: lairdm
'''

class FilterNotFound(LookupError):
        def __init__(self, value):
            self.value = value
        def __str__(self):
            return repr(self.value)

class AssemblyNotFound(Exception):
        def __init__(self, value):
            self.value = value
        def __str__(self):
            return repr(self.value)

class ReleaseNotFound(Exception):
        def __init__(self, value):
            self.value = value
        def __str__(self):
            return repr(self.value)

class FeatureNotFound(Exception):
        def __init__(self, value):
            self.value = value
        def __str__(self):
            return repr(self.value)

class MissingTranscriptError(LookupError):
    "No transcript has been associated with the mapper"

class FeatureNonCoding(LookupError):
    "This transcript doesn't have a translation"

class IncompatibleFeatureType(LookupError):
    "This feature isn't compatible with the operation requested"
    
class LocationNotInFeature(ValueError):
    "The location isn't contained in the given feature"
    
class BadLocationCoordinates(ValueError):
    "The location coordinates are either incorrectly formatted or otherwise invalid"

class UnknownCoordinateSystem(ValueError):
    "The coordinate system isn't known"
    
class WrongStrandError(TypeError):
    "The strand doesn't match between the given features"
