'''
Created on 15 Apr 2016

@author: lairdm
'''

class FilterNotFound(Exception):
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
