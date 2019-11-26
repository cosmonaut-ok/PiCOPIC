#!/usr/bin/env python

import os
import numpy
from numpy import fromfile

class TinyCache:
    def __init__(self, path):
        '''
        constructor with set of path to cached data
        '''
        self.__cache__ = {}
        self.path = path
        os.makedirs(path, exist_ok=True)

    def __cache_path__(self, key):
        return(os.path.join(self.path, key))


    def __contains__(self, key):
        cache_file = self.__cache_path__(key)
        res = False
        if key in self.__cache__ or os.path.isfile(cache_file):
            res = True
        
        return res


    def get_cache(self, key):
        cache_file = self.__cache_path__(key)
        cache_value = numpy.empty(0)
        
        if self.__contains__(key):
            if not key in self.__cache__:
                self.__cache__[key] = fromfile(cache_file, dtype='float')
                cache_value = self.__cache__[key]
            else:
                cache_value = self.__cache__[key]

        return(cache_value)


    def update_cache(self, key, value, is_file_bound=None):
        cache_file = self.__cache_path__(key)
        
        intvalue = numpy.asarray(value, dtype='float')
        self.__cache__[key] = intvalue
        intvalue.tofile(cache_file)


# In[42]:


tc = TinyCache("/tmp/tinycache")


# In[43]:


tc.update_cache("test", 150)


# In[44]:


tc.get_cache("test")

