#!/usr/bin/env python
"""
PiCOPIC
Copyright (C) 2020 Alexander Vynnyk

This file is part of picopic python data processing helper library.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

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

