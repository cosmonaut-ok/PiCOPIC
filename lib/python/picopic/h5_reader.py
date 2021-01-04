"""
PiCoPiC
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
import os, sys, re

import json
import h5py
import numpy as np


from os.path import join
from picopic.reader import Reader

class H5Reader (Reader):
    '''
    ds: dataset. All information, got from single data file
    component: data component (like 'E_z')
    rec: rec in dataset
    row: row in dataset (which consists of rows and columns)
    col: column in dataset
    dot: dot in frame, column, or row with coords
    '''
    def __init__(self, path, use_cache=False, verbose=False):
        real_path = ""
        if os.path.isfile(path):
            self.file = h5py.File(path, 'r')
            real_path = os.path.dirname(path)
        elif os.path.isdir(path):
            self.file = h5py.File(os.path.join(path, 'data.h5'), 'r')
            real_path = path
        else:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), path)

        config_json = json.loads(self.file['/'].attrs.get('metadata'))
        super(H5Reader, self).__init__(real_path, config_json, use_cache, verbose)

    def __del__(self):
        self.file.close()


    def __enter__(self):
        '''
        to use as
        ```
        with H5Reader('path/to/db.h5') as h5:
            ...
        ```
        '''
        return self


    def __exit__(self):
        self.__del__()


    def __ds_range__(self, p_component, p_type, shape):
        path = self.__path__(p_component, p_type, shape)
        frange = len(self.file[path].keys())

        return frange


    def rec(self, p_component, shape, number):
        '''get rectangle-shaped dataset by number'''
        p_type = 'rec'
        path = self.__path__(p_component, p_type, shape)
        size = self.__rec_size__(shape)
        if self.__verbose__:
            sys.stdout.write("Loading data set {}/{}:{}-{}_{}-{}/{}..."
                             .format(p_component, p_type, shape[0], shape[2], shape[1], shape[3], number))
            sys.stdout.flush()
        realnumber = re.sub("_.*$", "", number) if type(number) == str else number
        self.__validate_rec__(p_component, shape, number)
        rec = self.file[path][number]

        if self.__verbose__:
            sys.stdout.write('done\n')
            sys.stdout.flush()

        return rec


    def col(self, p_component, z, number): # longitude and number by time
        '''get column by number.'''
        p_type = 'col'
        shape = [0, 0, z, 0]
        path = self.__path__(p_component, p_type, shape)
        if self.__verbose__:
            sys.stdout.write("Loading data set {}/{}/{}..."
                             .format(p_component, p_type, shape[2], number))
            sys.stdout.flush()
        self.__validate_col__(p_component, shape, number)
        col = self.file[path][number]
        if self.__verbose__:
            sys.stdout.write('done\n')

        return col


    def row(self, p_component, r, number): # latitude and number by time
        '''get row by number'''
        p_type = 'row'
        shape = [0, 0, 0, r]
        path = self.__path__(p_component, p_type, shape)
        if self.__verbose__:
            sys.stdout.write("Loading data set {}/{}/{}..."
                             .format(p_component, p_type, shape[2], number))
            sys.stdout.flush()
        self.__validate_row__(p_component, shape, number)
        row = self.file[path][number]
        if self.__verbose__:
            sys.stdout.write('done\n')

        return row


    def dot(self, p_component, r, z, number): # latitude, longitude and number by time
        ''' get dot by number. Find required dataset automatically '''
        p_type = 'dot'
        shape = [0, 0, z, r]
        path = self.__path__(p_component, p_type, shape)
        self.__validate_dot__(p_component, shape, number)
        if self.__verbose__:
            sys.stdout.write("Loading data set {}/{}_{}..."
                             .format(p_component, p_type, r, z))
            sys.stdout.flush()
        dot = self.file[path][number]
        if self.__verbose__:
            sys.stdout.write('done\n')

        return dot


    def col_rec(self, p_component, z, rec_shape, number):
        '''get column from rectangle'''
        return self.rec(p_component, rec_shape, number)[:,z]


    def row_rec(self, p_component, r, rec_shape, number):
        '''get row from rectangle'''
        return self.rec(p_component, rec_shape, number)[r]


    def dot_rec(self, p_component, r, z, rec_shape, number):
        '''get dot from rectangle'''
        return(self.rec(p_component, rec_shape, number)[r,z])


    def rec_range(self, p_component, shape, start_number=0, end_number=None):
        '''get rectancles range'''
        rsize = self.__rec_size__(shape)
        if end_number == None:
            end_number = self.__ds_range__(p_component, 'rec', shape)
        rrange = end_number - start_number
        recs = np.zeros((rrange, rsize[0], rsize[1]))

        for i in range(start_number, end_number):
            self.__validate_rec__(p_component, shape, i)
            recs[i-start_number] = self.get_rec(p_component, shape, i)

        return(recs)
