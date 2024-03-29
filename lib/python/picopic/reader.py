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
import os, sys

import json
import h5py
import numpy as np


from os.path import join
from picopic.tinycache import TinyCache
from picopic.meta_reader import MetaReader

class Reader:
    '''
    ds: dataset. All information, got from single data file
    component: data component (like 'E_z')
    rec: rectangle in dataset
    row: row in dataset (which consists of rows and columns)
    col: column in dataset
    dot: dot in dataset, column, or row with coords
    '''
    def __init__(self, path, config_json, use_cache=False, verbose=False):
        config_json['base_path'] = path
        self.meta = MetaReader(config_json)
        self.__verbose__ = verbose
        self.__use_cache__ = use_cache
        if use_cache:
            self.__tiny_cache__ = TinyCache(join(self.meta.config_path, '.cache'))
        else:
            self.__tiny_cache__ = None


    def __path__(self, p_component, p_type, shape):
        if p_type == 'rec':
            path = "/{}/{}/{}-{}_{}-{}".format(p_component, p_type, shape[0], shape[2], shape[1], shape[3])
        elif p_type == 'col':
            path = "/{}/{}/{}".format(p_component, p_type, shape[3])
        elif p_type == 'row':
            path = "/{}/{}/{}".format(p_component, p_type, shape[2])
        elif p_type == 'dot':
            path = "/{}/{}/{}_{}".format(p_component, p_type, shape[3], shape[2])
        else:
            raise TypeError("Incorrect data type {}. Must be rec/col/row/dot".format(p_type))

        return path


    def __rec_size__(self, shape):
        '''
        return rec size by shape
        '''
        return [shape[2] - shape[0], shape[3] - shape[1]]


    def __ds_range__(self, p_component, p_type, shape):
        raise NotImplementedError


#### validators
    def __validate_ds__(self, p_component, p_type, shape, number):
        # drange = self.__ds_range__(p_component, p_type, shape)
        # if number < 0:
        #     raise IndexError('Dataset should not be less, than 0. The value was: {}'
        #                      .format(number))
        # elif drange < number:
        #     raise IndexError('Dataset should be less, than {}. The value was: {}.'
        #                      .format(drange, number))
        # else:
        return True


    def __validate_rec__(self, p_component, shape, number):
        # frange = self.__ds_range__(p_component, 'rec', shape)
        # if number < 0:
        #     raise IndexError('Rec should not be less, than 0. The value was: {}'
        #                      .format(number))
        # elif frange <= number:
        #     raise IndexError('Rec should be less, than {}. The value was: {}'
        #                      .format(frange - 1, number))
        # else:

        return True


    def __validate_col__(self, p_component, shape, number):
        if number < 0:
            raise IndexError('Column should not be less than 0. The value was {}'
                             .format(number))
        else:

            return True

    def __validate_row__(self, p_component, shape, number):
        if number < 0:
            raise IndexError('Column should not be less than 0. The value was {}'
                             .format(number))
        else:

            return True


    def __validate_dot__(self, p_component, shape, number):
        if number < 0:
            raise IndexError('Column should not be less than 0. The value was {}'
                             .format(number))
        else:

            return True


    def __validate_rec_range__(self, p_component, shape, from_rec, to_rec):
        if to_rec < from_rec:
            raise IndexError('from_rec should be less or equal, than to_rec. The values was: {} and {}'
                            .format(from_rec, to_rec))
        else:

            return(self.__validate_rec__(p_component, shape, to_rec)
                   and self.__validate_rec__(p_component, shape, to_rec))


    def __validate_row_range__(self, p_component, shape, from_row, to_row):
        if to_row < from_row:
            raise IndexError('from_row should be less than to_row. The values were {} and {}'
                            .format(from_row, to_row))
        else:

            return (self.__validate_row__(p_component, shape, from_row)
                    and self.__validate_row__(p_component, shape, to_row))


    def __validate_col_range__(self, p_component, shape, from_col, to_col):
        if to_col < from_col:
            raise IndexError('from_col should be less than to_col. The values were {} and {}'
                            .format(from_col, to_col))
        else:

            return (self.__validate_col__(p_component, shape, from_col)
                    and self.__validate_col__(p_component, shape, to_col))


    def __validate_dot_range__(self, p_component, shape, from_dot, to_dot):
        if to_dot < from_dot:
            raise IndexError('from_dot should be less than to_dot. The values were {} and {}'
                            .format(from_dot, to_dot))
        else:

            return (self.__validate_dot__(p_component, shape, from_dot)
                    and self.__validate_dot__(p_component, shape, to_dot))

    ################################################################################

    def rec(self, p_component, shape, number):
        '''get rec by number'''
        raise NotImplementedError


    def col(self, p_component, z, number): # longitude and number by time
        '''get rec by number'''
        raise NotImplementedError


    def row(self, p_component, r, number): # latitude and number by time
        '''get rec by number'''
        raise NotImplementedError


    def dot(self, p_component, r, z, number): # latitude, longitude and number by time
        '''get rec by number'''
        raise NotImplementedError


    def rec_range(self, p_component, shape, start_number=0, end_number=None):
        '''get rec range by numbers'''
        raise NotImplementedError


    def col_range(self, p_component, z, start_number=0, end_number=None):
        '''get col range by numbers'''
        raise NotImplementedError


    def row_range(self, p_component, r, start_number=0, end_number=None):
        '''get row range by numbers'''
        raise NotImplementedError


    def dot_range(self, p_component, r, z, start_number=0, end_number=None):
        '''get dot range by numbers'''
        raise NotImplementedError


    def col_rec(self, p_component, z, rec_shape, number):
        '''get dot range from rec range'''
        raise NotImplementedError


    def row_rec(self, p_component, r, rec_shape, number):
        '''get dot range from rec range'''
        raise NotImplementedError


    def dot_rec(self, p_component, r, z, rec_shape, number):
        '''get dot range from rec range'''
        raise NotImplementedError


    def col_range(self, p_component, z, start_number=0, end_number=None):
        '''get column range'''
        csize = self.meta.geometry_grid[0]
        if end_number == None:
            end_number = self.__ds_range__(p_component, 'col', [0, 0, 0, z])
        crange = start_number - end_number
        cols = np.zeros((crange, csize))

        for i in range(start_number, end_number):
            cols[i-start_number] = self.col(p_component, [0, 0, csize, z], i)

        return(cols)


    def row_range(self, p_component, r, start_number=0, end_number=None):
        '''get row range'''
        rsize = self.meta.geometry_grid[1]
        if end_number == None:
            end_number = self.__ds_range__(p_component, 'row', [0, 0, r, 0])
        rrange = end_number - start_number
        rows = np.zeros((rrange, rsize))

        for i in range(start_number, end_number):
            rows[i-start_number] = self.row(p_component, r, i)

        return(rows)


    def dot_range(self, p_component, r, z, start_number=0, end_number=None):
        '''get dot from rectangle range, using "rec" method'''
        if end_number == None:
            end_number = self.__ds_range__(p_component, 'dot', [0, 0, r, z])
        drange = end_number - start_number
        dots = np.zeros((drange))

        for i in range(start_number, end_number):
            dots[i-start_number] = self.dot(p_component, r, z, i)

        return(dots)


    def col_rec_range(self, p_component, z, rec_shape, start_number=0, end_number=None):
        '''get column range from rec range'''
        rsize = self.meta.geometry_grid[0]
        if end_number == None:
            end_number = self.__ds_range__(p_component, 'rec', shape)
        crange = end_number - start_number
        recs = np.zeros((crange, rsize))
        for i in range(start_number, end_number):
            self.__validate_rec__(p_component, rec_shape, i)
            recs[i-start_number] = self.col_rec(p_component, z, rec_shape, i)

        return(recs)


    def row_rec_range(self, p_component, r, rec_shape, start_number=0, end_number=None):
        '''get row range from rec range'''
        rsize = self.meta.geometry_grid[1]
        if end_number == None:
            end_number = self.__ds_range__(p_component, 'rec', shape)
        rrange = end_number - start_number
        recs = np.zeros((rrange, rsize))
        for i in range(start_number, end_number):
            self.__validate_rec__(p_component, rec_shape, i)
            recs[i-start_number] = self.row_rec(p_component, r, rec_shape, i)

        return(recs)


    def dot_rec_range(self, p_component, r, z, rec_shape, start_number=0, end_number=None):
        '''get dot range from rec range'''
        dsize = self.meta.geometry_grid[1]
        if end_number == None:
            end_number = self.__ds_range__(p_component, 'rec', shape)
        drange = end_number - start_number
        recs = np.zeros((drange))
        for i in range(start_number, end_number):
            self.__validate_rec__(p_component, rec_shape, i)
            recs[i-start_number] = self.dot_rec(p_component, r, z, rec_shape, i)

        return(recs)
