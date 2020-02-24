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
    frame: frame in single dataset (determined by fpds - frame-per-dataset)
    row: row in frame (which consists of rows and columns)
    col: column in frame
    dot: dot in frame, column, or row with coords
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


    def __get_path__(self, p_component, p_type, shape):
        if p_type == 'rec':
            path = "/{}/{}/{}-{}_{}-{}".format(p_component, p_type, shape[0], shape[2], shape[1], shape[3])
        elif p_type == 'col':
            path = "/{}/{}/{}".format(p_component, p_type, shape[3])
        elif p_type == 'row':
            path = "/{}/{}/{}".format(p_component, p_type, shape[2])
        elif p_type == 'dot':
            path = "/{}/{}/{}_{}".format(p_component, p_type, shape[0], shape[1])
        else:
            raise TypeError("Incorrect data type {}. Must be rec/col/row/dot".format(p_type))

        return path


    def __frame_size__(self, shape):
        '''
        return frame size by shape
        '''
        return [shape[2] - shape[0], shape[3] - shape[1]]


    def __get_ds_range__(self, p_component, p_type, shape):
        raise NotImplementedError


#### validators
    def __validate_ds__(self, p_component, p_type, shape, number):
        drange = self.__get_ds_range__(p_component, p_type, shape)
        if number < 0:
            raise IndexError('Dataset should not be less, than 0. The value was: {}'
                             .format(number))
        elif drange < number:
            raise IndexError('Dataset should be less, than {}. The value was: {}.'
                             .format(drange, number))
        else:
            return True


    def __validate_frame__(self, p_component, shape, number):
        frange = self.__get_ds_range__(p_component, 'rec', shape)
        if number < 0:
            raise IndexError('Frame should not be less, than 0. The value was: {}'
                             .format(number))
        elif frange <= number:
            raise IndexError('Frame should be less, than {}. The value was: {}'
                             .format(frange - 1, number))
        else:

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


    def __validate_frame_range__(self, p_component, shape, from_frame, to_frame):
        if to_frame < from_frame:
            raise IndexError('from_frame should be less or equal, than to_frame. The values was: {} and {}'
                            .format(from_frame, to_frame))
        else:

            return(self.__validate_frame__(p_component, shape, to_frame)
                   and self.__validate_frame__(p_component, shape, to_frame))


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


    def get_frame(self, p_component, shape, number):
        ''' get frame by number. Find required dataset automatically '''
        raise NotImplementedError


    def get_col(self, p_component, shape, number):
        ''' get column by number. Find required dataset automatically '''
        raise NotImplementedError


    def get_row(self, p_component, shape, number):
        ''' get row by number. Find required dataset automatically '''
        raise NotImplementedError


    def get_dot(self, p_component, shape, number):
        ''' get dot by number. Find required dataset automatically '''
        raise NotImplementedError


    def get_frame_range(self, p_component, shape, from_frame=0, to_frame=None):
        ''' get frame range. Find, using single get_frame method'''
        fsize = self.__frame_size__(shape)
        if to_frame == None:
            to_frame = self.__get_ds_range__(p_component, 'rec', shape)
        frange = to_frame - from_frame
        frames = np.zeros((frange, fsize[0], fsize[1]))

        for i in range(from_frame, to_frame):
            self.__validate_frame__(p_component, shape, i)
            frames[i-from_frame] = self.get_frame(p_component, shape, i)

        return(frames)

    def get_frame_range_col(self, p_component, number, from_frame=0, to_frame=None):
        ''' get column range. Find, using single get_frame method'''
        csize = self.meta.geometry_grid[0]
        if to_frame == None:
            to_frame = self.__get_ds_range__(p_component, 'col', [0, 0, 0, number])
        frange = to_frame - from_frame
        frames = np.zeros((frange, csize))

        for i in range(from_frame, to_frame):
            frames[i-from_frame] = self.get_col(p_component, [0, 0, csize, number], i)

        return(frames)


    def get_frame_range_row(self, p_component, number, from_frame=0, to_frame=None):
        ''' get row range. Find, using single get_frame method'''
        rsize = self.meta.geometry_grid[1]
        if to_frame == None:
            to_frame = self.__get_ds_range__(p_component, 'row', [0, 0, number, 0])
        frange = to_frame - from_frame
        frames = np.zeros((frange, rsize))

        for i in range(from_frame, to_frame):
            frames[i-from_frame] = self.get_row(p_component, [0, 0, number, rsize], i)

        return(frames)


    def get_frame_range_dot(self, p_component, row_number, col_number, from_frame=0, to_frame=None):
        ''' get frame range. Find, using single get_frame method'''
        if to_frame == None:
            to_frame = self.__get_ds_range__(p_component, 'dot', [row_number, col_number, 0, 0])
        frange = to_frame - from_frame
        frames = np.zeros((frange))

        for i in range(from_frame, to_frame):
            frames[i-from_frame] = self.get_dot(p_component, [row_number, col_number, 0, 0], i)

        return(frames)
