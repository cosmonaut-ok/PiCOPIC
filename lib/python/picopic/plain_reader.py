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
import fnmatch
from os.path import normpath

import numpy as np
from os.path import join # to use "join" for namespaces
from picopic.reader import Reader

class PlainReader (Reader):
    '''
    ds: dataset. All information, got from single data file
    component: data component (like 'E_z')
    frame: frame in single dataset (determined by fpds - frame-per-dataset)
    row: row in frame (which consists of rows and columns)
    col: column in frame
    dot: dot in frame, column, or row with coords
    '''
    def __init__(self, path, use_cache=False, verbose=False):
        self.__data_set_ranges__ = {}
        config_json = {}
        with open(join(path, 'metadata.json')) as json_file:
            config_json = json.loads(json_file.read())

        super(PlainReader, self).__init__(path, config_json, use_cache, verbose)

    #### special/service functions

    def __del__(self):
        return self


    def __enter__(self):
        '''
        to use as
        ```
        with PlainReader('path/to/db.h5') as h5:
            ...
        ```
        '''
        return self


    def __exit__(self, exception_type, exception_value, traceback):
        self.__del__()


    def __get_path__(self, p_component, p_type, shape):
        base_path = super(PlainReader, self).__get_path__(p_component, p_type, shape)
        path = join(self.meta.config_path, *base_path.split(os.path.sep))

        return path


    def __get_file_path__(self, p_component, p_type, shape, filenumber):
        path = self.__get_path__(p_component, p_type, shape)

        return(join(path, "{}.dat".format(filenumber)))


    def __get_file_range__(self, p_component, p_type, shape):
        path = self.__get_path__(p_component, p_type, shape)
        files = fnmatch.filter(os.listdir(path), '[0-9]*.dat')
        files_s = len(files)

        return files_s

    def __get_ds_range__(self, p_component, p_type, shape):
        path = self.__get_path__(p_component, p_type, shape)
        try:
            frange = self.__data_set_ranges__[path]
        except KeyError:
            frange = self.__get_file_range__(p_component, p_type, shape)
            self.__data_set_ranges__[path] = frange

        return frange


    #### external interface functions
    def get_frame(self, p_component, shape, number):
        ''' get frame by number. Find required dataset automatically '''
        p_type = 'rec'
        size = self.__frame_size__(shape)
        self.__validate_ds__(p_component, p_type, shape, number)
        self.__validate_frame__(p_component, shape, number)

        path = self.__get_file_path__(p_component, p_type, shape, number)
        size = self.__frame_size__(shape)

        frame = []
        cache_file_name = "{}_component:{}_shape:{}-{}-{}-{}_number:{}".format(
            p_type, p_component, shape[0], shape[2], shape[1], shape[3], number)

        if self.__use_cache__:
            frame = self.__tiny_cache__.get_cache(cache_file_name)

        if len(frame) != 0:
            frame = np.reshape(frame, [size[0], size[1]])
        else:

            if self.__verbose__:
                sys.stdout.write("Loading data set {}/{}:{}-{}_{}-{}/{}...".format(p_component, p_type, shape[0], shape[2], shape[1], shape[3], number))
                sys.stdout.flush()
            with open(path, 'r', encoding='utf-8') as datafile:
                frame = np.fromfile(datafile, dtype=float,
                                    count=size[0] * size[1],
                                    sep=' ')

            real_shape = [size[0], size[1]]
            frame = np.reshape(frame, real_shape)
            if self.__verbose__:
                sys.stdout.write('done\n')
                sys.stdout.flush()

        if self.__use_cache__:
            self.__tiny_cache__.update_cache(cache_file_name, frame)

        return frame


    def get_col(self, p_component, shape, number):
        ''' get col by number. Find required dataset automatically '''
        p_type = 'col'
        size = self.__frame_size__(shape)
        self.__validate_ds__(p_component, p_type, shape, number)
        self.__validate_col__(p_component, shape, number)
        col = []

        path = self.__get_file_path__(p_component, p_type, shape, number)
        size = self.__fullframe__[0]

        if self.__verbose__:
            sys.stdout.write("Loading data set {}/{}_{}/{}...".format(p_component, p_type, shape[0], number))
            sys.stdout.flush()

        cache_file_name = "{}_component:{}_shape:{}_number:{}".format(
            p_type, p_component, shape[0], number)

        if self.__use_cache__:
            col = self.__tiny_cache__.get_cache(cache_file_name)

        if len(col) != 0:
            col = np.reshape(col, [size[0]]) # TODO: is it needed here?
        else:
            with open(path, 'r', encoding='utf-8') as datafile:
                col = np.fromfile(datafile, dtype=float,
                                   count=size,
                                   sep=' ')

            real_shape = [size]
            col = np.reshape(cols, real_shape)
        if self.__verbose__:
            sys.stdout.write('done\n')
            sys.stdout.flush()

        if self.__use_cache__:
            self.__tiny_cache__.update_cache(cache_file_name, col)

        return col


    def get_row(self, p_component, shape, number):
        ''' get row by number. Find required dataset automatically '''
        p_type = 'row'
        size = self.__frame_size__(shape)
        self.__validate_ds__(p_component, p_type, shape, number)
        self.__validate_row__(p_component, shape, number)
        path = self.__get_file_path__(p_component, p_type, shape, number)
        size = self.__fullframe__[1]
        row = []
        cache_file_name = "{}_component:{}_type:{}_shape:{}_number:{}".format(
            p_type, p_component, shape[0], number)

        if self.__verbose__:
            sys.stdout.write("Loading data set {}/{}_{}/{}...".format(p_component, p_type, shape, number))
            sys.stdout.flush()

        if self.__use_cache__:
            row = self.__tiny_cache__.get_cache(cache_file_name)

        if len(row) != 0:
            row = np.reshape(row, [size[1]]) # TODO: is it needed here?
        else:
            with open(path, 'r', encoding='utf-8') as datafile:
                row = np.fromfile(datafile, dtype=float,
                                  count=size,
                                  sep=' ')

            real_shape = [size]
            row = np.reshape(row, real_shape)
        if self.__verbose__:
            sys.stdout.write('done\n')
            sys.stdout.flush()

        return rows_re


    def get_all_dots_in_ds(self, p_component, shape, number):
        ''' get all dots in dataset.
        helper function only for plain_reader
        to optimize reading from files set for
        external user '''
        p_type = 'dot'
        self.__validate_ds__(p_component, p_type, shape, number)

        path = self.__get_file_path__(p_component, p_type, shape, number)

        if self.__verbose__:
            sys.stdout.write("Loading data set {}/{}_{}/{}...".format(p_component, p_type, shape, number))
            sys.stdout.flush()
        with open(path, 'r', encoding='utf-8') as datafile:
            dots = np.fromfile(datafile, dtype=float,
                                 count=shape[0] * shape[1],
                                 sep=' ')

        if self.__verbose__:
            sys.stdout.write('done\n')
            sys.stdout.flush()

        return dots


    def get_dot(self, p_component, shape, number):
        ''' get dot by number. Find required dataset automatically '''
        p_type = 'dot'
        self.__validate_ds__(p_component, p_type, shape, number)
        self.__validate_dot__(p_component, shape, number)

        path = self.__get_file_path__(p_component, p_type, shape, number)
        dot = None

        if self.__verbose__:
            sys.stdout.write("Loading data set {}/{}_{}/{}...".format(p_component, p_type, shape, number))
            sys.stdout.flush()

        with open(path, 'r', encoding='utf-8') as datafile:
            dot = np.fromfile(datafile, dtype=float,
                              count=shape[0] * shape[1],
                              sep=' ')

        if self.__verbose__:
            sys.stdout.write('done\n')
            sys.stdout.flush()

        return dot


    ###################################################################
    def get_frame_range(self, p_component, shape, from_frame=0, to_frame=None):
        ''' get frame range. Find, using single get_frame method'''
        fsize = self.__frame_size__(shape)
        self.__validate_frame_range__(p_component, shape, from_frame, to_frame)

        if to_frame == None:
            to_frame = self.__get_ds_range__(p_component, 'rec', shape)
        frange = to_frame - from_frame
        frames = np.zeros((frange, fsize[0], fsize[1]))

        for i in range(from_frame, to_frame):
            self.__validate_frame__(p_component, shape, i)
            frames[i-from_frame] = self.get_frame(p_component, shape, i)

        return(frames)


    def get_col_range(self, p_component, number, from_col=0, to_col=None):
        ''' get column range. Find, using single get_col method'''
        csize = self.meta.geometry_grid[0]
        shape =  [0, 0, csize, number]
        self.__validate_col_range__(p_component, shape, from_row, to_row)
        if to_col == None:
            to_col = self.__get_ds_range__(p_component, 'col', [0, 0, 0, number])
        frange = to_col - from_col
        cols = np.zeros((frange, csize))

        for i in range(from_col, to_col):
            cols[i-from_col] = self.get_col(p_component, shape, i)

        return(cols)


    def get_row_range(self, p_component, number, from_row=0, to_row=None):
        ''' get row range. Find, using single get_row method'''
        rsize = self.meta.geometry_grid[1]
        shape =  [0, 0, number, rsize]
        self.__validate_row_range__(p_component, shape, from_row, to_row)

        if to_row == None:
            to_row = self.__get_ds_range__(p_component, 'row', [0, 0, number, 0])
        frange = to_row - from_row
        rows = np.zeros((frange, rsize))

        for i in range(from_row, to_row):
            rows[i-from_row] = self.get_row(p_component, shape, i)

        return(rows)


    def get_dot_range(self, p_component, row_number, col_number, from_dot=0, to_dot=None):
        ''' get dot range. Find, using single get_dot method'''
        if to_dot == None:
            to_dot = self.__get_ds_range__(p_component, 'dot', [row_number, col_number, 0, 0])
        drange = to_dot - from_dot
        shape = [row_number, col_number, 0, 0]
        self.__validate_dot_range__(p_component, shape, from_dot, to_dot)
        dots = np.zeros((drange))

        for i in range(from_dot, to_dot):
            dots[i-from_dot] = self.get_dot(p_component, shape, i)

        return(dots)
