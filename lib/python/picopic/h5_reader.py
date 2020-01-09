import os, sys

import json
import h5py
import numpy as np


from os.path import join
from picopic.reader import Reader

class H5Reader (Reader):
    '''
    ds: dataset. All information, got from single data file
    component: data component (like 'E_z')
    frame: frame in single dataset (determined by fpds - frame-per-dataset)
    row: row in frame (which consists of rows and columns)
    col: column in frame
    dot: dot in frame, column, or row with coords
    '''
    def __init__(self, path, use_cache=False, verbose=False):

        self.file = h5py.File(os.path.join(path, 'data.h5'), 'r')
        config_json = json.loads(self.file['/metadata'].attrs.get('metadata'))
        super(H5Reader, self).__init__(path, config_json, use_cache, verbose)

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


    def __get_ds_range__(self, p_component, p_type, shape):
        path = self.__get_path__(p_component, p_type, shape)
        frange = len(self.file[path].keys())

        return frange


    def get_frame(self, p_component, shape, number):
        ''' get frame by number. Find required dataset automatically '''
        p_type = 'rec'
        path = self.__get_path__(p_component, p_type, shape)
        path += "/{}".format(number)
        size = self.__frame_size__(shape)
        self.__validate_frame__(p_component, shape, number)
        if self.__verbose__:
            sys.stdout.write("Loading data set {}/{}:{}-{}_{}-{}/{}...".format(p_component, p_type, shape[0], shape[2], shape[1], shape[3], number))
            sys.stdout.flush()

        frame = self.file[path][shape[0]:shape[2],shape[1]:shape[3]]

        if self.__verbose__:
            sys.stdout.write('done\n')
            sys.stdout.flush()

        return frame


    def get_col(self, p_component, shape, number):
        ''' get column by number. Find required dataset automatically '''
        p_type = 'col'
        path = self.__get_path__(p_component, p_type, shape)
        path += "/{}".format(number)
        self.__validate_col__(p_component, shape, number)
        col = self.file[path][shape[0]:shape[2]]

        return col


    def get_row(self, p_component, shape, number):
        ''' get row by number. Find required dataset automatically '''
        p_type = 'row'
        path = self.__get_path__(p_component, p_type, shape)
        path += "/{}".format(number)
        self.__validate_row__(p_component, shape, number)
        row = self.file[path][shape[1]:shape[3]]

        return row


    def get_dot(self, p_component, shape, number):
        ''' get dot by number. Find required dataset automatically '''
        p_type = 'dot'
        path = self.__get_path__(p_component, p_type, shape)
        path += "/{}".format(number)
        self.__validate_dot__(p_component, shape, number)
        dot = self.file[path][0]

        return dot


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
