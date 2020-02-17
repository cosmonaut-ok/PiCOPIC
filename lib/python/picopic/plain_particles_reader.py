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
import sys
import fnmatch
from os.path import normpath

import numpy as np
from os.path import join # to use "join" for namespaces
from lib.tinycache import TinyCache

class PlainParticlesReader:
    '''
    ds: dataset. All information, got from single data file
    component: data component (like 'E_z')
    frame: frame in single dataset (determined by fpds - frame-per-dataset)
    row: row in frame (which consists of rows and columns)
    col: column in frame
    dot: dot in frame, column, or row with coords
    '''
    def __init__(self, path, fullframe_size=[-1, -1], fpds=1, use_cache=False, verbose=False):
        '''
        fullframe_size is [r_size, z_size]
        '''
        self.__data_path__ = path
        self.__fpds__ = fpds
        self.__data_set_ranges__ = {}
        self.use_cache = use_cache
        self.__tiny_cache__ = TinyCache(join(self.__data_path__, '.cache'))
        self.__verbose__ = verbose

        if (fullframe_size[0] < 0 or fullframe_size[1] < 0):
            raise IndexError("fullframe_size values should be more, than 0. The value was: {}".format(frame_file))

        self.__fullframe__ = fullframe_size # set size of full simulation frame


    #### special/service functions

    def __del__(self):
        return self


    def __enter__(self):
        '''
        to use as
        ```
        with PlainParticlesReader('path/to/db.h5') as h5:
            ...
        ```
        '''
        return self


    def __exit__(self, exception_type, exception_value, traceback):
        self.__del__()


    def get_ds_frame_by_frame(self, frame):
        ''' get dataset and frame in this file by absolute frame number '''
        frame_file = frame // self.__fpds__
        frame_in_last_file = frame % self.__fpds__

        return frame_file, frame_in_last_file


    def __get_path__(self, p_component, p_type, shape):
        if p_type == 'rec':
            path = join(self.__data_path__, "{}/{}/{}-{}_{}-{}".format(p_component, p_type, shape[0], shape[2], shape[1], shape[3]))
        elif p_type == 'col' or p_type == 'row':
            path = join(self.__data_path__, "{}/{}/{}".format(p_component, p_type, shape[0]))
        elif p_type == 'dot':
            path = join(self.__data_path__, "{}/{}/{}_{}".format(p_component, p_type, shape[1], shape[0])) # use shape[1] then [0] because xy translates to zr :)
        else:
            raise TypeError("Incorrect data type {}. Must be rec/col/row/dot".format(p_type))

        return path

    def __get_file_path__(self, p_component, p_type, shape, filenumber):
        path = self.__get_path__(p_component, p_type, shape)

        print(join(path, "{}.dat".format(filenumber)))
        return(join(path, "{}.dat".format(filenumber)))


    def __frame_size__(self, shape):
        '''
        return frame size by shape
        '''
        return [shape[2] - shape[0], shape[3] - shape[1]]


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
        frange = self.__get_ds_range__(p_component, 'rec', shape) * self.__fpds__
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
            raise IndexError('Row should not be less than 0. The value was {}'
                             .format(number))
        else:

            return True


    def __validate_dot__(self, p_component, shape, number):
        if number < 0:
            raise IndexError('Row should not be less than 0. The value was {}'
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


    #### external interface functions

    def get_all_frames_in_ds(self, p_component, shape, number):
        ''' get all frames in dataset.
        helper function only for plain_reader
        to optimize reading from files set for
        external user '''
        p_type = 'rec'
        self.__validate_ds__(p_component, p_type, shape, number)

        path = self.__get_file_path__(p_component, p_type, shape, number)
        size = self.__frame_size__(shape)

        if self.__verbose__:
            sys.stdout.write("Loading data set {}/{}:{}-{}_{}-{}/{}...".format(p_component, p_type, shape[0], shape[2], shape[1], shape[3], number))
            sys.stdout.flush()
        with open(path, 'r', encoding='utf-8') as datafile:
            frames = np.fromfile(datafile, dtype=float,
                                 count=size[0] * size[1] * self.__fpds__,
                                 sep=' ')

        real_shape = [self.__fpds__, size[0], size[1]]
        frames_re = np.reshape(frames, real_shape)
        if self.__verbose__:
            sys.stdout.write('done\n')
            sys.stdout.flush()

        return frames_re


    def get_all_cols_in_ds(self, p_component, shape, number):
        ''' get all cols in dataset.
        helper function only for plain_reader
        to optimize reading from files set for
        external user '''
        p_type = 'col'
        self.__validate_ds__(p_component, p_type, shape, number)

        path = self.__get_file_path__(p_component, p_type, shape, number)
        size = self.__fullframe__[0]

        if self.__verbose__:
            sys.stdout.write("Loading data set {}/{}_{}/{}...".format(p_component, p_type, shape[0], number))
            sys.stdout.flush()
            #  print(path)
        with open(path, 'r', encoding='utf-8') as datafile:
            cols = np.fromfile(datafile, dtype=float,
                                 count=size * self.__fpds__,
                                 sep=' ')

        real_shape = [self.__fpds__, size]
        cols_re = np.reshape(cols, real_shape)
        if self.__verbose__:
            sys.stdout.write('done\n')
            sys.stdout.flush()

        return cols_re


    def get_all_rows_in_ds(self, p_component, shape, number):
        ''' get all rows in dataset.
        helper function only for plain_reader
        to optimize reading from files set for
        external user '''
        p_type = 'row'
        self.__validate_ds__(p_component, p_type, shape, number)

        path = self.__get_file_path__(p_component, p_type, shape, number)
        size = self.__fullframe__[1]

        if self.__verbose__:
            sys.stdout.write("Loading data set {}/{}_{}/{}...".format(p_component, p_type, shape, number))
            sys.stdout.flush()
        with open(path, 'r', encoding='utf-8') as datafile:
            rows = np.fromfile(datafile, dtype=float,
                                 count=size * self.__fpds__,
                                 sep=' ')

        real_shape = [self.__fpds__, size]
        rows_re = np.reshape(rows, real_shape)
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
                                 count=shape[0] * shape[1] * self.__fpds__,
                                 sep=' ')

        if self.__verbose__:
            sys.stdout.write('done\n')
            sys.stdout.flush()

        return dots


    def get_frame(self, p_component, shape, number):
        ''' get frame by number. Find required dataset automatically '''
        p_type = 'rec'
        size = self.__frame_size__(shape)
        self.__validate_frame__(p_component, shape, number)
        frame = []
        cache_file_name = "{}_component:{}_shape:{}-{}-{}-{}_number:{}".format(
            p_type, p_component, shape[0], shape[2], shape[1], shape[3], number)

        if self.use_cache:
            frame = self.__tiny_cache__.get_cache(cache_file_name)

        if len(frame) != 0:
            frame = np.reshape(frame, [size[0], size[1]])
        else:
            frameds, frame_in_ds = self.get_ds_frame_by_frame(number)
            self.__validate_ds__(p_component, p_type, shape, frameds)
            frame = self.get_all_frames_in_ds(p_component, shape, frameds)[frame_in_ds]

        if self.use_cache:
            self.__tiny_cache__.update_cache(cache_file_name, frame)

        return frame


    def get_col(self, p_component, shape, number):
        ''' get col by number. Find required dataset automatically '''
        p_type = 'col'
        size = self.__frame_size__(shape)
        self.__validate_col__(p_component, shape, number)
        col = []
        cache_file_name = "{}_component:{}_shape:{}_number:{}".format(
            p_type, p_component, shape[0], number)

        if self.use_cache:
            col = self.__tiny_cache__.get_cache(cache_file_name)

        if len(col) != 0:
            col = np.reshape(col, [size[0]]) # TODO: is it needed here?
        else:
            colds, col_in_ds = self.get_ds_frame_by_frame(number)
            self.__validate_ds__(p_component, p_type, shape, colds)
            col = self.get_all_cols_in_ds(p_component, shape, colds)[col_in_ds]

        if self.use_cache:
            self.__tiny_cache__.update_cache(cache_file_name, col)

        return col


    def get_row(self, p_component, shape, number):
        ''' get row by number. Find required dataset automatically '''
        p_type = 'row'
        size = self.__frame_size__(shape)
        self.__validate_row__(p_component, shape, number)
        row = []
        cache_file_name = "{}_component:{}_type:{}_shape:{}_number:{}".format(
            p_type, p_component, shape[0], number)

        if self.use_cache:
            row = self.__tiny_cache__.get_cache(cache_file_name)

        if len(row) != 0:
            row = np.reshape(row, [size[1]]) # TODO: is it needed here?
        else:
            rowds, row_in_ds = self.get_ds_frame_by_frame(number)
            self.__validate_ds__(p_component, p_type, shape, rowds)
            row = self.get_all_rows_in_ds(p_component, shape, rowds)[row_in_ds]

        if self.use_cache:
            self.__tiny_cache__.update_cache(cache_file_name, row)

        return row


    def get_dot(self, p_component, shape, number):
        ''' get dot by number. Find required dataset automatically '''
        p_type = 'dot'
        self.__validate_dot__(p_component, shape, number)

        dotds, dot_in_ds = self.get_ds_frame_by_frame(number)
        self.__validate_ds__(p_component, p_type, shape, dotds)
        dot = self.get_all_dots_in_ds(p_component, shape, dotds)[dot_in_ds]

        return dot

    ###################################################################

    def get_frame_range(self, p_component, shape, from_frame=0, to_frame=None):
        ''' get frame range to 3D array '''
        p_type = 'rec'
        size = self.__frame_size__(shape)
        self.__validate_frame_range__(p_component, shape, from_frame, to_frame)

        if not to_frame:
            to_frame = self.__get_ds_range__(space)

        frame_range = to_frame - from_frame
        from_ds, from_frame_in_ds = self.get_ds_frame_by_frame(from_frame)
        to_ds, to_frame_in_ds = self.get_ds_frame_by_frame(to_frame)
        frames = np.empty([frame_range, size[0], size[1]])

        if frame_range >= self.__fpds__:
            # first ds
            frames[0:self.__fpds__ - from_frame_in_ds] = self.get_all_frames_in_ds(p_component, shape, from_ds)[from_frame_in_ds:]
            # last ds
            if to_frame_in_ds > 0: frames[frame_range - to_frame_in_ds - 1:frame_range - 1] = self.get_all_frames_in_ds(p_component, shape, to_ds)[:to_frame_in_ds]

            shifted_frame = self.__fpds__ - from_frame_in_ds # + 1
            for i in range(from_ds + 1, to_ds):
                i_shifted = i - from_ds - 1   # from_ds + 1 -> 0
                k = i_shifted * self.__fpds__ # from_ds + 1 -> 0
                k_1 = (i_shifted + 1) * self.__fpds__ # from_ds + 1 -> 10
                frames[shifted_frame + k:shifted_frame + k_1] = self.get_all_frames_in_ds(p_component, shape, i)[0:self.__fpds__]
        else:
            frames[0:to_frame_in_ds] = self.get_all_frames_in_ds(p_component, shape, from_ds)[from_frame_in_ds:to_frame_in_ds]

        return frames


    def get_frame_range_col(self, p_component, number, from_frame=0, to_frame=None):
        p_type = 'col'
        size = self.__fullframe__[0]
        self.__validate_col_range__(p_component, number, from_frame, to_frame)

        if not to_frame:
            to_frame = self.__get_ds_range__(space)

        cache_file_name = "{}_component:{}_from:{}_to:{}_number:{}".format(
            p_type, p_component, from_frame, to_frame, number)
        frames = np.empty(0)

        if self.use_cache:
            frames = self.__tiny_cache__.get_cache(cache_file_name)

        if len(frames) != 0:
            frames = np.reshape(frames, [to_frame - from_frame + 1, size])
        else:
            frame_range = to_frame - from_frame
            from_ds, from_frame_in_ds = self.get_ds_frame_by_frame(from_frame)
            to_ds, to_frame_in_ds = self.get_ds_frame_by_frame(to_frame)

            frames = np.empty([frame_range, size])

            # first ds
            frames[0:self.__fpds__ - from_frame_in_ds] = self.get_all_cols_in_ds(p_component, [number], from_ds)[from_frame_in_ds:self.__fpds__]
            # last ds
            frames[frame_range - to_frame_in_ds:frame_range] = self.get_all_cols_in_ds(p_component, [number], to_ds)[:to_frame_in_ds]

            shifted_frame = self.__fpds__ - from_frame_in_ds
            for i in range(from_ds + 1, to_ds):
                i_shifted = i - from_ds - 1
                k = i_shifted * self.__fpds__
                k_1 = (i_shifted + 1) * self.__fpds__
                frames[shifted_frame + k:shifted_frame + k_1] = self.get_all_cols_in_ds(p_component, [number], i)[0:self.__fpds__]

            if self.use_cache:
                    self.__tiny_cache__.update_cache(cache_file_name, frames)

        return frames


    def get_frame_range_row(self, p_component, number, from_frame=0, to_frame=None):
        p_type = 'row'
        size = self.__fullframe__[1]
        self.__validate_row_range__(p_component, number, from_frame, to_frame)

        if not to_frame:
            to_frame = self.__get_ds_range__(space)

        cache_file_name = "{}_component:{}_from:{}_to:{}_number:{}".format(
            p_type, p_component, from_frame, to_frame, number)
        frames = np.empty(0)

        if self.use_cache:
            frames = self.__tiny_cache__.get_cache(cache_file_name)

        if len(frames) != 0:
            frames = np.reshape(frames, [to_frame - from_frame, size])
        else:
            frame_range = to_frame - from_frame
            from_ds, from_frame_in_ds = self.get_ds_frame_by_frame(from_frame)
            to_ds, to_frame_in_ds = self.get_ds_frame_by_frame(to_frame)

            frames = np.empty([frame_range, size])

            # first ds
            frames[0:self.__fpds__ - from_frame_in_ds] = self.get_all_rows_in_ds(p_component, [number], from_ds)[from_frame_in_ds:self.__fpds__]
            # last ds
            frames[frame_range - to_frame_in_ds:frame_range] = self.get_all_rows_in_ds(p_component, [number], to_ds)[:to_frame_in_ds]

            shifted_frame = self.__fpds__ - from_frame_in_ds
            for i in range(from_ds + 1, to_ds):
                i_shifted = i - from_ds - 1
                k = i_shifted * self.__fpds__
                k_1 = (i_shifted + 1) * self.__fpds__
                frames[shifted_frame + k:shifted_frame + k_1] = self.get_all_rows_in_ds(p_component, [number], i)[0:self.__fpds__]

            if self.use_cache:
                    self.__tiny_cache__.update_cache(cache_file_name, frames)

        return frames


    def get_frame_range_dot(self, p_component, row_number, col_number, from_frame=0, to_frame=None):
        p_type = 'dot'
        shape = [col_number, row_number]
        self.__validate_dot_range__(p_component, shape, from_frame, to_frame)

        if not to_frame:
            to_frame = self.__get_ds_range__(space)

        cache_file_name = "{}_component:{}_from:{}_to:{}_row:{}_col:{}".format(
            p_type, p_component, p_type, from_frame, to_frame, row_number, col_number)
        frames = np.empty(0)

        if self.use_cache:
            frames = self.__tiny_cache__.get_cache(cache_file_name)

        if len(frames) == 0:
            frame_range = to_frame - from_frame
            frames = np.empty(frame_range)
            from_ds, from_frame_in_ds = self.get_ds_frame_by_frame(from_frame)
            to_ds, to_frame_in_ds = self.get_ds_frame_by_frame(to_frame)

            # first ds
            frames[0:self.__fpds__ - from_frame_in_ds] = self.get_all_dots_in_ds(p_component, shape, from_ds)[from_frame_in_ds:self.__fpds__]

            # last ds
            frames[frame_range - to_frame_in_ds:frame_range] = self.get_all_dots_in_ds(p_component, shape, to_ds)[:to_frame_in_ds]

            shifted_frame = self.__fpds__ - from_frame_in_ds
            for i in range(from_ds + 1, to_ds):
                i_shifted = i - from_ds - 1
                k = i_shifted * self.__fpds__
                k_1 = (i_shifted + 1) * self.__fpds__
                frames[shifted_frame + k:shifted_frame + k_1] = self.get_all_dots_in_ds(p_component, shape, i)[0:self.__fpds__]
            if self.use_cache:
                self.__tiny_cache__.update_cache(cache_file_name, frames)

        return frames
