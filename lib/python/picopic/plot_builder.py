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
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data


class PlotBuilder:
    def __init__(self, x_plot_size, y_plot_size, z_plot_size=0,
                 fig_color=None, fig_width=10.5, fig_height=7, fig_dpi=100,
                 font_family='sans-serif', font_name='DejaVu Sans', font_size=14,
                 x_ticklabel_start=0, y_ticklabel_start=0, z_ticklabel_start=0,
                 x_ticklabel_end=None, y_ticklabel_end=None, z_ticklabel_end=None,
                 number_x_ticks=None, number_y_ticks=None, number_z_ticks=None, number_cbar_ticks=3,
                 tickbox=False, grid=False, is_invert_y_axe=True, aspect='auto',
                 image_interpolation='nearest', image_cmap='terrain', image_clim=[-1, 1],
                 guess_number_ticks=40):

        self.__subplots__ = {}
        self.__images__ = {}
        self.__figure__ = plt.figure(figsize=[fig_width, fig_height], dpi=fig_dpi, facecolor=fig_color)
        #
        self.__font_size__ = font_size
        self.__font_family__ = font_family
        self.__font_name__ = font_name
        rc('font', **{'family':font_family, font_family:[font_name], 'size':font_size})
        rc('text', usetex=True)

        self.x_plot_size = x_plot_size
        self.y_plot_size = y_plot_size
        self.z_plot_size = z_plot_size

        self.x_ticklabel_start = x_ticklabel_start
        self.y_ticklabel_start = y_ticklabel_start
        self.z_ticklabel_start = z_ticklabel_start

        self.x_ticklabel_end = x_ticklabel_end or x_plot_size
        self.y_ticklabel_end = y_ticklabel_end or y_plot_size
        self.z_ticklabel_end = z_ticklabel_end or z_plot_size

        # calculate accurate ticks number automatically (recommended)
        if x_plot_size > 0 and y_plot_size > 0:
            tick_value = self.__guess_tick_size__(guess_number_ticks)
            self.number_x_ticks = number_x_ticks or int(np.around((self.x_ticklabel_end - self.x_ticklabel_start) / tick_value, 0))
            self.number_y_ticks = number_y_ticks or int(np.around((self.y_ticklabel_end - self.y_ticklabel_start) / tick_value, 0))
            self.number_z_ticks = number_z_ticks or int(np.around((self.z_ticklabel_end - self.z_ticklabel_start) / tick_value, 0))
        else:
            self.number_x_ticks = guess_number_ticks
            self.number_y_ticks = guess_number_ticks
            self.number_z_ticks = guess_number_ticks

        self.tickbox = tickbox
        self.grid = grid
        self.is_invert_y_axe = is_invert_y_axe
        self.aspect = aspect

        self.number_cbar_ticks = number_cbar_ticks
        self.image_interpolation = image_interpolation
        self.image_clim = image_clim
        self.image_cmap = image_cmap


    def get_figure(self):
        return self.__figure__

    def get_subplot(self, name):
        ''' Object collects axes (subplots) and images for quick access
        this function allows to quick access to subplot by name
        '''
        return self.__subplots__[name]

    def get_image(self, name):
        ''' Object collects axes (subplots) and images for quick access
        this function allows to quick access to image by name
        (same, as subplot, to what it bount)
        '''
        return self.__images__[name]

################################################################################

    def __guess_tick_size__(self, guess_number_ticks):
        sizes = [self.x_ticklabel_end - self.x_ticklabel_start,
                 self.y_ticklabel_end - self.y_ticklabel_start,
                 self.z_ticklabel_end - self.z_ticklabel_start]
        max_size = np.max(sizes)
        exp_number = np.floor(np.log10(max_size))
        base_number = int(max_size / pow(10, exp_number))
        unit = (base_number / guess_number_ticks)
        unit = unit * pow(10, exp_number)

        return unit

    def __add_subplot_common__(self,
                               name, number, projection, title=None,
                               x_plot_size=None, y_plot_size=None, z_plot_size=None,
                               x_axe_label='X', y_axe_label='Y', z_axe_label='Z',
                               x_ticklabel_start=None, y_ticklabel_start=None, z_ticklabel_start=None,
                               x_ticklabel_end=None, y_ticklabel_end=None, z_ticklabel_end=None,
                               number_x_ticks=None, number_y_ticks=None, number_z_ticks=None,
                               tickbox=False, grid=None, is_invert_y_axe=None,
                               position=None, aspect=None):
        ''' common method to add any'''
        subplot = self.__figure__.add_subplot(number, projection=projection)
        self.__subplots__[name] = subplot

        #  set aspect
        subplot.set_aspect(aspect or self.aspect)

        # set position of needed
        invert_y_axe = self.is_invert_y_axe if is_invert_y_axe is None else is_invert_y_axe
        if invert_y_axe: subplot.invert_yaxis()

        # set axis properties
        subplot.set_title(title or name)

        # set grid
        subplot.grid(self.grid if grid is None else grid)

        # make all plot numbers scientifically
        subplot.ticklabel_format(style='sci', scilimits=(0, 0))

        # set position if needed
        if position: subplot.set_position(position)

        # set axes labels
        yz_rotation_angle = 45
        subplot.set_xlabel(x_axe_label)
        subplot.set_ylabel(y_axe_label, rotation=yz_rotation_angle)
        if projection == '3d': subplot.set_zlabel(z_axe_label, rotation=yz_rotation_angle)

        # set subplot tickbox
        subplot.spines['top'].set_visible(self.tickbox if tickbox is None else tickbox)
        subplot.spines['right'].set_visible(self.tickbox if tickbox is None else tickbox)

        # ticks, that sets grid dimensions, required for data placement
        x_size = x_plot_size or self.x_plot_size
        y_size = y_plot_size or self.y_plot_size
        z_size = z_plot_size or self.z_plot_size
        x_ticks = number_x_ticks or self.number_x_ticks
        y_ticks = number_y_ticks or self.number_y_ticks
        z_ticks = number_z_ticks or self.number_z_ticks

        x_tlabel_start = x_ticklabel_start or self.x_ticklabel_start
        y_tlabel_start = y_ticklabel_start or self.y_ticklabel_start
        z_tlabel_start = z_ticklabel_start or self.z_ticklabel_start
        x_tlabel_end = x_ticklabel_end or self.x_ticklabel_end
        y_tlabel_end = y_ticklabel_end or self.y_ticklabel_end
        z_tlabel_end = z_ticklabel_end or self.z_ticklabel_end

        x_tick_grid_size = np.linspace(0, x_size, x_ticks + 1)
        y_tick_grid_size = np.linspace(0, y_size, y_ticks + 1)
        if projection == '3d': z_tick_grid_size = np.linspace(0, z_size, z_ticks + 1)

        if x_size > 0: subplot.set_xticks(x_tick_grid_size)
        if y_size > 0: subplot.set_yticks(y_tick_grid_size)
        if projection == '3d' and z_size > 0: subplot.set_zticks(z_tick_grid_size)

        # tick labels, that shows __real__ model space dimensions
        # translates from grid_size
        x_tlabel_range = list(map(lambda a: float("{:.2e}".format(a)), np.linspace(x_tlabel_start, x_tlabel_end, x_ticks + 1)))
        y_tlabel_range = list(map(lambda a: float("{:.2e}".format(a)), np.linspace(y_tlabel_start, y_tlabel_end, y_ticks + 1)))

        if projection == '3d': z_tlabel_range = np.around(np.linspace(z_tlabel_start, z_tlabel_end, z_ticks + 1), 5)

        if x_size > 0: subplot.set_xticklabels(x_tlabel_range)
        if y_size > 0: subplot.set_yticklabels(y_tlabel_range)
        if projection == '3d' and z_size > 0: subplot.set_zticklabels(z_tlabel_range)

        # set label on every 2nd grid
        for label in [x for i, x in enumerate(subplot.xaxis.get_ticklabels()) if i%2 != 0]:
            label.set_visible(False)
        for label in [x for i, x in enumerate(subplot.yaxis.get_ticklabels()) if i%2 != 0]:
            label.set_visible(False)
        if projection == '3d':
            for label in [x for i, x in enumerate(subplot.zaxis.get_ticklabels()) if i%2 != 0]:
                label.set_visible(False)

        return subplot


    def add_subplot_cartesian_2d(self, name, number, title=None,
                                 x_plot_size=None, y_plot_size=None,
                                 x_axe_label='X', y_axe_label='Y',
                                 x_ticklabel_start=None, y_ticklabel_start=None,
                                 x_ticklabel_end=None, y_ticklabel_end=None,
                                 number_x_ticks=None, number_y_ticks=None,
                                 tickbox=None, grid=None, is_invert_y_axe=None,
                                 position=None, aspect=None):
        ''' add 2D subplot, cartesian projection '''

        subplot = self.__add_subplot_common__(name, number, projection=None,
                                              title=title,
                                              x_plot_size=x_plot_size,
                                              y_plot_size=y_plot_size,
                                              x_axe_label=x_axe_label,
                                              y_axe_label=y_axe_label,
                                              x_ticklabel_start=x_ticklabel_start,
                                              y_ticklabel_start=y_ticklabel_start,
                                              x_ticklabel_end=x_ticklabel_end,
                                              y_ticklabel_end=y_ticklabel_end,
                                              number_x_ticks=number_x_ticks,
                                              number_y_ticks=number_y_ticks,
                                              tickbox=tickbox,
                                              grid=grid,
                                              is_invert_y_axe=is_invert_y_axe,
                                              position=position,
                                              aspect=aspect)

        return subplot


#     def add_subplot_cartesian_3d(self, name, number, title=None, x_axe_label='X', y_axe_label='Y', z_axe_label='Z',
#                       tickbox=False, grid=False, position=None):
#         ''' https://matplotlib.org/gallery/mplot3d/subplot3d.html '''
#         subplot = self.__add_subplot_common__(name, number, title or name,
#                                               x_axe_label, y_axe_label, z_axe_label,
#                                               tickbox, grid, position, True,
#                                               projection='3d')

#         return(subplot)


    def add_image(self, subplot_name, data, cmap=None, clim=None, interpolation=None):
        '''
        cmap reference: https://matplotlib.org/examples/color/colormaps_reference.html
        '''
        subplot = self.get_subplot(subplot_name)

        if subplot:
            ####
            data_len = len(data[0])
            data_height = len(data)
            if data_len != self.x_plot_size and self.x_plot_size > 0:
                raise ValueError('data array length is not equal to grid X-dimension {}. The value was {}.'.format(self.x_plot_size, data_len))
            elif data_height != self.y_plot_size and self.y_plot_size > 0:
                raise ValueError('data array height is not equal to grid Y-dimension {}. The value was {}.'.format(self.y_plot_size, data_height))
            else:
                image = subplot.imshow(data,
                                       cmap=cmap or self.image_cmap,
                                       origin='lower',
                                       interpolation=interpolation or self.image_interpolation)
                image.set_clim(clim or self.image_clim)
                image.set_data(data)
                self.__images__[subplot_name] = image

                return image
        else:
            raise 'There is no subplot, named {}'.format(subplot_name)


    def add_colorbar(self, subplot_name, image_name=None, title=None,
                     ticks=[-1, 1], ticklabels=None,
                     font_size=None, size="2%", position="right"):
        image = self.get_image(image_name or subplot_name)
        if not font_size: font_size = self.__font_size__

        if not title: title = subplot_name
        if image:
            ax = self.get_subplot(subplot_name)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes(position, size=size, pad=0.05)

            cbar = self.__figure__.colorbar(image, cax=cax) # , ax = axes)

            __ticks = np.linspace(ticks[0], ticks[1], self.number_cbar_ticks)
            __ticklabels = np.linspace(ticklabels[0], ticklabels[1], self.number_cbar_ticks) if ticklabels else __ticks

            def format_s(x):
                return('%.1e' % x)

            __ticklabels = list(map(format_s, __ticklabels))

            cbar.set_label(title, rotation=45)
            cbar.set_ticks(__ticks)
            cbar.set_ticklabels(__ticklabels)
            cbar.ax.tick_params(labelsize=font_size)

    def show(self):
        self.__figure__.show()

    def redraw(self):
        ''' Redraw figure (can be used for animation and video writing) '''
        self.__figure__.canvas.draw_idle()
        self.__figure__.canvas.flush_events()

    def save(self, file):
        self.__figure__.savefig(file)
