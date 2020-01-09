#!/usr/bin/env python3

import os,sys
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "lib/python"))

import argparse
import numpy as np

import matplotlib.animation as ani

from picopic.plot_builder import PlotBuilder
from picopic.h5_reader import H5Reader
from picopic.plain_reader import PlainReader


def run(config_path, clim_e_r, clim_e_z, rho_beam_scale, video_file=None,
        time_range=None, cmap=None, frame_step=1, dry_run=False, view=False, use_grid=False):

    # set reader
    if os.path.isfile(os.path.join(config_path, "metadata.json")):
        reader = PlainReader(path = config_path, use_cache=use_cache, verbose=False)
    elif os.path.isfile(os.path.join(config_path, "data.h5")):
        reader = H5Reader(path = config_path, use_cache=use_cache, verbose=False)
    else:
        raise EnvironmentError("There is no corresponding data/metadata files in the path " + config_path + ". Can not continue.")

    ##  configuration options
    x_axis_label = r'$\mathit{Z (m)}$'
    y_axis_label = r'$\mathit{R (m)}$'
    cbar_axis_label = r'$\frac{V}{m}$'
    cbar_bunch_density_axis_label = r'$m^{-3}$'

    e_r_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Radial\enspace Component}\enspace(E_r)$'
    e_z_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Longitudal\enspace Component}\enspace(E_z)$'
    rho_beam_plot_name = r'$\mathbf{Electron\enspace Bunch\enspace Density}\enspace (\rho_{bunch})$'

    # calculate/update video file path
    video_file = os.path.join(config_path, 'field_movie.avi') if not video_file else video_file

    # define reader (plain reader used)
    autoselect = True
    use_cache = False

    clim_rho_beam = [0, 1e16] # [-(cfg.bunch_density * el_charge * rho_beam_scale), 0]
    clim_estimation = reader.meta.get_clim_estimation()

    if not clim_e_r: clim_e_r = [-clim_estimation, clim_estimation]
    if not clim_e_z: clim_e_z = [-clim_estimation, clim_estimation]
    if not rho_beam_scale: rho_beam_scale = 1

    # define plot builder
    plot = PlotBuilder(reader.meta.geometry_grid[1],
                       reader.meta.geometry_grid[0],
                       fig_color=reader.meta.figure_color,
                       fig_width=reader.meta.figure_width,
                       fig_height=reader.meta.figure_height,
                       fig_dpi=reader.meta.figure_dpi,
                       font_family=reader.meta.figure_font_family,
                       font_name=reader.meta.figure_font_name,
                       font_size=reader.meta.figure_font_size,

                       x_ticklabel_end=reader.meta.geometry_size[0],
                       y_ticklabel_end=reader.meta.geometry_size[1],
                       tickbox=True, grid=use_grid, is_invert_y_axe=False,
                       aspect='equal', image_interpolation='nearest')

    # add subplots
    plot.add_subplot_cartesian_2d(e_r_plot_name, 312, x_axe_label=x_axis_label, y_axe_label=y_axis_label)
    plot.add_subplot_cartesian_2d(e_z_plot_name, 311, x_axe_label=x_axis_label, y_axe_label=y_axis_label)
    plot.add_subplot_cartesian_2d(rho_beam_plot_name, 313, x_axe_label=x_axis_label, y_axe_label=y_axis_label)

    # add initial image with zeros and colorbar
    initial_image = np.zeros([reader.meta.geometry_grid[0], reader.meta.geometry_grid[1]])

    plot.add_image(e_r_plot_name, initial_image, cmap=cmap, clim=clim_e_r)
    plot.add_colorbar(e_r_plot_name, ticks=clim_e_r, title=cbar_axis_label)

    plot.add_image(e_z_plot_name, initial_image, cmap=cmap, clim=clim_e_z)
    plot.add_colorbar(e_z_plot_name, ticks=clim_e_z, title=cbar_axis_label)

    plot.add_image(rho_beam_plot_name, initial_image, cmap=cmap, clim=clim_rho_beam)
    plot.add_colorbar(rho_beam_plot_name, ticks=clim_rho_beam, title=cbar_bunch_density_axis_label)

    if view: plot.show()

    # dirty hack
    for p in reader.meta.probes:
        if p.component == 'E_r' or p.component == 'E_z' or (p.component == 'rho_beam' and p.specie == 'beam_electrons'):
            dump_interval = p.schedule
            break

    if not time_range:
        start_frame = reader.meta.get_frame_number_by_timestamp(reader.meta.time[0], dump_interval)
        end_frame = reader.meta.get_frame_number_by_timestamp(reader.meta.time[1], dump_interval)
    else:
        if time_range[0] > time_range[1]: raise ValueError("End time should not be less, than start time. The values were: {}, {}".format(time_range[0], time_range[1]))
        if time_range[1] > reader.meta.end_time: raise IndexError("End time is out of simulation range {}. The value was {}".format(reader.meta.end_time, time_range[1]))

        start_frame = reader.meta.get_frame_number_by_timestamp(time_range[0], dump_interval)
        end_frame = reader.meta.get_frame_number_by_timestamp(time_range[1], dump_interval)

    FFMpegWriter = ani.writers['ffmpeg']
    metadata = dict(title='Field Components and Beam Evolution Movie', artist='Matplotlib',
                    comment='Movie support!')
    writer = FFMpegWriter(fps = reader.meta.video_fps,
                          metadata = metadata,
                          codec = reader.meta.video_codec,
                          bitrate = reader.meta.video_bitrate)

    if dry_run: video_file = '/dev/null'
    fig = plot.get_figure()

    with writer.saving(fig, video_file, reader.meta.figure_dpi):
        for i in range(start_frame, end_frame):
            sys.stdout.write('Loading dataset ' + str(i) + ' ... ')
            sys.stdout.flush()
            shape = [0, 0, reader.meta.geometry_grid[0], reader.meta.geometry_grid[1]]
            data_r = reader.get_frame('E_r', shape, i)
            data_z = reader.get_frame('E_z', shape, i)
            data_beam = reader.get_frame('density/beam_electrons', shape, i)

            # add timestamp to each frame
            timestamp = reader.meta.get_timestamp_by_frame_number(i, dump_interval)
            fig.suptitle("Time: {:.2e} s".format(timestamp), x=.85, y=.95)
            if i % frame_step == 0:
                plot.add_image(e_r_plot_name, data_r, cmap=cmap, clim=clim_e_r)
                plot.add_image(e_z_plot_name, data_z, cmap=cmap, clim=clim_e_z)
                plot.add_image(rho_beam_plot_name, data_beam, cmap=cmap, clim=clim_rho_beam)

                if view: plot.redraw()
                if not dry_run: writer.grab_frame()
                print('done')
            else:
                print('skip')


def main():
    ## configure RC properties
    # plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

    ####
    parser = argparse.ArgumentParser(description='Tool for creating movie file from PDP3 modelled data.')
    parser.add_argument('properties_path', metavar='properties_path', type=str,
                        help='Full path to properties.xml')

    parser.add_argument('--video-file', type=str,
                        help='Full path to output video file. Default <path/to/parameters.xml>/field_movie.avi',
                        default=None)

    parser.add_argument('--time-range', type=str, help='Time range', default=None)

    parser.add_argument('--cmap', type=str,
                        help='''Use custom colormap. Default %s.
                        Reference: https://matplotlib.org/examples/color/colormaps_reference.html''' % 'gray',
                        default=None)

    parser.add_argument('--beam-scale-factor', type=float,
                        help='''Beam density setting automatically, but you can set scale factor to sets,
                        where initial bunch density should be placed in color range''',
                        default=1)

    parser.add_argument('--frame-step', type=int,
                        help='Use only every Nth frame for movie creation. Default: 1',
                        default=1)

    parser.add_argument('--clim-e-r', type=str,
                        help='Color limit range for Electrical field longitual component', default=None)

    parser.add_argument('--clim-e-z', type=str,
                        help='Color limit range for Electrical field radial component', default=None)

    parser.add_argument('--dry-run', action='store_true', help='Do not write anything. Just for debug', default=False)

    parser.add_argument('--view', action='store_true', default=False,
                        help='View animation as well as write it to file')

    parser.add_argument('--with-grid', action='store_true', help='Use tick grid for plots', default=False)

    args = parser.parse_args()

    time_range=None

    # check if config file exists
    if args.time_range:
        time_range = list(map(float, args.time_range.split(':')))

    clim_e_r = list(map(float, args.clim_e_r.split(':'))) if args.clim_e_r else None
    clim_e_z = list(map(float, args.clim_e_z.split(':'))) if args.clim_e_r else None
    run(args.properties_path,
        clim_e_r=clim_e_r,
        clim_e_z=clim_e_z,
        rho_beam_scale=args.beam_scale_factor,
        video_file=args.video_file,
        time_range=time_range,
        frame_step=args.frame_step,
        cmap=args.cmap,
        dry_run=args.dry_run,
        view=args.view,
        use_grid=args.with_grid)
    if args.view:
        input("Press 'Return' to exit ")


## call main function
if __name__ == "__main__":
    main()

# run(config_file, clim_e_r, clim_e_z, clim_rho_beam, video_file='field_movie.avi',
#     time_range=None, cmap=None, dry-run=False, view=False, use_grid=False):
# input("Press 'Return' to exit ")
