#!/usr/bin/env python3

import os,sys
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "lib/python"))

import argparse
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as ani

from picopic.plot_builder import PlotBuilder
from picopic.h5_reader import H5Reader
from picopic.plain_reader import PlainReader


def run(config_path, specie, clim, video_file=None, time_range=None, cmap=None, frame_step=1, frame_size=None, dry_run=False, view=False, use_grid=False, plot_image_interpolation='bicubic'):

    ##  configuration options
    x_axis_label = r'$\mathit{Z (m)}$'
    y_axis_label = r'$\mathit{R (m)}$'
    cbar_axis_label = r'$m^{-3}$'

    plot_name_1 = r'$\mathbf{Beam\enspace Electrons\enspace Density}$'
    plot_name_2 = r'$\mathbf{Ions\enspace Density}$'
    plot_name_3 = r'$\mathbf{Electrons\enspace Density}$'
    plot_name_4 = r'$\mathbf{Ions\enspace Temperature}$'
    plot_name_5 = r'$\mathbf{Electrons\enspace Temperature}$'


    # calculate/update video file path
    video_file = os.path.join(config_path, 'multiplot_movie.avi') if not video_file else video_file

    # define reader (plain reader used)
    autoselect = True
    use_cache = False

    # set reader
    if os.path.isfile(os.path.join(config_path, "metadata.json")):
        reader = PlainReader(path = config_path, use_cache=use_cache, verbose=False)
    elif os.path.isfile(os.path.join(config_path, "data.h5")):
        reader = H5Reader(path = config_path, use_cache=use_cache, verbose=False)
    else:
        raise EnvironmentError("There is no corresponding data/metadata files in the path " + config_path + ". Can not continue.")

    # TODO: estimate clim
    # if not clim: clim = [0,2]
    clim_1 = [0,1e15]
    clim_2 = [4e19,3.5e20]
    clim_3 = [5e16,1.5e17]
    clim_4 = [0,25]
    clim_5 = [0,250]

    if not frame_size: frame_size = [0, 0, reader.meta.geometry_grid[0], reader.meta.geometry_grid[1]]
    frame_src_size=[-1, -1, -1, -1]

    # detect probe shape
    for probe in reader.meta.probes:
        if (probe.shape == 'rec') and (probe.specie == specie) and (probe.size[0] == frame_size[0]) and (probe.size[1] == frame_size[1]) and(probe.size[2] == frame_size[2]) and(probe.size[3] == frame_size[3]):
            frame_src_size = probe.size
    # try bigger frames, if autoselect enabled
    if frame_src_size[0] == -1 and frame_src_size[1] == -1 and frame_src_size[2] == -1 and frame_src_size[3] == -1 and autoselect:
        for probe in reader.meta.probes:
            if (probe.shape == 'rec') and (probe.specie == specie) and (probe.size[0] <= frame_size[0]) and (probe.size[1] <= frame_size[1]) and(probe.size[2] >= frame_size[2]) and(probe.size[3] >= frame_size[3]):
                frame_src_size = probe.size

    r_scale = (frame_size[2] - frame_size[0]) / reader.meta.geometry_grid[0]
    z_scale = (frame_size[3] - frame_size[1]) / reader.meta.geometry_grid[1]
    # define plot builder
    plot = PlotBuilder(frame_size[3] - frame_size[1],
                       frame_size[2] - frame_size[0],
                       fig_color=reader.meta.figure_color,
                       fig_width=reader.meta.figure_width * 2,
                       fig_height=reader.meta.figure_height * 2,
                       fig_dpi=reader.meta.figure_dpi,
                       font_family=reader.meta.figure_font_family,
                       font_name=reader.meta.figure_font_name,
                       font_size=12, # reader.meta.figure_font_size,
                       x_ticklabel_end=reader.meta.geometry_size[1] * z_scale,
                       y_ticklabel_end=reader.meta.geometry_size[0] * r_scale,
                       tickbox=True, grid=use_grid, is_invert_y_axe=False,
                       aspect='equal', image_interpolation=plot_image_interpolation,
                       number_y_ticks=3, guess_number_ticks=20)

    # plt.rcParams['axes.titlesize'] = 10
    # plt.rcParams["axes.titlepad"] = 3

    # add subplots
    plot.add_subplot_cartesian_2d(plot_name_1, 511, x_axe_label="",
                                  y_axe_label=y_axis_label)
    plot.get_subplot(plot_name_1).set_xticklabels([])
    for tic in plot.get_subplot(plot_name_1).xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    plot.add_subplot_cartesian_2d(plot_name_2, 512, x_axe_label="",
                                  y_axe_label=y_axis_label)
    plot.get_subplot(plot_name_2).set_xticklabels([])
    for tic in plot.get_subplot(plot_name_2).xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    plot.add_subplot_cartesian_2d(plot_name_3, 513, x_axe_label="",
                                  y_axe_label=y_axis_label)
    plot.get_subplot(plot_name_3).set_xticklabels([])
    for tic in plot.get_subplot(plot_name_3).xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    plot.add_subplot_cartesian_2d(plot_name_4, 514, x_axe_label="",
                                  y_axe_label=y_axis_label)
    plot.get_subplot(plot_name_4).set_xticklabels([])
    for tic in plot.get_subplot(plot_name_4).xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    plot.add_subplot_cartesian_2d(plot_name_5, 515, x_axe_label=x_axis_label, y_axe_label=y_axis_label)

    # add initial image with zeros and colorbar
    initial_image = np.zeros([frame_size[2] - frame_size[0], frame_size[3] - frame_size[1]])

    cmap_r = "{}_r".format(cmap)
    print(cmap_r)
    plot.add_image(plot_name_1, initial_image, cmap=cmap_r, clim=clim_1)
    plot.add_colorbar(plot_name_1, ticks=clim_1, title=cbar_axis_label)

    plot.add_image(plot_name_2, initial_image, cmap=cmap, clim=clim_2)
    plot.add_colorbar(plot_name_2, ticks=clim_2, title=cbar_axis_label)

    plot.add_image(plot_name_3, initial_image, cmap=cmap, clim=clim_3)
    plot.add_colorbar(plot_name_3, ticks=clim_3, title=cbar_axis_label)

    plot.add_image(plot_name_4, initial_image, cmap=cmap, clim=clim_4)
    plot.add_colorbar(plot_name_4, ticks=clim_4, title=cbar_axis_label)

    plot.add_image(plot_name_5, initial_image, cmap=cmap, clim=clim_5)
    plot.add_colorbar(plot_name_5, ticks=clim_5, title=cbar_axis_label)

    if view: plot.show()

    # dirty hack
    for p in reader.meta.probes:
        if p.component == 'density':
            dump_interval = p.schedule
            break

    if not time_range:
        start_frame = reader.meta.get_frame_number_by_timestamp(reader.meta.time[0], dump_interval)
        end_frame = reader.meta.get_frame_number_by_timestamp(reader.meta.time[1], dump_interval)
    else:
        if time_range[0] > time_range[1]: raise ValueError("End time should not be less, than start time. The values were: {}, {}".format(time_range[0], time_range[1]))
        if time_range[1] > reader.meta.time[1]: raise IndexError("End time is out of simulation range {}. The value was {}".format(reader.meta.end_time, time_range[1]))

        start_frame = reader.meta.get_frame_number_by_timestamp(time_range[0], dump_interval)
        end_frame = reader.meta.get_frame_number_by_timestamp(time_range[1], dump_interval)

    FFMpegWriter = ani.writers['ffmpeg']
    metadata = dict(title='Plasma Multiplot Movie', artist='Matplotlib',
                    comment='Movie support!')
    writer = FFMpegWriter(fps = reader.meta.video_fps,
                          metadata = metadata,
                          codec = reader.meta.video_codec,
                          bitrate = reader.meta.video_bitrate)

    if dry_run: video_file = '/dev/null'
    fig = plot.get_figure()

    with writer.saving(fig, video_file, reader.meta.figure_dpi):
        for i in range(start_frame, end_frame):
            if i % frame_step == 0:
                sys.stdout.write('Loading dataset ' + str(i) + '... ')
                sys.stdout.flush()
                data_1 = reader.rec("density/beam_electrons",
                                    frame_src_size, i)[frame_size[0]:frame_size[2],
                                                       frame_size[1]:frame_size[3]]
                data_2 = reader.rec("density/ions",
                                    frame_src_size, i)[frame_size[0]:frame_size[2],
                                                       frame_size[1]:frame_size[3]]
                data_3 = reader.rec("density/electrons",
                                    frame_src_size, i)[frame_size[0]:frame_size[2],
                                                       frame_size[1]:frame_size[3]]
                data_4 = reader.rec("temperature/ions",
                                    frame_src_size, i)[frame_size[0]:frame_size[2],
                                                       frame_size[1]:frame_size[3]]
                data_5 = reader.rec("temperature/electrons",
                                    frame_src_size, i)[frame_size[0]:frame_size[2],
                                                       frame_size[1]:frame_size[3]]


                # add timestamp to each frame
                timestamp = reader.meta.get_timestamp_by_frame_number(i, dump_interval)
                fig.suptitle("Time: {:.2e} s".format(timestamp), x=.85, y=.95)

                plot.add_image(plot_name_1, data_1, cmap=cmap_r, clim=clim_1)
                plot.add_image(plot_name_2, data_2, cmap=cmap, clim=clim_2)
                plot.add_image(plot_name_3, data_3, cmap=cmap, clim=clim_3)
                plot.add_image(plot_name_4, data_4, cmap=cmap, clim=clim_4)
                plot.add_image(plot_name_5, data_5, cmap=cmap, clim=clim_5)

                if view: plot.redraw()
                if not dry_run: writer.grab_frame()
                print('DONE')
                # input()


def main():
    ## configure RC properties
    # plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

    ####
    parser = argparse.ArgumentParser(description='Tool for creating movie file from PDP3 modelled data.')
    parser.add_argument('data_path', metavar='data_path', type=str,
                        help='Full path to properties.xml')

    parser.add_argument('--video-file', type=str,
                        help='Full path to output video file. Default <path/to/parameters.xml>/field_movie.avi',
                        default=None)

    parser.add_argument('--specie', type=str,
                        help='''Use custom colormap. Default %s.
                        Reference: https://matplotlib.org/examples/color/colormaps_reference.html''' % 'electrons',
                        default='electrons')

    parser.add_argument('--time-range', type=str, help='Time range', default=None)

    parser.add_argument('--cmap', type=str,
                        help='''Use custom colormap. Default %s.
                        Reference: https://matplotlib.org/examples/color/colormaps_reference.html''' % 'gray',
                        default='terrain')

    parser.add_argument('--frame-step', type=int,
                        help='Use only every Nth frame for movie creation. Default: 1',
                        default=1)

    parser.add_argument('--clim', type=str,
                        help='Color limit range for density. Format: "bottom:top"', default=None)

    parser.add_argument('--frame-size', type=str,
                        help='Size of movie images frame. Format: "bottom_r:bottom_z:top_r:top_z"', default=None)

    parser.add_argument('--dry-run', action='store_true', help='Do not write anything. Just for debug', default=False)

    parser.add_argument('--view', action='store_true', default=False,
                        help='View animation as well as write it to file')

    parser.add_argument('--with-grid', action='store_true', help='Use tick grid for plots', default=False)

    args = parser.parse_args()

    time_range=None

    # check if config file exists
    if args.time_range:
        time_range = list(map(float, args.time_range.split(':')))

    clim = list(map(float, args.clim.split(':'))) if args.clim else None
    frame_size = list(map(int, args.frame_size.split(':'))) if args.frame_size else None
    run(args.data_path,
        args.specie,
        clim,
        video_file=args.video_file,
        time_range=time_range,
        cmap=args.cmap,
        frame_step=args.frame_step,
        frame_size=frame_size,
        dry_run=args.dry_run,
        view=args.view,
        use_grid=args.with_grid)
    if args.view:
        input("Press 'Return' to exit ")


## call main function
if __name__ == "__main__":
    main()
