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
        time_range=None, cmap=None, frame_step=1, frame_size=None, dry_run=False, view=False, use_grid=False):

    ##  configuration options
    x_axis_label = r'$\mathit{Z (m)}$'
    y_axis_label = r'$\mathit{R (m)}$'
    cbar_axis_label = r'$\frac{V}{m}$'
    cbar_bunch_density_axis_label = r'$eV$'

    e_r_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Radial\enspace Component}\enspace(E_r)$'
    e_z_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Longitudal\enspace Component}\enspace(E_z)$'
    rho_beam_plot_name = r'$\mathbf{Temperature\enspace T\enspace (eV)}$'

    # calculate/update video file path
    video_file = os.path.join(config_path, 'field_movie.avi') if not video_file else video_file

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

    clim_rho_beam = [0, 2.5] # 0-2.5:0-5e-10, 0-5:5e-10-1e-9, 0-40:1e-9-1.5e-9, 0-100:1.5e-9-2e-9, 0-200:2e-9-2.5e-9, 0-400:2.5e-9-3e-9 # [-(cfg.bunch_density * el_charge * rho_beam_scale), 0]
    clim_estimation = reader.meta.get_clim_estimation()

    if not clim_e_r: clim_e_r = [-clim_estimation, clim_estimation]
    if not clim_e_z: clim_e_z = [-clim_estimation, clim_estimation]
    if not rho_beam_scale: rho_beam_scale = 1

    if not frame_size: frame_size = [0, 0, reader.meta.geometry_grid[0], reader.meta.geometry_grid[1]]
    frame_src_size=[-1, -1, -1, -1]
    
    # detect probe shape
    for probe in reader.meta.probes:
        if (probe.shape == 'rec') and (probe.size[0] == frame_size[0]) and (probe.size[1] == frame_size[1]) and(probe.size[2] == frame_size[2]) and(probe.size[3] == frame_size[3]):
            frame_src_size = probe.size
    # try bigger frames, if autoselect enabled
    if frame_src_size[0] == -1 and frame_src_size[1] == -1 and frame_src_size[2] == -1 and frame_src_size[3] == -1 and autoselect:
        for probe in reader.meta.probes:
            if (probe.shape == 'rec') and (probe.size[0] <= frame_size[0]) and (probe.size[1] <= frame_size[1]) and(probe.size[2] >= frame_size[2]) and(probe.size[3] >= frame_size[3]):
                frame_src_size = probe.size

    r_scale = (frame_size[2] - frame_size[0]) / reader.meta.geometry_grid[0]
    z_scale = (frame_size[3] - frame_size[1]) / reader.meta.geometry_grid[1]
    # define plot builder
    plot = PlotBuilder(frame_size[3] - frame_size[1],
                       frame_size[2] - frame_size[0],
                       fig_color=reader.meta.figure_color,
                       fig_width=reader.meta.figure_width,
                       fig_height=reader.meta.figure_height,
                       fig_dpi=reader.meta.figure_dpi,
                       font_family=reader.meta.figure_font_family,
                       font_name=reader.meta.figure_font_name,
                       font_size=reader.meta.figure_font_size,
                       x_ticklabel_end=reader.meta.geometry_size[1] * z_scale,
                       y_ticklabel_end=reader.meta.geometry_size[0] * r_scale,
                       tickbox=True, grid=use_grid, is_invert_y_axe=False,
                       aspect='equal', image_interpolation=plot_image_interpolation,
                       guess_number_ticks = 20)


    # add subplots
    plot.add_subplot_cartesian_2d(e_r_plot_name, 312, x_axe_label=x_axis_label, y_axe_label=y_axis_label)
    plot.add_subplot_cartesian_2d(e_z_plot_name, 311, x_axe_label=x_axis_label, y_axe_label=y_axis_label)
    plot.add_subplot_cartesian_2d(rho_beam_plot_name, 313, x_axe_label=x_axis_label, y_axe_label=y_axis_label)

    # add initial image with zeros and colorbar
    initial_image = np.zeros([frame_size[2] - frame_size[0], frame_size[3] - frame_size[1]])

    plot.add_image(e_r_plot_name, initial_image, cmap=cmap, clim=clim_e_r)
    plot.add_colorbar(e_r_plot_name, ticks=clim_e_r, title=cbar_axis_label)

    plot.add_image(e_z_plot_name, initial_image, cmap=cmap, clim=clim_e_z)
    plot.add_colorbar(e_z_plot_name, ticks=clim_e_z, title=cbar_axis_label)

    plot.add_image(rho_beam_plot_name, initial_image, cmap="{}_r".format(cmap), clim=clim_rho_beam)
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
        if time_range[1] > reader.meta.time[1]: raise IndexError("End time is out of simulation range {}. The value was {}".format(reader.meta.time[1], time_range[1]))

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
            if i % frame_step == 0:
                sys.stdout.write('Loading dataset ' + str(i) + '... ')
                sys.stdout.flush()
                data_r = reader.get_frame('E_r', frame_src_size, i)[frame_size[0]:frame_size[2], frame_size[1]:frame_size[3]]
                data_z = reader.get_frame('E_z', frame_src_size, i)[frame_size[0]:frame_size[2], frame_size[1]:frame_size[3]]
                data_beam = reader.get_frame('temperature/electrons', frame_src_size, i)[frame_size[0]:frame_size[2], frame_size[1]:frame_size[3]]

                # add timestamp to each frame
                timestamp = reader.meta.get_timestamp_by_frame_number(i, dump_interval)
                fig.suptitle("Time: {:.2e} s".format(timestamp), x=.85, y=.95)
                plot.add_image(e_r_plot_name, data_r, cmap=cmap, clim=clim_e_r)
                plot.add_image(e_z_plot_name, data_z, cmap=cmap, clim=clim_e_z)
                plot.add_image(rho_beam_plot_name, data_beam, cmap="{}_r".format(cmap), clim=clim_rho_beam)

                if view: plot.redraw()
                if not dry_run: writer.grab_frame()
                print('DONE')


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

    parser.add_argument('--time-range', type=str, help='Time range', default=None)

    parser.add_argument('--cmap', type=str,
                        help='''Use custom colormap. Default %s.
                        Reference: https://matplotlib.org/examples/color/colormaps_reference.html''' % 'gray',
                        default=None)

    parser.add_argument('--beam-scale-factor', type=float,
                        help='''Beam density setting automatically, but you can set scale factor to sets,
                        where initial bunch density should be placed in color range''',
                        default=1)

    parser.add_argument('--frame-size', type=str,
                        help='Size of movie images frame. Format: "bottom_r:bottom_z:top_r:top_z"', default=None)

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
    frame_size = list(map(int, args.frame_size.split(':'))) if args.frame_size else None
    run(args.data_path,
        clim_e_r=clim_e_r,
        clim_e_z=clim_e_z,
        rho_beam_scale=args.beam_scale_factor,
        video_file=args.video_file,
        time_range=time_range,
        frame_size=frame_size,
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
