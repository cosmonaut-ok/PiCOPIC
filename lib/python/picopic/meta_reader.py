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
import math
import json

class Probe:
    def __init__(self, shape, component, schedule, r_start=-1, z_start=-1, r_end=-1, z_end=-1, specie=None):
        if (shape == 'rec' or shape == 'col' or shape == 'row' or shape == 'dot'):
            self.shape = shape
        else:
            raise ValueError("Wrong probe shape " + shape)

        self.component = component
        self.schedule = schedule
        self.specie = specie
        self.size = [r_start, z_start, r_end, z_end]


class SpecieP:
    def __init__(self, name, macro_amount,
                 left_density, right_density, temperature,
                 charge, mass):
        self.name = name
        self.charge = charge
        self.mass = mass
        self.macro_amount = macro_amount
        self.density = [left_density, right_density]
        self.temperature = temperature


class BeamP(SpecieP):
    def __init__(self, name, macro_amount, velocity, charge, mass,
                 start_time, bunch_density, bunch_amount, bunch_length, bunch_radius, bunches_distance):
        super(BeamP, self).__init__(name, macro_amount, 0, 0, 0, charge, mass)
        self.velocity = velocity
        self.start_time = start_time
        self.bunch_density = bunch_density
        self.bunch_amount = bunch_amount
        self.bunch_length = bunch_length
        self.bunch_radius = bunch_radius
        self.bunches_distance = bunches_distance


class MetaReader:
    def __init__(self, config):
        if isinstance(config, str):
            with open(config) as f:
                cfg = json.load(f)
        elif isinstance(config, dict):
            cfg = config

        self.json = cfg
        self.geometry_size = [cfg['geometry']['size']['radius'], cfg['geometry']['size']['longitude']]
        self.geometry_grid = [cfg['geometry']['grid']['radius'], cfg['geometry']['grid']['longitude']]
        self.geometry_domains = [cfg['geometry']['domains_amount']['radius'], cfg['geometry']['domains_amount']['longitude']]
        self.geometry_pml = [cfg['geometry']['pml']['left_wall'], cfg['geometry']['pml']['right_wall'], cfg['geometry']['pml']['outer_wall'], cfg['geometry']['pml']['sigma_1'], cfg['geometry']['pml']['sigma_2']]
        self.time = [cfg['time']['start'], cfg['time']['end'], cfg['time']['step']]

        self.species = []
        for i in cfg['particles']:
            pp = SpecieP(i['name'], 1e6, i['density'],
                         i['density'], i['temperature'],
                         i['charge'], i['mass'])
            self.species.append(pp)

        self.beams = []
        for i in cfg['beams']:
            bm = BeamP(i['name'], 1e6, i['velocity'], i['charge'], i['mass'],
                       i['start_time'], i['bunch']['density'], i['bunch']['amount'], i['bunch']['length'],
                       i['bunch']['radius'], i['bunch']['distance'])
            self.beams.append(bm)

        if isinstance(config, str):
            self.config_path = os.path.normpath(os.path.dirname(config))
        else:
            self.config_path = cfg['base_path']

        # calculate data path
        if str.startswith(cfg['data']['data_root'], '/'):
            self.data_path = cfg['data']['data_root']
        else:
            self.data_path = os.path.normpath(os.path.join(self.config_path, cfg['data']['data_root']))

        self.use_compression = cfg['data']['compression']['use']
        self.compression_level = cfg['data']['compression']['level']

        self.probes = []
        probe_size = [0, 0, -1, -1]

        for i in cfg['data']['probes']:
            try:
                probe_specie = i['specie']
            except:
                probe_specie = None
            if i['shape'] == 'rec':
                probe_size = [i['start']['r'], i['start']['z'], i['end']['r'], i['end']['z']]
            elif i['shape'] == 'row':
                probe_size = [i['r'], 0, -1, -1]
            elif i['shape'] == 'col':
                probe_size = [0, i['z'], -1, -1]
            elif i['shape'] == 'dot':
                probe_size = [i['r'], i['z'], -1, -1]

            pb = Probe(i['shape'], i['component'], i['schedule'],
                       probe_size[0], probe_size[1], probe_size[2], probe_size[3],
                       probe_specie)
            self.probes.append(pb)

        self.plot_cmap = cfg['data']['plot']['subplot']['cmap']
        self.plot_ticks_x = cfg['data']['plot']['subplot']['ticks']['x_amount']
        self.plot_ticks_y = cfg['data']['plot']['subplot']['ticks']['y_amount']
        self.plot_ticks_z = cfg['data']['plot']['subplot']['ticks']['z_amount']
        self.plot_cbar = cfg['data']['plot']['subplot']['ticks']['cbar_amount']
        self.figure_color = cfg['data']['plot']['figure']['color']
        self.figure_width = cfg['data']['plot']['figure']['width']
        self.figure_height = cfg['data']['plot']['figure']['height']
        self.figure_dpi = cfg['data']['plot']['figure']['dpi']
        self.figure_font_family = cfg['data']['plot']['figure']['font']['family']
        self.figure_font_size = cfg['data']['plot']['figure']['font']['size']
        self.figure_font_name = cfg['data']['plot']['figure']['font']['name']
        self.video_codec = cfg['data']['plot']['video']['codec']
        self.video_fps = cfg['data']['plot']['video']['fps']
        self.video_bitrate = cfg['data']['plot']['video']['bitrate']


    def get_frame_number_by_timestamp(self, timestamp, dump_interval):
        number_frames = round(timestamp / self.time[2] / dump_interval)
        return number_frames


    def get_timestamp_by_frame_number(self, frame_number, dump_interval):
        timestamp = frame_number * self.time[2] * dump_interval
        return timestamp


    def get_clim_estimation(self):
        p_factor_array = []
        factor = 1e-5 # very empirical coefitient :)
        counter = 0

        for i in self.species:
            n_l = i.density[0]
            n_r = i.density[1]
            t = i.temperature

            n = n_l if n_l > n_r else n_r

            p_factor_array.append(factor * n * t * math.sqrt(self.time[2]))
            counter = counter+1

        return sum(p_factor_array)


    def get_row_by_radius(self, radius):
        return(int(round(radius * self.geometry_grid[0] / self.geometry_size[1])))


    def get_col_by_longitude(self, longitude):
        return(int(round(longitude * self.geometry_grid[1] / self.geometry_size[0])))


    def get_radius_by_row(self, row):
        return(row * self.geometry_size[0] / self.geometry_grid[0])


    def get_longitude_by_col(self, col):
        return(col * self.geometry_size[1] / self.geometry_grid[1])
