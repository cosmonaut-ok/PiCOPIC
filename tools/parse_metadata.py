#!/usr/bin/env python3

import os,sys
import math
import argparse
import json

from pygments import highlight
from pygments.lexers import JsonLexer
from pygments.formatters import Terminal256Formatter
from pygments.styles import STYLE_MAP

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "lib/python"))

from picopic.h5_reader import H5Reader
from picopic.plain_reader import PlainReader

## for phys-info
def langmur_freq(density):
    # pi = 3.1415
    charge_el = -1.6e-19 # coul
    mass_el = 9.1e-31 # kg
    epsilon_0 = 8.8e-12
    ##
    w_p = math.sqrt(density * math.pow(charge_el, 2) / (mass_el * epsilon_0))
    return w_p

def debye_length(density, temperature):
    charge_el = -1.6e-19 # coul
    mass_el = 9.1e-31 # kg
    epsilon_0 = 8.8e-12
    boltzman_const = 1.38e-23

    l_d = epsilon_0 * boltzman_const * temperature / (4 * math.pi * density * math.pow(charge_el, 2))
    l_d = math.sqrt(l_d)
    return l_d

def eV2kelvin(ev):
    res = ev * 1.16045221e4
    return res

def wake_length(density, beam_velocity):
    w_p = langmur_freq(density)
    lambda_ = 2* math.pi * beam_velocity / w_p
    return lambda_
## for phys-info

def phys_info(properties_path):
    # set reader
    if os.path.isfile(os.path.join(properties_path, "metadata.json")):
        reader = PlainReader(path = properties_path, use_cache=False, verbose=False)
    elif os.path.isfile(os.path.join(properties_path, "data.h5")):
        reader = H5Reader(path = properties_path, use_cache=False, verbose=False)
    else:
        reader = PlainReader(path = properties_path, use_cache=False, verbose=False)

    config = reader.meta

    for i in config.species:
        if i.name.lower() == 'electrons':
            el_density = (i.density[0] + i.density[1]) / 2
            el_temperature = i.temperature

    if not el_density:
        print("Incorrect parameters.xml file (no particles->particle_kind->electrons section)")
    else:
        bunch_duration = config.beams[0].bunch_length / config.beams[0].velocity
        interval_bunches_duration = config.beams[0].bunches_distance / config.beams[0].velocity
        beam_length = config.beams[0].bunch_length * config.beams[0].bunch_amount + config.beams[0].bunches_distance * (config.beams[0].bunch_amount - 1)
        beam_duration = beam_length / config.beams[0].velocity
        w_p = langmur_freq(el_density)
        wake_len = wake_length(el_density, config.beams[0].velocity)
        el_temperature = eV2kelvin(el_temperature)
        debye = debye_length(el_density, el_temperature)
        bunch_part_number = math.pi * math.pow(config.beams[0].bunch_radius, 2) * bunch_duration * config.beams[0].velocity * config.beams[0].bunch_density

        print("General:")
        print("Expected plasma frequency:\t %.4g Hz"%(w_p/(2 * math.pi)))
        print("Expected wake wavelength:\t %.2g m"%(wake_len))
        print("Expected Debye length:\t\t %.4g m"%(debye))
        print("Estimated beam passage time:\t %.4g s"%(config.geometry_size[1] / config.beams[0].velocity + beam_duration))
        print("Number of calculation steps is\t %.4g"%((config.time[1] - config.time[0]) / float(config.time[2])))
        print()

        print("Beam:")
        print("Initial beam velocity:\t\t\t %.4g m/s"%(config.beams[0].velocity))
        if config.beams[0].bunch_amount > 1: print("Single bunch duration:\t\t\t %.4g s"%bunch_duration)
        print("Whole beam duration:\t\t\t %.4g s"%beam_duration)
        if config.beams[0].bunch_amount > 1:
            print("Bunch dimensions:\t\t\t length: %.4g m radius: %.4g m"%(config.beams[0].bunch_length, config.beams[0].bunch_radius))

        print("Beam dimensions:\t\t\t length: %.4g m radius: %.4g m"%(beam_length, config.beams[0].bunch_radius))
        if config.beams[0].bunch_amount > 1: print("Distance between bunch beginings:\t %.4g m"%(config.beams[0].bunches_distance + config.beams[0].bunch_length))
        if config.beams[0].bunch_amount > 1: print("Number of bunches in beam:\t\t %.4g"%config.beams[0].bunch_amount)
        if config.beams[0].bunch_amount > 1: print("Number of particles in single bunch:\t %.2g"%(bunch_part_number))
        print("Number of particles in beam:\t\t %.2g"%(bunch_part_number * config.beams[0].bunch_amount))
        print()























def main():
    parser = argparse.ArgumentParser(description='Parser of output files metadata for PiCOPIC.')
    parser.add_argument('properties_path', metavar='properties_path', type=str,
                        help='Full path to properties.xml')
    parser.add_argument('--subtree', type=str,
                        help='Config`s subtree', default=None)
    parser.add_argument('--color-style', type=str,
                        help='Style of metadata colorization (works only with "--colorize" option)', default='autumn')
    parser.add_argument('--help-colors', action='store_true',
                        help='Output list of available color style names and exits', default=False)
    parser.add_argument('--colorize', action='store_true',
                        help='Make metadata output colorful',
                        default=False)
    parser.add_argument('--bo', action='store_true',
                        help='Display build options of PiCOPIC package, generated the data. Shortcut for `--subtree=[build_options]`',
                        default=False)
    parser.add_argument('--phys-info', action='store_true',
                        help='Display useful physical info about background plasma and particle beams',
                        default=False)

    args = parser.parse_args()

    ## display physical info and exit
    if args.phys_info:
        phys_info(args.properties_path)
        sys.exit(0)

    ## display help about color schemes and exit
    if args.help_colors:
        print(list(STYLE_MAP.keys()))
        sys.exit(0)

    # set reader
    if os.path.isfile(os.path.join(args.properties_path, "metadata.json")):
        reader = PlainReader(path = args.properties_path, use_cache=False, verbose=False)
    elif os.path.isfile(os.path.join(args.properties_path, "data.h5")):
        reader = H5Reader(path = args.properties_path, use_cache=False, verbose=False)
    else:
        reader = PlainReader(path = args.properties_path, use_cache=False, verbose=False)
        # raise EnvironmentError("There is no corresponding data/metadata files in the path " + args.properties_path + ". Can not continue.")

    config = reader.meta.json
    if args.bo: args.subtree = False
    if args.subtree:
        try:
            config = eval('config' + str(args.subtree).replace('[', '["').replace(']', '"]'))
        except KeyError:
            print("Invalid metadata subtree key >>", args.subtree, "<< can not continue")
            sys.exit(1)

    if args.bo:
        try:
            config = config['build_options']
        except KeyError:
            print("ERROR: Build options are absent. Seems, metadata is incomplete (is it configfile?).")
            sys.exit(1)

    json_dumped = json.dumps(config, indent=2, sort_keys=True)
    if args.colorize:
        print(highlight(json_dumped, JsonLexer(), Terminal256Formatter(style=args.color_style)))
    else:
        print(json_dumped)

## call main function
if __name__ == "__main__":
    main()
