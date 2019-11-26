#!/usr/bin/env python3

import os,sys
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "lib/python"))

import math

# import sys
# import os

import argparse

# from numpy import *
# from pylab import *

# import matplotlib.pyplot as plt
# import matplotlib.animation as ani

from picopic.cfg import Cfg

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

def main():
    parser = argparse.ArgumentParser(description='Quick calculator for some PDP3 model parameters.')
    parser.add_argument('properties_path', metavar='properties_path', type=str,
                        help='Full path to properties.xml')

    args = parser.parse_args()

    config = Cfg(args.properties_path)

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

## call main function
if __name__ == "__main__":
    main()
