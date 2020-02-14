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

    args = parser.parse_args()

    if args.help_colors:
        print(list(STYLE_MAP.keys()))
        sys.exit(0)

    # set reader
    if os.path.isfile(os.path.join(args.properties_path, "metadata.json")):
        reader = PlainReader(path = args.properties_path, use_cache=False, verbose=False)
    elif os.path.isfile(os.path.join(args.properties_path, "data.h5")):
        reader = H5Reader(path = args.properties_path, use_cache=False, verbose=False)
    else:
        raise EnvironmentError("There is no corresponding data/metadata files in the path " + config_path + ". Can not continue.")

    config = reader.meta.json
    if args.subtree:
        try:
            config = eval('config' + str(args.subtree).replace('[', '["').replace(']', '"]'))
        except KeyError:
            print("Invalid metadata subtree key >>", args.subtree, "<< can not continue")
            sys.exit(1)

    json_dumped = json.dumps(config, indent=2, sort_keys=True)
    if args.colorize:
        print(highlight(json_dumped, JsonLexer(), Terminal256Formatter(style=args.color_style)))
    else:
        print(json_dumped)

## call main function
if __name__ == "__main__":
    main()
