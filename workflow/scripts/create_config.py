#!/usr/bin/env python3
"""to create config file in yaml format"""

import yaml
import argparse


def get_args():
    """
    Get command line arguments
    """

    parser = argparse.ArgumentParser(
        description="Script for creating config file in YAML format.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("strains_file", metavar="<strains file>",
                        help="file with strain names, one name per line")
    parser.add_argument("output_file", metavar="<config.yaml>",
                        help="config file name")

    return parser.parse_args()


def main(strains_file):
    """create config dictionary"""
    # read file
    with open(strains_file, 'r') as f:
        strains = [line.rstrip() for line in f.readlines()]

    # create dict
    config = dict({"strains": {}})

    # fill in the dict
    for strain in strains:
        config["strains"][strain] = strain

    return config


if __name__ == '__main__':
    args = get_args()
    input_strains_file = args.strains_file
    output_filename = args.output_file

    config_dict = main(strains_file=input_strains_file)

    # write as yaml
    with open(output_filename, 'w') as outfile:
        yaml.dump(config_dict, outfile, default_flow_style=False)
