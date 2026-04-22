# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# detray json schema definitions
from json_schema import geometry_schema
from json_schema import homogeneous_material_schema
from json_schema import material_map_schema
from json_schema import surface_grid_schema

# python includes
import argparse
import json
import os
import sys
from jsonschema import validate


def __main__():
    # ---------------------------------------------------------------arg parsing

    parser = argparse.ArgumentParser(description="Detray File Validation")

    parser.add_argument(
        "--geometry_file", help=("Input geometry json file."), default="", type=str
    )
    parser.add_argument(
        "--homogeneous_material_file",
        help=("Input homogeneous material json file."),
        default="",
        type=str,
    )
    parser.add_argument(
        "--material_map_file",
        help=("Input material map json file."),
        default="",
        type=str,
    )
    parser.add_argument(
        "--grid_file", help=("Surface grid json file."), default="", type=str
    )

    args = parser.parse_args()

    # Check input json files
    filename_dict = {}

    geo_file = args.geometry_file
    if geo_file != "":
        if not os.path.isfile(geo_file):
            print(f"Geometry file does not exist! ({geo_file})")
            sys.exit(1)
        else:
            filename_dict[geo_file] = geometry_schema

    hom_mat_file = args.homogeneous_material_file
    if hom_mat_file != "":
        if not os.path.isfile(hom_mat_file):
            print(f"Homogeneous material file does not exist! ({hom_mat_file})")
            sys.exit(1)
        else:
            filename_dict[hom_mat_file] = homogeneous_material_schema

    mat_map_file = args.material_map_file
    if mat_map_file != "":
        if not os.path.isfile(mat_map_file):
            print(f"Material map file does not exist! ({mat_map_file})")
            sys.exit(1)
        else:
            filename_dict[mat_map_file] = material_map_schema

    grid_file = args.grid_file
    if grid_file != "":
        if not os.path.isfile(grid_file):
            print(f"Surface grid file does not exist! ({grid_file})")
            sys.exit(1)
        else:
            filename_dict[grid_file] = surface_grid_schema

    # -----------------------------------------------------------------------run

    for filename, schema in filename_dict.items():
        with open(filename) as file:
            try:
                input_json = json.load(file)
            except json.decoder.JSONDecodeError:
                print(f"Invalid json file: {filename}")
            else:
                validate(instance=input_json, schema=schema)
                print(f"{filename}: OK")


# ------------------------------------------------------------------------------

if __name__ == "__main__":
    __main__()

# ------------------------------------------------------------------------------
