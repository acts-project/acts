# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# detray json schema definitions
from impl import merge_surfaces, update_grids, update_material

from options import (
    common_options,
    parse_common_options,
    detector_io_options,
    parse_detector_io_options,
)

# python includes
import argparse
import json
import os

""" Update the file name for the output json files"""


def __update_filename(input_path):
    _, filename = os.path.split(input_path)
    stem, extension = os.path.splitext(filename)

    return f"{stem}_merged{extension}"


def __main__():
    # ---------------------------------------------------------------arg parsing

    descr = "Detray Surface Merging Service"

    # Define options
    parent_parsers = [
        common_options(descr),
        detector_io_options(),
    ]

    parser = argparse.ArgumentParser(description=descr, parents=parent_parsers)
    args = parser.parse_args()

    logging = parse_common_options(args, descr)
    parse_detector_io_options(args, logging)

    if args.grid_file == "" or args.material_file == "":
        logging.warning(
            "If there are grid or material files, these need to be updated as well!!!\n"
        )

    # -----------------------------------------------------------------------run

    with open(args.geometry_file) as geo_file:
        try:
            geo_json = json.load(geo_file)
        except json.decoder.JSONDecodeError:
            logging.error(f"Invalid json file: {geo_file}")
        else:
            # Merge the portal surfaces
            logging.info(f"Read geometry file {args.geometry_file}")
            out_geo_json, sf_index_dict = merge_surfaces(logging, geo_json)

            out_file_name = __update_filename(args.geometry_file)
            with open(out_file_name, "w") as out_file:
                json.dump(out_geo_json, out_file, indent=4)

            logging.info(f"Geometry: Finished, written to: {out_file_name}\n")

            # Update the grids
            if args.grid_file != "":
                with open(args.grid_file) as grid_file:
                    try:
                        grid_json = json.load(grid_file)
                    except json.decoder.JSONDecodeError:
                        logging.error(f"Invalid json file: {grid_file}")
                    else:
                        logging.info(f"Read grid file {args.grid_file}")
                        out_grid_json = update_grids(logging, grid_json, sf_index_dict)

                        out_file_name = __update_filename(args.grid_file)
                        with open(out_file_name, "w") as out_file:
                            json.dump(out_grid_json, out_file, indent=4)

                        logging.info(f"Grids: Finished, written to: {out_file_name}\n")

            # Update the material
            if args.material_file != "":
                with open(args.material_file) as mat_file:
                    try:
                        mat_json = json.load(mat_file)
                    except json.decoder.JSONDecodeError:
                        logging.error(f"Invalid json file: {mat_file}")
                    else:
                        logging.info(f"Read material file {args.material_file}")
                        out_mat_json = update_material(logging, mat_json, sf_index_dict)

                        out_file_name = __update_filename(args.material_file)
                        with open(out_file_name, "w") as out_file:
                            json.dump(out_mat_json, out_file, indent=4)

                        logging.info(
                            f"Material: Finished, written to: {out_file_name}\n"
                        )


# ------------------------------------------------------------------------------

if __name__ == "__main__":
    __main__()

# ------------------------------------------------------------------------------
