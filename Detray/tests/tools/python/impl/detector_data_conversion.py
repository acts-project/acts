# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# python includes
import copy
import json
import math
import sys

"""Append a surface to the volume and save the new surface index"""


def __append_surface(volume, surface, extent, new_idx, idx_dict):
    # If the 'mask' entry has not been updated, remove it now
    if "mask" in surface:
        surface["masks"] = [surface["mask"]]
        del surface["mask"]

    # Map the original surface index to the new one after merging
    idx_dict[surface["index_in_coll"]] = (new_idx, extent)
    surface["index_in_coll"] = new_idx

    # Add the updated surface to the volume
    volume["surfaces"].append(surface)

    # Update and return the current surface index
    new_idx += 1
    return new_idx


""" Add a new value to a dictionary of lists. Append if the key exists"""


def __add_to_dict(key, value, dictionary):
    if key in dictionary:
        dictionary[key].append(value)
    else:
        dictionary[key] = [value]


""" Read geometry data from json and merge surfaces """


def merge_surfaces(logging, in_json):

    # Copy identical data
    out_json = {}
    out_json["header"] = in_json["header"]
    logging.info("Geometry: Converted header")

    out_json["data"] = {}
    out_json["data"]["volumes"] = []

    # Map old to new indices per volume
    sf_index_per_volume = {}
    # Total number of surfaces after merging
    n_surfaces = 0

    # Go through each volume
    for volume in in_json["data"]["volumes"]:

        sf_index_per_volume[volume["index"]] = {}
        sf_index_dict = sf_index_per_volume[volume["index"]]

        # Save the input surfaces
        old_surfaces = volume["surfaces"]
        volume["surfaces"] = []

        # Sort the portals by shape for merging
        shape_dict = {}

        # Go through the surfaces
        sf_idx = 0
        for surface in old_surfaces:
            # Do NOT merge sensitive or passive surfaces
            if surface["type"] == 1 or surface["type"] == 2:
                sf_idx = __append_surface(
                    volume, surface, [0, 0], sf_idx, sf_index_dict
                )
                continue

            # Add portals to shape dict
            if "mask" in surface:
                shape_id = surface["mask"]["shape"]
                __add_to_dict(shape_id, surface, shape_dict)
            elif "masks" in surface:
                shape_id = surface["masks"][0]["shape"]
                __add_to_dict(shape_id, surface, shape_dict)

        # Create new surfaces of same shape by merging masks into list

        # Identify the rings that can be merged to a disc, by comparing the
        # z-position (portals do not have rotations or x, y translation)
        ring_dict = {}
        for ring in shape_dict[6]:
            z_pos = ring["transform"]["translation"][2]
            __add_to_dict(z_pos, ring, ring_dict)

        for rings in ring_dict.values():
            # Nothing to merge, just append to output
            if len(rings) == 1:
                boundaries = (
                    rings[0]["mask"]["boundaries"]
                    if "mask" in rings[0]
                    else rings[0]["masks"][0]["boundaries"]
                )
                extent = [boundaries[0], boundaries[1]]
                sf_idx = __append_surface(
                    volume, rings[0], extent, sf_idx, sf_index_dict
                )

                continue

            # The merged disc
            disc = copy.deepcopy(rings[0])

            # Remove old mask entry and replace by masks list
            if "mask" in disc:
                del disc["mask"]
                disc["masks"] = []
            elif "masks" in disc and len(disc["masks"]) > 1:
                logging.error(f"Geometry: surface already converted {disc}")

            for ring in rings:
                disc["masks"].append(
                    ring["mask"] if "mask" in ring else ring["masks"][0]
                )

            # Copy disc into volume
            boundaries = (
                rings[0]["mask"]["boundaries"]
                if "mask" in ring
                else ring["masks"][0]["boundaries"]
            )
            extent = [boundaries[0], boundaries[1]]
            sf_idx = __append_surface(volume, disc, extent, sf_idx, sf_index_dict)

            # Add also the other surfaces to the dict
            for i in range(1, len(rings)):
                boundaries = (
                    rings[i]["mask"]["boundaries"]
                    if "mask" in ring
                    else ring["masks"][0]["boundaries"]
                )
                extent = [boundaries[0], boundaries[1]]
                sf_index_dict[rings[i]["index_in_coll"]] = (sf_idx - 1, extent)

        # Identify the rings that can be merged to a disc, by comparing the z-position (portals do not have rotations or x, y translation)
        cyl_dict = {}
        for cyl in shape_dict[4]:
            rad = (
                cyl["mask"]["boundaries"][0]
                if "mask" in cyl
                else cyl["masks"][0]["boundaries"][0]
            )
            __add_to_dict(rad, cyl, cyl_dict)

        for sub_cyls in cyl_dict.values():
            # Set the cylinder boundaries correctly
            for sub_cyl in sub_cyls:
                # Apply the translation of the surface to the mask z-boundaries
                z_shift = sub_cyl["transform"]["translation"][2]
                if "mask" in sub_cyl:
                    sub_cyl["mask"]["boundaries"][1] += z_shift
                    sub_cyl["mask"]["boundaries"][2] += z_shift
                else:
                    sub_cyl["masks"][0]["boundaries"][1] += z_shift
                    sub_cyl["masks"][0]["boundaries"][2] += z_shift

                sub_cyl["transform"]["translation"] = [0.0, 0.0, 0.0]

            # Nothing to merge, just append to output
            if len(sub_cyls) == 1:
                boundaries = (
                    sub_cyls[0]["mask"]["boundaries"]
                    if "mask" in sub_cyl
                    else sub_cyls[0]["masks"][0]["boundaries"]
                )
                extent = [boundaries[1], boundaries[2]]
                sf_idx = __append_surface(
                    volume, sub_cyls[0], extent, sf_idx, sf_index_dict
                )
                continue

            # The merged disc
            cyl = copy.deepcopy(sub_cyls[0])

            # Remove old mask entry and replace by masks list
            if "mask" in cyl:
                del cyl["mask"]
                cyl["masks"] = []
            elif "masks" in cyl and len(cyl["masks"]) > 1:
                logging.error(f"Geometry: surface already converted {cyl}")

            # Remove translation of the merged cylinder
            cyl["transform"]["translation"] = [0.0, 0.0, 0.0]

            # Add the sub cylinders as masks
            for sub_cyl in sub_cyls:
                cyl["masks"].append(
                    sub_cyl["mask"] if "mask" in sub_cyl else sub_cyl["masks"][0]
                )

            # Extent of the first cylinder
            boundaries = (
                sub_cyls[0]["mask"]["boundaries"]
                if "mask" in sub_cyl
                else sub_cyl["masks"][0]["boundaries"]
            )
            extent = [boundaries[1], boundaries[2]]

            # Copy disc into volume
            sf_idx = __append_surface(volume, cyl, extent, sf_idx, sf_index_dict)

            # Add also the other surfaces to the dict
            for i in range(1, len(sub_cyls)):
                boundaries = (
                    sub_cyls[i]["mask"]["boundaries"]
                    if "mask" in sub_cyl
                    else sub_cyl["masks"][0]["boundaries"]
                )
                extent = [boundaries[1], boundaries[2]]
                sf_index_dict[sub_cyls[i]["index_in_coll"]] = (sf_idx - 1, extent)

        out_json["data"]["volumes"].append(volume)
        n_surfaces += sf_idx

    logging.info(f"Geometry: Converted {len(out_json['data']['volumes'])} volumes")

    # Copy the volume search data structure
    if "volume_grid" in in_json["data"]:
        out_json["data"]["volume_grid"] = in_json["data"]["volume_grid"]

    logging.info("Geometry: Converted volume acceleration structure")

    # Update the metadata
    out_json["header"]["surface_count"] = n_surfaces
    out_json["header"]["volume_count"] = len(out_json["data"]["volumes"])

    return out_json, sf_index_per_volume


""" Update the surface indices in accelerator grids after surface merging"""


def update_grids(logging, in_json, sf_index_per_volume):

    # Copy identical data
    out_json = {}
    out_json["header"] = in_json["header"]
    logging.info("Grids: Converted header")

    out_json["data"] = {}
    out_json["data"]["grids"] = []

    # Go through each volume
    for grids in in_json["data"]["grids"]:
        # Go through all grids in a volume
        for grid in grids["grid_data"]:
            volume_idx = grid["owner_link"]
            sf_idx_map = sf_index_per_volume[volume_idx]

            # Update the surface indices in the grid bins
            for grid_bin in grid["bins"]:
                for i, sf_idx in enumerate(grid_bin["content"]):
                    grid_bin["content"][i] = sf_idx_map[sf_idx][0]

        out_json["data"]["grids"].append(grids)

    # Update metadata (helps debugging)
    out_json["header"]["grid_count"] = len(out_json["data"]["grids"])

    logging.info(f"Grids: Converted {len(out_json['data']['grids'])} grids")

    return out_json


""" Translate the grid bins to global position ranges, observing the translation of the second local coordinate t_y"""


def __to_glob_grid(logging, mmap):
    axes = mmap["axes"]
    step_x = math.fabs(axes[0]["edges"][1] - axes[0]["edges"][0]) / axes[0]["bins"]
    step_y = math.fabs(axes[1]["edges"][1] - axes[1]["edges"][0]) / axes[1]["bins"]

    # Create a grid of the bin corner positions in global (phi, z) coordinates
    glob_grid = {}
    if mmap["grid_link"]["type"] == 3:
        for b_x in range(0, axes[0]["bins"]):
            x = b_x * step_x + axes[0]["edges"][0]
            for b_y in range(0, axes[1]["bins"]):
                y = b_y * step_y + axes[1]["edges"][0]
                bin_data = mmap["bins"][b_x * axes[1]["bins"] + b_y]
                assert bin_data["loc_index"] == [b_x, b_y]

                # Only one entry in material map bin
                glob_grid[(x, y)] = bin_data["content"][0]
    # Create a grid of the bin corner positions in global (r, phi) coordinates
    else:
        for b_y in range(0, axes[1]["bins"]):
            y = b_y * step_y + axes[1]["edges"][0]
            for b_x in range(0, axes[0]["bins"]):
                x = b_x * step_x + axes[0]["edges"][0]
                bin_data = mmap["bins"][b_y * axes[0]["bins"] + b_x]

                assert bin_data["loc_index"] == [b_x, b_y]

                # Only one entry in material map bin
                glob_grid[(x, y)] = bin_data["content"][0]

    return glob_grid


""" Mix material of overlapping bins, as done in the material mapping in ACTS.
    See https://github.com/acts-project/acts/blob/4d3d89dd3c949cb9addd1bd507d42d1b54e58ad9/Core/src/Material/AverageMaterials.cpp"""


def __mix_material(mat1, mat2):
    result = copy.deepcopy(mat1)

    # Calculate the weights according to the thickness
    total_t = mat1["thickness"] + mat2["thickness"]
    w1 = mat1["thickness"] / total_t
    w2 = mat2["thickness"] / total_t

    # Add the material slab thickness
    result["thickness"] = total_t

    # Material parametrization
    # 0: X0
    # 1: L0
    # 2: Ar
    # 3: Z
    # 4: Mass density
    # 5: Molar densitty (unused)
    # 6: @c material_state enum (solid, liquid, gaseous)
    params1 = mat1["material"]["params"]
    params2 = mat2["material"]["params"]

    if params1 != params2:
        # Average X0
        X0_inv = w1 / params1[0] + w2 / params1[0]
        result["material"]["params"][0] = 1 / X0_inv

        # Average L0
        L0_inv = w1 / params1[1] + w2 / params1[1]
        result["material"]["params"][1] = 1 / L0_inv

        molar_density1 = params1[5]
        molar_density2 = params1[5]

        molar_amount1 = molar_density1 * mat1["thickness"]
        molar_amount2 = molar_density2 * mat2["thickness"]
        molar_amount = molar_amount1 + molar_amount2

        molar_weight1 = molar_amount1 / molar_amount
        molar_weight2 = molar_amount2 / molar_amount

        # Average Ar
        result["material"]["params"][2] = (
            molar_weight1 * params1[2] + molar_weight2 * params2[2]
        )

        # Average Z
        z1 = params1[3]
        z2 = params2[3]
        z = 0.0
        if z1 > 0.0 and z2 > 0.0:
            z = math.exp(z1 * math.log(z1) + z2 * math.log(z2))
        else:
            z = w1 * z1 + w2 * z2

        result["material"]["params"][3] = z

        # Average mass density
        mass1 = params1[4] * mat1["thickness"]
        mass2 = params2[4] * mat2["thickness"]
        result["material"]["params"][4] = (mass1 + mass2) / total_t

        # Average molar density
        result["material"]["params"][5] = molar_amount / total_t

    return result


""" Merge material maps of adjacent portals"""


def __merge_material_maps(logging, mmap, new_mmap, extent, major_dir):

    new_glob_mmap = __to_glob_grid(logging, new_mmap)

    # Bin width of the merged grid
    axes = mmap["axes"]
    step_x = math.fabs(axes[0]["edges"][1] - axes[0]["edges"][0]) / axes[0]["bins"]
    step_y = math.fabs(axes[1]["edges"][1] - axes[1]["edges"][0]) / axes[1]["bins"]

    # Bin width of the grid that should be added
    new_axes = mmap["axes"]
    new_step_x = (
        math.fabs(new_axes[0]["edges"][1] - new_axes[0]["edges"][0])
        / new_axes[0]["bins"]
    )
    new_step_y = (
        math.fabs(new_axes[1]["edges"][1] - new_axes[1]["edges"][0])
        / new_axes[1]["bins"]
    )

    major_step = new_step_x if major_dir == 0 else new_step_y

    for pos, bin_data in new_glob_mmap.items():
        # Compensate the lower edge boundary check by also checking
        # the upper bin edge in the major direction (r for disc, z for cylinder)
        if (
            bin_data["thickness"] != 0
            and (
                extent[0] <= pos[major_dir] or extent[0] < (pos[major_dir] + major_step)
            )
            and pos[major_dir] <= extent[1]
        ):
            # Bin index of the new position in the merged grid
            b_x = round((pos[0] - axes[0]["edges"][0]) / step_x + 1) - 1
            b_y = round((pos[1] - axes[1]["edges"][0]) / step_y + 1) - 1

            gbin = 0
            if mmap["grid_link"]["type"] == 3:
                gbin = b_x * axes[1]["bins"] + b_y
            else:
                gbin = b_y * axes[0]["bins"] + b_x

            if 0 <= gbin and gbin < len(mmap["bins"]):
                orig_bin = mmap["bins"][gbin]["content"][0]

                # Bin clash
                if orig_bin["thickness"] != 0:
                    logging.warning(
                        f"Bin clash (sf {new_mmap['owner_link']}, bin {mmap['bins'][gbin]['loc_index']}): Averaging the material"
                    )
                    orig_bin = __mix_material(mat1=orig_bin, mat2=bin_data)
                else:
                    orig_bin["thickness"] = bin_data["thickness"]
                    orig_bin["material"]["params"] = bin_data["material"]["params"]


""" Update the surface owner indices of the material maps after surface merging"""


def update_material(logging, in_json, sf_index_per_volume):

    # Copy identical data
    out_json = {}
    out_json["header"] = in_json["header"]
    logging.info("Material: Converted header")

    if out_json["header"]["common"]["tag"] != "material_maps":
        logging.error("Material: Only material maps conversion is available!")
        return in_json

    out_json["data"] = {}
    out_json["data"]["grids"] = []

    # Go through each volume grid collection
    for grids in in_json["data"]["grids"]:
        volume_idx = grids["volume_link"]
        sf_idx_map = sf_index_per_volume[volume_idx]

        # Make sure that each surface gets only one material map
        surface_grid_map = {}

        # Output grid data collection
        merged_grid_data = []

        # Set a new maximum extent and bin number for each merged grid
        for grid in grids["grid_data"]:

            # Get the surface index after merging
            old_sf_idx = grid["owner_link"]
            new_sf_idx = sf_idx_map[old_sf_idx][0]

            axes = grid["axes"]

            if new_sf_idx not in surface_grid_map:
                step_x = (
                    math.fabs(axes[0]["edges"][1] - axes[0]["edges"][0])
                    / axes[0]["bins"]
                )
                step_y = (
                    math.fabs(axes[1]["edges"][1] - axes[1]["edges"][0])
                    / axes[1]["bins"]
                )

                # Set this as the merged grid
                surface_grid_map[new_sf_idx] = [
                    copy.deepcopy(grid),
                    [step_x, step_y],
                    1,
                ]

                # Set the new surface index as owner link
                surface_grid_map[new_sf_idx][0]["owner_link"] = new_sf_idx
            else:
                # Increase the counter of grids to be merged
                surface_grid_map[new_sf_idx][2] += 1
                merged_grid = surface_grid_map[new_sf_idx][0]
                merged_axes = merged_grid["axes"]

                # Update the axis span and number of bins of the merged grid
                merged_axes[0]["edges"][0] = min(
                    merged_axes[0]["edges"][0], axes[0]["edges"][0]
                )
                merged_axes[0]["edges"][1] = max(
                    merged_axes[0]["edges"][1], axes[0]["edges"][1]
                )
                merged_axes[1]["edges"][0] = min(
                    merged_axes[1]["edges"][0], axes[1]["edges"][0]
                )
                merged_axes[1]["edges"][1] = max(
                    merged_axes[1]["edges"][1], axes[1]["edges"][1]
                )

        # For all merged grids, adapt the number of bins and set all bins to 0
        for sf_idx, grid_tuple in surface_grid_map.items():

            # Only one grid, no merging required
            if grid_tuple[2] == 1:
                continue

            grid = grid_tuple[0]
            steps = grid_tuple[1]
            axes = grid["axes"]

            nbins_x = math.ceil(
                math.fabs(axes[0]["edges"][1] - axes[0]["edges"][0]) / steps[0]
            )
            nbins_y = math.ceil(
                math.fabs(axes[1]["edges"][1] - axes[1]["edges"][0]) / steps[1]
            )

            axes[0]["bins"] = nbins_x
            axes[1]["bins"] = nbins_y

            zero_bin = copy.deepcopy(grid["bins"][0])
            bin_content = zero_bin["content"][0]
            bin_content["surface_idx"] = sf_idx
            bin_content["thickness"] = 0

            # Reset bin collection and fill it with zero bins
            grid["bins"] = []
            if grid["grid_link"]["type"] == 3:
                for b_x in range(0, axes[0]["bins"]):
                    for b_y in range(0, axes[1]["bins"]):
                        zero_bin["loc_index"] = [b_x, b_y]
                        grid["bins"].append(copy.deepcopy(zero_bin))
            else:
                for b_y in range(0, axes[1]["bins"]):
                    for b_x in range(0, axes[0]["bins"]):
                        zero_bin["loc_index"] = [b_x, b_y]
                        grid["bins"].append(copy.deepcopy(zero_bin))

        # Add every grid to the merged grid
        for grid in grids["grid_data"]:

            old_sf_idx = grid["owner_link"]

            new_sf_idx = sf_idx_map[old_sf_idx][0]
            extent = sf_idx_map[old_sf_idx][1]

            # Only one grid, no merging required
            if surface_grid_map[new_sf_idx][2] == 1:
                continue

            merged_grid = surface_grid_map[new_sf_idx][0]

            # Direction index along which grids have to be merged (cyl vs. disc)
            grid_type = merged_grid["grid_link"]["type"]
            assert grid["grid_link"]["type"] == grid_type

            major_index = 1 if grid_type == 3 else 0

            # Merge the grid by transforming to global coordinate and clipping
            # to surface extent
            __merge_material_maps(logging, merged_grid, grid, extent, major_index)

        # Append the merged grids
        for grid_tuple in surface_grid_map.values():
            merged_grid_data.append(grid_tuple[0])

        # Add the new grid collection to volume
        grids["grid_data"] = merged_grid_data

        assert len(grids["grid_data"]) <= len(surface_grid_map)
        out_json["data"]["grids"].append(grids)

    # Update metadata (helps debugging)
    out_json["header"]["grid_count"] = len(out_json["data"]["grids"])

    logging.info(f"Material: Converted {len(out_json['data']['grids'])} grids")

    return out_json
