# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# detray includes
from options import (
    common_options,
    detector_io_options,
    display_options,
)
from options import (
    parse_common_options,
    parse_detector_io_options,
    parse_display_options,
    fill_reader_config,
)

# detray python bindings
import detray.core
import detray.io
import detray.svgtools

# python includes
import argparse
import os


def run_detector_display(args, outdir):

    # Read the detector, but also display incorrect geometries for debugging
    reader_cfg = fill_reader_config(args, detray.io.DetectorReaderConfig())
    reader_cfg.doCheck = False
    det, names = detray.io.readDetector(detray.core.HostMemoryResource(), reader_cfg)

    # Style settings for the illustrator
    style = detray.svgtools.tableauColorblindStyle()
    style.fontSize = args.font_size
    # Align the material gradient box with the major axis
    style.materialGradientPos = [args.material_gradient_axis + 100.0, -200.0]

    # Creating the svg generator for the detector
    il = detray.svgtools.Illustrator(det, names, style)
    il.showInfo(args.show_info)
    il.hideEtaLines(args.hide_eta_lines)
    il.hidePortals(args.hide_portals)
    il.hidePassives(args.hide_passives)
    il.hideMaterial(not args.material_file or args.hide_material)
    il.hideGrids(not args.grid_file)
    il.searchWindow(args.search_window)

    # The geometry context to be displayed
    gctx = detray.core.GeometryContext(args.context)

    # x-y, z-r and z-phi axes
    xy_axis = detray.svgtools.drawAxes(
        "axes",
        [-args.x_axis, args.x_axis],
        [-args.y_axis, args.y_axis],
        "x",
        "y",
        args.font_size,
    )
    zr_axis = detray.svgtools.drawAxes(
        "axes",
        [-args.z_axis, args.z_axis],
        [-5.0, args.r_axis],
        "z",
        "r",
        args.font_size,
    )
    zphi_axis = detray.svgtools.drawAxes(
        "axes",
        [-args.z_axis, args.z_axis],
        [-args.r_axis, args.r_axis],
        "z",
        "phi",
        args.font_size,
    )

    # Creating the views
    xy = detray.svgtools.ViewXY()
    zr = detray.svgtools.ViewZR()
    zphi = detray.svgtools.ViewZPhi()
    zrphi = detray.svgtools.ViewZRPhi()

    def write_svg(name, svgs):
        detray.svgtools.writeSvg(os.path.join(outdir, name), svgs)

    # Draw the volumes from a collection of volume identifiers
    def draw_volumes(vol_ids):
        vol_xy_svg, xy_grids = il.drawVolumes(vol_ids, xy, gctx)
        write_svg(vol_xy_svg.id, [xy_axis, vol_xy_svg])
        for grid in xy_grids:
            write_svg(grid.id, [grid])

        vol_zr_svg, _ = il.drawVolumes(vol_ids, zr, gctx)
        write_svg(vol_zr_svg.id, [zr_axis, vol_zr_svg])

        _, zrphi_grids = il.drawVolumes(vol_ids, zrphi, gctx)
        for grid in zrphi_grids:
            write_svg(grid.id, [grid])

    # Display the volumes
    if args.volume_indices:
        draw_volumes(args.volume_indices)
    if args.volume_names:
        draw_volumes(args.volume_names)

    # Display the surfaces
    if args.surfaces:
        sf_xy_svg, mat_xy_svg = il.drawSurfaces(args.surfaces, xy, gctx)
        write_svg(sf_xy_svg.id, [xy_axis, sf_xy_svg])
        write_svg(mat_xy_svg.id, [xy_axis, mat_xy_svg])

        sf_zr_svg, _ = il.drawSurfaces(args.surfaces, zr, gctx)
        write_svg(sf_zr_svg.id, [zr_axis, sf_zr_svg])

        _, mat_zphi_svg = il.drawSurfaces(args.surfaces, zphi, gctx)
        write_svg(mat_zphi_svg.id, [zphi_axis, mat_zphi_svg])

    # If nothing was specified, display the whole detector
    if not args.volume_indices and not args.volume_names and not args.surfaces:
        det_xy_svg = il.drawDetector(xy, gctx)
        write_svg(det_xy_svg.id, [xy_axis, det_xy_svg])

        det_zr_svg = il.drawDetector(zr, gctx, args.r_axis, args.z_axis)
        write_svg(det_zr_svg.id, [zr_axis, det_zr_svg])

    # Display the detector volume graph
    if args.write_volume_graph:
        graph = detray.core.VolumeGraph(det)

        graph_file = os.path.join(outdir, il.detName + "_volume_graph.dot")
        with open(graph_file, "w") as out_file:
            out_file.write(graph.toDotString())


def __main__():

    # ---------------------------------------------------------------arg parsing

    descr = "Detray Detector Display"

    # Define options
    parent_parsers = [
        common_options(descr),
        detector_io_options(),
        display_options(),
    ]

    parser = argparse.ArgumentParser(description=descr, parents=parent_parsers)

    args = parser.parse_args()

    logging = parse_common_options(args, descr)
    parse_detector_io_options(args, logging)
    outdir = parse_display_options(args, logging)

    # -----------------------------------------------------------------------run

    logging.debug("Drawing the detector")
    run_detector_display(args, outdir)

    logging.info(f"Wrote the svg files to '{outdir}'")


# ------------------------------------------------------------------------------

if __name__ == "__main__":
    __main__()

# ------------------------------------------------------------------------------
