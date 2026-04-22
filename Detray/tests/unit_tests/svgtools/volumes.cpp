// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"

// Detray plugin include(s)
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/plugins/svgtools/writer.hpp"

// Detray test include(s)
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include <actsvg/core.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s)
#include <array>
#include <string>

GTEST_TEST(svgtools, volumes) {
  // This test creates the visualization using the illustrator class.
  // However, for full control over the process, it is also possible to use
  // the tools in svgstools::conversion, svgstools::display, and
  // actsvg::display by converting the object to a proto object, optionally
  // styling it, and then displaying it.

  // Axes.
  const auto axes = actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                                           actsvg::style::stroke());

  // Creating the views.
  const actsvg::views::x_y xy;
  const actsvg::views::z_r zr;

  // Creating the detector and geomentry context.
  vecmem::host_memory_resource host_mr;
  const auto [det, names] =
      detray::build_toy_detector<detray::test::algebra>(host_mr);

  // Creating the svg generator for the detector.
  detray::svgtools::illustrator il{det, names};
  // Only test volume conversion here
  il.hide_grids(true);

  // Indexes of the volumes in the detector to be visualized.
  std::array indices{0u,  1u,  2u,  3u,  4u,  5u,  6u,  7u,  8u,  9u,
                     10u, 11u, 12u, 13u, 14u, 15u, 16u, 17u, 18u, 19u};

  for (detray::dindex i : indices) {
    // Visualization of volume i:
    const auto [vol_svg_xy, xy_grid] = il.draw_volume(i, xy);
    detray::svgtools::write_svg("test_svgtools_" + vol_svg_xy._id,
                                {axes, vol_svg_xy});

    const auto [vol_svg_zr, zr_grid] = il.draw_volume(i, zr);
    detray::svgtools::write_svg("test_svgtools_" + vol_svg_zr._id,
                                {axes, vol_svg_zr});
  }

  // Names of volumes to be visualized
  std::array vol_names{"beampipe_0", "gap_14", "endcap_1", "connector_gap_2",
                       "barrel_7"};

  for (std::string_view name : vol_names) {
    // Visualization of volume i:
    const auto [vol_svg_xy, xy_grid] = il.draw_volume(name, xy);
    detray::svgtools::write_svg(
        "test_svgtools_" + vol_svg_xy._id + "_fetched_by_name",
        {axes, vol_svg_xy});

    const auto [vol_svg_zr, zr_grid] = il.draw_volume(name, zr);
    detray::svgtools::write_svg(
        "test_svgtools_" + vol_svg_zr._id + "_fetched_by_name",
        {axes, vol_svg_zr});
  }
}
