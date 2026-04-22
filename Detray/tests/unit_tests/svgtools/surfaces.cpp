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

GTEST_TEST(svgtools, surfaces) {
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
  const detray::svgtools::illustrator il{det, names};

  // Indexes of the surfaces in the detector to be visualized.
  std::array indices{200u, 201u, 202u, 203u, 204u};

  for (detray::dindex i : indices) {
    std::string name = "test_svgtools_surface" + std::to_string(i);
    // Visualization of surface/portal i:
    const auto [svg_xy, mat_xy] = il.draw_surface(i, xy);
    detray::svgtools::write_svg(name + "_xy", {axes, svg_xy});
    detray::svgtools::write_svg(name + "mat_xy", {axes, mat_xy});
    const auto [svg_zr, mat_zr] = il.draw_surface(i, zr);
    detray::svgtools::write_svg(name + "_zr", {axes, svg_zr});
    detray::svgtools::write_svg(name + "mat_zr", {axes, mat_zr});
  }
}
