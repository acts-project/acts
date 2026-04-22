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

// GTest include(s).
#include <gtest/gtest.h>

GTEST_TEST(svgtools, grids) {
  // This test creates the visualization using the illustrator class.
  // However, for full control over the process, it is also possible to use
  // the tools in svgstools::conversion, svgstools::display, and
  // actsvg::display by converting the object to a proto object, optionally
  // styling it, and then displaying it.

  // Creating the detector and geomentry context.
  vecmem::host_memory_resource host_mr;
  const auto [det, names] =
      detray::build_toy_detector<detray::test::algebra>(host_mr);

  // Creating the view.
  const actsvg::views::x_y view;

  // Creating the svg generator for the detector.
  detray::svgtools::illustrator il{det, names};

  // In this example we want to draw the grids of the volumes with indices 0,
  // 1, ... in the detector.
  std::vector<detray::dindex> indices = {3u, 5u, 7u, 9u};

  for (const detray::dindex i : indices) {
    // Draw volume i.
    il.hide_grids(false);
    const auto [volume_svg, grid] = il.draw_volume(i, view);

    // Write volume i and its grid
    detray::svgtools::write_svg("test_svgtools_" + volume_svg._id, volume_svg);
    detray::svgtools::write_svg("test_svgtools_" + grid._id, grid);
  }
}
