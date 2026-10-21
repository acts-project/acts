// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/plugins/svgtools/utils/groups.hpp"

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

GTEST_TEST(svgtools, groups) {
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

  // Visualisation of a group of surfaces.
  const std::array surface_group_indices{1u, 100u, 10u, 200u};

  const auto [svg_surface_group_xy, mat_group_xy] =
      il.draw_surfaces(surface_group_indices, xy);
  detray::svgtools::write_svg("test_svgtools_surface_group_xy",
                              {axes, svg_surface_group_xy});

  const auto [svg_surface_group_zr, mat_group_zr] =
      il.draw_surfaces(surface_group_indices, zr);
  detray::svgtools::write_svg("test_svgtools_surface_group_zr.svg",
                              {axes, svg_surface_group_zr});

  // Visualisation of a group of volumes.
  const std::array volume_group_indices{3u, 5u};

  const auto [vol_group_xy, gr1] = il.draw_volumes(volume_group_indices, xy);
  detray::svgtools::write_svg("test_svgtools_volume_group_xy",
                              {axes, vol_group_xy});

  const auto [vol_group_zr, gr2] = il.draw_volumes(volume_group_indices, zr);
  detray::svgtools::write_svg("test_svgtools_volume_group_zr",
                              {axes, vol_group_zr});

  // We can also use the svgtools::utils to group actsvg::svg::objects into
  // one.
  auto svg_combined_group = detray::svgtools::utils::group("combined_group");
  svg_combined_group.add_object(svg_surface_group_xy);
  svg_combined_group.add_object(vol_group_xy);

  detray::svgtools::write_svg("test_svgtools_combined_group1.svg",
                              {axes, svg_combined_group});

  // Alternatively, this is equivalent to:
  detray::svgtools::write_svg("test_svgtools_combined_group2.svg",
                              {axes, svg_surface_group_xy, vol_group_zr});

  // NOTE: The all svg object's identification must be unique in the
  // file!
}
