// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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

GTEST_TEST(svgtools, material) {
  // This test creates the visualization using the illustrator class.

  // Axes.
  const auto axes = actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                                           actsvg::style::stroke());

  // Creating the views.
  const actsvg::views::x_y xy;
  const actsvg::views::z_r zr;
  const actsvg::views::z_phi zphi;

  // Creating the detector and geomentry context.
  vecmem::host_memory_resource host_mr;
  detray::toy_det_config<detray::test::scalar> toy_cfg{};
  toy_cfg.use_material_maps(true).cyl_map_bins(20, 20).disc_map_bins(5, 20);
  const auto [det, names] =
      detray::build_toy_detector<detray::test::algebra>(host_mr, toy_cfg);

  // Creating the svg generator for the detector.
  detray::svgtools::illustrator il{det, names};
  il.hide_material(false);

  // Indexes of the surfaces in the detector to be visualized.
  std::array indices{0u, 346u, 347u, 578u, 579u, 2992u, 2993u};

  auto& portal_mat_style = il.style()
                               ._detector_style._volume_style._portal_style
                               ._surface_style._material_style;

  portal_mat_style._gradient_color_scale =
      detray::svgtools::styling::colors::gradient::viridis_scale;
  for (detray::dindex i : indices) {
    std::string name = "test_svgtools_material_" + std::to_string(i);
    // Visualization of material map of portal i:
    const auto svg_xy = il.draw_surface_material(i, xy);
    detray::svgtools::write_svg(name + "_xy", {axes, svg_xy});
    const auto svg_zphi = il.draw_surface_material(i, zphi);
    detray::svgtools::write_svg(name + "_zphi", {axes, svg_zphi});
  }

  // Draw multiple material maps together
  portal_mat_style._gradient_color_scale =
      detray::svgtools::styling::colors::gradient::plasma_scale;

  std::vector indices2{578u, 579u, 1034u, 1035u, 1770u, 1771u, 2870u, 2871u};
  std::string name = "test_svgtools_cyl_materials";

  const auto svg_xy = il.draw_surface_materials(indices2, xy);
  detray::svgtools::write_svg(name + "_xy", {axes, svg_xy});
}
