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
#include "detray/tracks/trajectories.hpp"

// Detray plugin include(s)
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/plugins/svgtools/utils/groups.hpp"
#include "detray/plugins/svgtools/writer.hpp"

// Detray test include(s)
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/test/validation/detector_scanner.hpp"
#include "detray/test/validation/svg_display.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include <actsvg/core.hpp>
#include <actsvg/web/web_builder.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s)
#include <array>
#include <filesystem>
#include <string>

GTEST_TEST(svgtools, web) {
  // In this test we will create a web page to show the detector geometry and
  // more. We will start by creating the svgs we want to include in the web
  // page.

  // Axes.
  const auto axes = actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                                           actsvg::style::stroke());

  // Creating the views.
  const actsvg::views::x_y view;

  // Creating the detector and geomentry context.
  vecmem::host_memory_resource host_mr;
  const auto [det, names] =
      detray::build_toy_detector<detray::test::algebra>(host_mr);
  using detector_t = decltype(det);

  using test_algebra = typename detector_t::algebra_type;
  using scalar = detray::dscalar<test_algebra>;
  using vector3 = detray::dvector3D<test_algebra>;

  detector_t::geometry_context gctx{};

  // Creating the svg generator for the detector.
  const detray::svgtools::illustrator il{det, names};

  // The vector of svgs that we want to include on the webpage.
  std::vector<actsvg::svg::object> svgs;

  // Indexes of the volumes in the detector to be visualized.
  std::array indices{0u,  1u,  2u,  3u,  4u,  5u,  6u,  7u,  8u,  9u,
                     10u, 11u, 12u, 13u, 14u, 15u, 16u, 17u, 18u, 19u};

  // Draw the volumes and include them in the svg vector.
  for (detray::dindex i : indices) {
    const auto [svg, _] = il.draw_volume(i, view);
    svgs.push_back(svg);
  }

  // Draw some example trajectories and include them in the svg vector (along
  // with their intersections).
  for (const auto qop : std::vector{-4, -8, -16}) {
    std::string name = "Helix_qop_" + std::to_string(qop) + ")";

    const typename detector_t::point3_type ori{0.f, 0.f, 80.f};
    const typename detector_t::point3_type dir{0.f, 1.f, 1.f};

    // Create the helix trajectory.
    // Constant magnetic field
    vector3 B{0.f * detray::unit<scalar>::T, 0.f * detray::unit<scalar>::T,
              1.f * detray::unit<scalar>::T};

    const detray::detail::helix<test_algebra> helix(
        ori, 0.f, detray::vector::normalize(dir), static_cast<float>(qop), B);
    const auto helix_ir =
        detray::detector_scanner::run<detray::helix_scan>(gctx, det, helix);

    // Draw the helix trajectory.
    const auto svg_helix =
        il.draw_trajectory(name + "_trajectory", helix, view);

    // Draw the intersection record.
    auto helix_intersections =
        detray::detail::transcribe_intersections(helix_ir);
    const auto svg_helix_ir = il.draw_intersections(
        name + "_record", helix_intersections, helix.dir(), view);

    // We one the trajectory and intersection record to be considered as one
    // svg. Thus we group them together before adding the group to the svg
    // vector.
    const auto svg_group = detray::svgtools::utils::group(
        name, std::vector{svg_helix, svg_helix_ir});
    svgs.push_back(svg_group);
  }

  // The output directory for the web page.
  const auto current_directory = std::filesystem::current_path();

  // Create the web page builder.
  actsvg::web::web_builder builder{};

  // To visualize the svg objects in a specific order we need to pass a
  // comparator before we build the page. For instance we might want to
  // display helices on top on the volumes (rather than the other way around).
  // In this example we will alphanumeric_compare which renders the svgs such
  // that the id of the svg object with the greatest id alphanumerically will
  // be displayed on top. Build the web page.
  auto alphanum_cmp = actsvg::web::compare::alphanumeric;
  builder.build(current_directory / "test_svgtools_website", svgs,
                alphanum_cmp);

  // Once the direcroy has been created, run the server using "python3 -m
  // http.server" in the directory. Subsequently connect to localhost using
  // the respective port. On the web page, click and drag to move the view
  // box. Use the scroll wheel to zoom.
}
