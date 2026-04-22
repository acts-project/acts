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
#include "detray/tracks/ray.hpp"

// Detray plugin include(s)
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/plugins/svgtools/writer.hpp"

// Detray test include(s)
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/common/track_generators.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/test/validation/detector_scanner.hpp"
#include "detray/test/validation/svg_display.hpp"
#include "detray/tracks/tracks.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include <actsvg/core.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s)
#include <array>
#include <string>

GTEST_TEST(svgtools, intersections) {
  // This test creates the visualization using the illustrator class.
  // However, for full control over the process, it is also possible to use
  // the tools in svgstools::conversion, svgstools::display, and
  // actsvg::display by converting the object to a proto object, optionally
  // styling it, and then displaying it.

  // Axes.
  const auto axes =
      actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                             actsvg::style::stroke(), "axis1", "axis2");

  // Creating the view.
  const actsvg::views::z_r view;

  // Creating the detector and geomentry context.
  vecmem::host_memory_resource host_mr;
  const auto [det, names] =
      detray::build_toy_detector<detray::test::algebra>(host_mr);

  using detector_t = decltype(det);
  using test_algebra = typename detector_t::algebra_type;

  detector_t::geometry_context gctx{};

  // Creating the illustrator.
  const detray::svgtools::illustrator il{det, names};

  // Drawing the detector.
  const auto svg_det = il.draw_detector(view);

  // Creating the rays.
  using generator_t =
      detray::uniform_track_generator<detray::detail::ray<test_algebra>>;
  auto trk_gen_cfg = generator_t::configuration{};
  trk_gen_cfg.origin(0.f, 0.f, 100.f).phi_steps(10u).theta_steps(10u);

  std::size_t index = 0;
  // Iterate through uniformly distributed momentum directions with ray
  for (const auto test_ray : generator_t{trk_gen_cfg}) {
    // Record all intersections and objects along the ray
    const auto intersection_record =
        detray::detector_scanner::run<detray::ray_scan>(gctx, det, test_ray);

    const std::string name =
        "test_svgtools_intersection_record" + std::to_string(index);

    // Drawing the intersections.
    auto intersections =
        detray::detail::transcribe_intersections(intersection_record);
    const auto svg_ir =
        il.draw_intersections(name, intersections, test_ray.dir(), view);

    detray::svgtools::write_svg(name, {axes, svg_det, svg_ir});

    index++;
  }
}
