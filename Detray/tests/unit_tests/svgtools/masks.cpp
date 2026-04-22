// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"

// Detray plugin include(s)
#include "detray/plugins/svgtools/conversion/surface.hpp"
#include "detray/plugins/svgtools/writer.hpp"

// Detray test include(s)
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Actsvg include(s)
#include <actsvg/core.hpp>
#include <actsvg/meta.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s)
#include <array>
#include <string>

GTEST_TEST(svgtools, masks) {
  // This tests demonstrate the different masks that can be visualized.

  // Create the axes.
  const auto axes = actsvg::draw::x_y_axes("axes", {-250, 250}, {-250, 250},
                                           actsvg::style::stroke());

  using toy_detector_t = detray::detector<detray::test::toy_metadata>;
  using algebra_t = typename toy_detector_t::algebra_type;
  using transform3_t = detray::dtransform3D<algebra_t>;

  const detray::dvector3D<algebra_t> tr{50.f, 100.f, 0.f};
  const transform3_t transform(tr);
  const actsvg::views::x_y view{};

  // Visualize a 2D annulus.
  // e_min_r, e_max_r, e_min_phi_rel, e_max_phi_rel, e_average_phi, e_shift_x,
  // e_shift_y
  detray::mask<detray::annulus2D, algebra_t> ann2D{0u,   100.f, 200.f, -0.5f,
                                                   0.5f, 0.f,   4.f,   30.f};
  const auto ann2D_proto =
      detray::svgtools::conversion::surface(transform, ann2D);
  const auto ann2D_svg = actsvg::display::surface("", ann2D_proto, view);
  detray::svgtools::write_svg("test_svgtools_annulus2D", {axes, ann2D_svg});

  // Visualize a 2D cylinder.
  // e_r, e_lower_z, e_upper_z
  detray::mask<detray::cylinder2D, algebra_t> cyl2D{0u, 100.f, -10.f, 10.f};
  const auto cyl2D_proto =
      detray::svgtools::conversion::surface(transform, cyl2D);
  const auto cyl2D_svg = actsvg::display::surface("", cyl2D_proto, view);
  detray::svgtools::write_svg("test_svgtools_cylinder2D", {axes, cyl2D_svg});

  // Visualize a 2D rectangle.
  // e_half_x, e_half_y
  detray::mask<detray::rectangle2D, algebra_t> rec2D{0u, 100.f, 100.f};
  const auto rec2D_proto =
      detray::svgtools::conversion::surface(transform, rec2D);
  const auto rec2D_svg = actsvg::display::surface("", rec2D_proto, view);
  detray::svgtools::write_svg("test_svgtools_rectangle2D", {axes, rec2D_svg});

  // Visualize a 2D ring.
  // e_inner_r, e_outer_r
  detray::mask<detray::ring2D, algebra_t> rin2D{0u, 50.f, 100.f};
  const auto rin2D_proto =
      detray::svgtools::conversion::surface(transform, rin2D);
  const auto rin2D_svg = actsvg::display::surface("", rin2D_proto, view);
  detray::svgtools::write_svg("test_svgtools_ring2D", {axes, rin2D_svg});

  // Visualize a 2D trapezoid.
  // e_half_length_0, e_half_length_1, e_half_length_2, e_divisor
  detray::mask<detray::trapezoid2D, algebra_t> tra2D{0u, 100.f, 50.f, 200.f,
                                                     1.f / 400.f};
  const auto tra2D_proto =
      detray::svgtools::conversion::surface(transform, tra2D);
  const auto tra2D_svg = actsvg::display::surface("", tra2D_proto, view);
  detray::svgtools::write_svg("test_svgtools_trapezoid2D", {axes, tra2D_svg});

  // Visualize a line.
  // e_cross_section, e_half_z
  detray::mask<detray::line_circular, algebra_t> lin2D{0u, 10.f, 100.f};
  const auto lin2D_proto =
      detray::svgtools::conversion::surface(transform, lin2D);
  const auto lin2D_svg = actsvg::display::surface("", lin2D_proto, view);
  detray::svgtools::write_svg("test_svgtools_line2D", {axes, lin2D_svg});
}
