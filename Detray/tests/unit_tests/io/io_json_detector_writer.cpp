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
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/annulus2D.hpp"

// Detray IO include(s)
#include "detray/io/backend/geometry_writer.hpp"
#include "detray/io/backend/homogeneous_material_writer.hpp"
#include "detray/io/backend/material_map_writer.hpp"
#include "detray/io/backend/surface_grid_writer.hpp"
#include "detray/io/frontend/detector_writer.hpp"
#include "detray/io/json/json_converter.hpp"

// Detray test include(s)
#include "detray/test/common/build_telescope_detector.hpp"
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/common/build_wire_chamber.hpp"
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <ios>

using namespace detray;

namespace {

using test_algebra = test::algebra;
using scalar = test::scalar;

// Use stereo annulus surfaces
constexpr scalar minR{7.2f * unit<scalar>::mm};
constexpr scalar maxR{12.0f * unit<scalar>::mm};
constexpr scalar minPhi{0.74195f};
constexpr scalar maxPhi{1.33970f};
constexpr scalar offset_x{-2.f * unit<scalar>::mm};
constexpr scalar offset_y{-2.f * unit<scalar>::mm};
mask<annulus2D, test_algebra> ann2{0u,     minR,     maxR,     minPhi,
                                   maxPhi, offset_x, offset_y, 0.f};

// Surface positions
std::vector<scalar> positions = {0.f,   50.f,  100.f, 150.f, 200.f, 250.f,
                                 300.f, 350.f, 400.f, 450.f, 500.f};

tel_det_config tel_cfg{ann2};

}  // anonymous namespace

/// Test the writing of a telescope detector geometry to json
GTEST_TEST(io, json_telescope_geometry_writer) {
  using detector_t = detector<telescope_metadata<test_algebra, annulus2D>>;

  // Telescope detector
  vecmem::host_memory_resource host_mr;
  auto [det, names] = build_telescope_detector<test_algebra>(
      host_mr, tel_cfg.positions(positions));

  io::json_converter<detector_t, io::geometry_writer> geo_writer;
  geo_writer.write(det, names);
}

/// Test the writing of the toy detector material to json
GTEST_TEST(io, json_telescope_material_writer) {
  using detector_t = detector<telescope_metadata<test_algebra, annulus2D>>;

  // Telescope detector
  vecmem::host_memory_resource host_mr;

  auto [det, names] = build_telescope_detector<test_algebra>(
      host_mr, tel_cfg.positions(positions));

  io::json_converter<detector_t, io::homogeneous_material_writer> mat_writer;
  mat_writer.write(det, names);
}

/// Test the writing of the toy detector grids to json
GTEST_TEST(io, json_toy_material_maps_writer) {
  using detector_t = detector<test::toy_metadata>;

  // Toy detector
  vecmem::host_memory_resource host_mr;
  toy_det_config<scalar> toy_cfg{};
  toy_cfg.use_material_maps(true);
  auto [det, names] = build_toy_detector<test_algebra>(host_mr, toy_cfg);

  io::json_converter<detector_t, io::material_map_writer> map_writer;
  map_writer.write(det, names,
                   std::ios::out | std::ios::binary | std::ios::trunc);
}

/// Test the writing of the toy detector grids to json
GTEST_TEST(io, json_toy_grid_writer) {
  using detector_t = detector<test::toy_metadata>;

  // Toy detector
  vecmem::host_memory_resource host_mr;
  auto [det, names] = build_toy_detector<test_algebra>(host_mr);

  io::json_converter<detector_t, io::surface_grid_writer> grid_writer;
  grid_writer.write(det, names,
                    std::ios::out | std::ios::binary | std::ios::trunc);
}

/// Test the writing of the entire toy detector to json
GTEST_TEST(io, json_toy_detector_writer) {
  // Toy detector
  vecmem::host_memory_resource host_mr;
  toy_det_config<scalar> toy_cfg{};
  toy_cfg.use_material_maps(true);
  const auto [det, names] = build_toy_detector<test_algebra>(host_mr, toy_cfg);

  auto writer_cfg =
      io::detector_writer_config{}.format(io::format::json).replace_files(true);
  io::write_detector(det, names, writer_cfg);
}

/// Test the writing of the entire wire chamber to json
GTEST_TEST(io, json_wire_chamber_writer) {
  // Wire chamber
  vecmem::host_memory_resource host_mr;
  wire_chamber_config<scalar> wire_cfg{};
  auto [det, names] = build_wire_chamber<test_algebra>(host_mr, wire_cfg);

  auto writer_cfg =
      io::detector_writer_config{}.format(io::format::json).replace_files(true);
  io::write_detector(det, names, writer_cfg);
}
