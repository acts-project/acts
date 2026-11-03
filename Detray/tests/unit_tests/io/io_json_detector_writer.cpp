// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/annulus2D.hpp"
#include "detray/utils/detector_statistics.hpp"
#include "detray/utils/type_list.hpp"

// Detray IO include(s)
#include "detray/io/backend/geometry_writer.hpp"
#include "detray/io/backend/homogeneous_material_writer.hpp"
#include "detray/io/backend/material_map_writer.hpp"
#include "detray/io/backend/surface_grid_writer.hpp"
#include "detray/io/json/detector_writer.hpp"
#include "detray/io/utils/create_path.hpp"

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
#include <stdexcept>
#include <string>
#include <vector>

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

// Detector writer config
const std::string file_path{"./detray_writer_test"};
detray::io::detector_writer_config writer_cfg =
    detray::io::detector_writer_config{}
        .format(detray::io::format::json)
        .path(file_path)
        .source("detray unit test")
        .write_material(true)
        .write_grids(true)
        .replace_files(true);

/// Write a single detector component to json file
/// @param payload the detector payload containing the component
/// @param mode the file open mode
/// @returns the names of the written files
template <detray::concepts::detector detector_t>
void write_component(const detector_t& det,
                     const typename detector_t::name_map& names,
                     io::detector_payload& payload,
                     const detray::io::detector_writer_config& cfg) {
  // Fill the common information
  payload.header = io::generate_common_header(det, cfg.source(), names);

  EXPECT_EQ(payload.header.version, io::detail::get_detray_version());
  EXPECT_FALSE(payload.header.date.empty());
  EXPECT_EQ(payload.header.metadata,
            detray::types::get_name<typename detector_t::metadata>());
  EXPECT_EQ(payload.header.detector, det.name(names));
  EXPECT_EQ(payload.header.source, "detray unit test");

  io::convert_payload_to_json(payload, cfg);
}

}  // anonymous namespace

/// Test the writing of a telescope detector geometry to json
GTEST_TEST(io, json_telescope_geometry_writer) {
  using detector_t =
      host::detector<telescope_metadata<test_algebra, annulus2D>>;

  // Telescope detector
  vecmem::host_memory_resource host_mr;
  auto [det, names] = build_telescope_detector<test_algebra>(
      host_mr, tel_cfg.positions(positions));

  io::detector_payload payload{};
  io::geometry_writer<detector_t>{}.to_payload(det, names, payload);

  write_component(det, names, payload, writer_cfg);

  // Check that files exist
  EXPECT_TRUE(io::file_exists(file_path + "/telescope_detector_geometry.json"));

  // Check sub-header contents
  ASSERT_TRUE(payload.header.geometry.has_value());
  EXPECT_EQ(payload.header.geometry->n_volumes, det.volumes().size());
  EXPECT_EQ(payload.header.geometry->n_surfaces, det.surfaces().size());
}

/// Test the writing of the toy detector material to json
GTEST_TEST(io, json_telescope_material_writer) {
  using detector_t =
      host::detector<telescope_metadata<test_algebra, annulus2D>>;

  // Telescope detector
  vecmem::host_memory_resource host_mr;

  auto [det, names] = build_telescope_detector<test_algebra>(
      host_mr, tel_cfg.positions(positions));

  io::detector_payload payload{};
  io::homogeneous_material_writer<detector_t>{}.to_payload(det, names, payload);

  write_component(det, names, payload, writer_cfg);

  // Check that files exist
  EXPECT_TRUE(io::file_exists(file_path +
                              "/telescope_detector_homogeneous_material.json"));

  // Check sub-header contents
  EXPECT_TRUE(payload.header.homogeneous_material.has_value());
  EXPECT_EQ(payload.header.homogeneous_material->n_slabs,
            detray::n_material_slabs(det));
  EXPECT_EQ(payload.header.homogeneous_material->n_slab_surfaces,
            detray::n_material_slabs(det));
  EXPECT_EQ(payload.header.homogeneous_material->n_rods,
            detray::n_material_rods(det));
}

/// Test the writing of the toy detector grids to json
GTEST_TEST(io, json_toy_material_maps_writer) {
  using detector_t = host::detector<test::toy_metadata>;

  // Toy detector
  vecmem::host_memory_resource host_mr;
  toy_det_config<scalar> toy_cfg{};
  toy_cfg.use_material_maps(true);
  auto [det, names] = build_toy_detector<test_algebra>(host_mr, toy_cfg);

  io::detector_payload payload{};
  io::material_map_writer<detector_t>{}.to_payload(det, names, payload);

  write_component(det, names, payload, writer_cfg);

  // Check that files exist
  EXPECT_TRUE(io::file_exists(file_path + "/toy_detector_material_maps.json"));

  // Check sub-header contents
  EXPECT_TRUE(payload.header.material_maps.has_value());
  EXPECT_EQ(payload.header.material_maps->n_grids,
            detray::n_material_maps(det));
}

/// Test the writing of the toy detector grids to json
GTEST_TEST(io, json_toy_grid_writer) {
  using detector_t = host::detector<test::toy_metadata>;

  // Toy detector
  vecmem::host_memory_resource host_mr;
  auto [det, names] = build_toy_detector<test_algebra>(host_mr);

  io::detector_payload payload{};
  io::surface_grid_writer<detector_t>{}.to_payload(det, names, payload);

  write_component(det, names, payload, writer_cfg);

  // Check that files exist
  EXPECT_TRUE(io::file_exists(file_path + "/toy_detector_surface_grids.json"));

  // Check sub-header contents
  EXPECT_TRUE(payload.header.surface_grids.has_value());
  EXPECT_EQ(payload.header.surface_grids->n_grids,
            detray::n_surface_grids(det));
}

/// Test the writing of the entire toy detector to json
GTEST_TEST(io, json_toy_detector_writer) {
  // Toy detector
  vecmem::host_memory_resource host_mr;
  toy_det_config<scalar> toy_cfg{};
  toy_cfg.use_material_maps(true);
  const auto [det, names] = build_toy_detector<test_algebra>(host_mr, toy_cfg);

  io::write_detector_json(det, names, writer_cfg);

  // Check that files exist
  EXPECT_TRUE(io::file_exists(file_path + "/toy_detector_geometry.json"));
  EXPECT_TRUE(io::file_exists(file_path + "/toy_detector_material_maps.json"));
  EXPECT_FALSE(
      io::file_exists(file_path + "/toy_detector_homogeneous_material.json"));
  EXPECT_TRUE(io::file_exists(file_path + "/toy_detector_surface_grids.json"));
}

/// Test the writing of the entire wire chamber to json
GTEST_TEST(io, json_wire_chamber_writer) {
  // Wire chamber
  vecmem::host_memory_resource host_mr;
  wire_chamber_config<scalar> wire_cfg{};
  auto [det, names] = build_wire_chamber<test_algebra>(host_mr, wire_cfg);

  io::write_detector_json(det, names, writer_cfg);

  // Check that files exist
  EXPECT_TRUE(io::file_exists(file_path + "/wire_chamber_geometry.json"));
  EXPECT_FALSE(io::file_exists(file_path + "/wire_chamber_material_maps.json"));
  EXPECT_TRUE(
      io::file_exists(file_path + "/wire_chamber_homogeneous_material.json"));
  EXPECT_TRUE(io::file_exists(file_path + "/wire_chamber_surface_grids.json"));
}

/// Test the generation of the the detector payload
GTEST_TEST(io, json_detector_header) {
  using detector_t = host::detector<test::toy_metadata>;

  vecmem::host_memory_resource host_mr;
  auto toy_writer_cfg =
      io::detector_writer_config{}.source("detray unit test").write_grids(true);

  // Toy detector with homogeneous material
  toy_det_config<scalar> toy_cfg{};
  const auto [det, names] = build_toy_detector<test_algebra>(host_mr, toy_cfg);

  // Only fill the payload: No output converters
  io::detector_payload payload{};
  io::write_detector(det, names, toy_writer_cfg, payload);

  // Common header
  const io::header_payload& header = payload.header;
  EXPECT_EQ(header.version, io::detail::get_detray_version());
  EXPECT_FALSE(header.date.empty());
  EXPECT_EQ(header.metadata,
            detray::types::get_name<typename detector_t::metadata>());
  EXPECT_EQ(payload.header.detector, det.name(names));
  EXPECT_EQ(header.source, "detray unit test");

  // Geometry header
  ASSERT_TRUE(header.geometry.has_value());
  EXPECT_EQ(header.geometry->n_volumes, det.volumes().size());
  EXPECT_EQ(header.geometry->n_surfaces, det.surfaces().size());

  // Surface grids header
  ASSERT_TRUE(header.surface_grids.has_value());
  EXPECT_EQ(header.surface_grids->n_grids, detray::n_surface_grids(det));

  // Homogeneous material header
  ASSERT_TRUE(header.homogeneous_material.has_value());
  EXPECT_EQ(header.homogeneous_material->n_slabs,
            detray::n_material_slabs(det));
  EXPECT_EQ(header.homogeneous_material->n_slab_surfaces,
            detray::n_material_slabs(det));
  EXPECT_EQ(header.homogeneous_material->n_rods, detray::n_material_rods(det));

  // Material maps header (no material maps in the detector)
  ASSERT_TRUE(header.material_maps.has_value());
  EXPECT_EQ(header.material_maps->n_grids, detray::n_material_maps(det));

  // No volume grids in the detector payload
  EXPECT_FALSE(header.volume_grids.has_value());

  // Toy detector with material maps
  toy_cfg.use_material_maps(true);
  const auto [map_det, map_names] =
      build_toy_detector<test_algebra>(host_mr, toy_cfg);

  io::detector_payload map_payload{};
  io::write_detector<detector_t>(map_det, map_names, toy_writer_cfg,
                                 map_payload);

  ASSERT_TRUE(map_payload.header.material_maps.has_value());
  EXPECT_EQ(map_payload.header.material_maps->n_grids,
            detray::n_material_maps(map_det));

  // Write the payload to json files
  toy_writer_cfg.replace_files(true);
  io::convert_payload_to_json(map_payload, toy_writer_cfg);

  // Writing a payload without geometry data fails
  io::detector_payload empty_payload{};
  EXPECT_THROW(io::convert_payload_to_json(empty_payload, toy_writer_cfg),
               std::runtime_error);
}
