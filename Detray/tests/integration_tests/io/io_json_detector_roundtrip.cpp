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
#include "detray/utils/consistency_checker.hpp"
#include "detray/utils/logging.hpp"

// Detray IO include(s)
#include "detray/io/backend/geometry_reader.hpp"
#include "detray/io/backend/geometry_writer.hpp"
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/io/frontend/detector_writer.hpp"
#include "detray/io/json/json_converter.hpp"

// Detray test include(s)
#include "detray/test/common/build_telescope_detector.hpp"
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/common/build_wire_chamber.hpp"
#include "detray/test/cpu/toy_detector_test.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <filesystem>
#include <ios>

using namespace detray;

namespace {

/// Compare two files with names @param file_name1 and @param file_name2 for
/// equality, while skipping the first @param skip lines (header part)
bool compare_files(const std::string& file_name1, const std::string& file_name2,
                   std::size_t skip = 15u) {
  auto file1 =
      io::file_handle(file_name1, std::ios_base::in | std::ios_base::binary);
  auto file2 =
      io::file_handle(file_name2, std::ios_base::in | std::ios_base::binary);

  std::string line1;
  std::string line2;

  // Check files line by line
  std::size_t i{1u};
  while (std::getline(*file1, line1)) {
    if (std::getline(*file2, line2)) {
      if (skip < i && line1 != line2) {
        DETRAY_ERROR_HOST("In line " << i << ":" << std::endl
                                     << line1 << std::endl
                                     << line2);
        return false;
      }
    } else {
      DETRAY_ERROR_HOST("Could not read next line from file 2:"
                        << std::endl
                        << "In line " << i << ":" << std::endl
                        << line1);
      return false;
    }
    ++i;
  }

  // Are there more lines in file2 than file1?
  if (std::getline(*file2, line2)) {
    DETRAY_ERROR_HOST("Could not read next line from file 1:"
                      << std::endl
                      << "In line " << i << ":" << std::endl
                      << line2);
    return false;
  }

  // Passed
  return true;
}

/// Full IO round trip for a given detector
/// @returns a detector read back in from the writer files
template <std::size_t CAP = 0u, typename detector_t>
auto test_detector_json_io(
    const detector_t& det, const typename detector_t::name_map& names,
    std::map<std::string, std::string, std::less<>>& file_names,
    vecmem::host_memory_resource& host_mr) {
  auto writer_cfg = io::detector_writer_config{}
                        .format(io::format::json)
                        .replace_files(true)
                        .write_grids(true)
                        .write_material(true);
  io::write_detector(det, names, writer_cfg);

  // Read the detector back in
  io::detector_reader_config reader_cfg{};
  reader_cfg.verbose_check(true);
  for (auto& [_, name] : file_names) {
    reader_cfg.add_file(name);
  }

  auto [det2, names2] = io::read_detector<detector_t, CAP>(host_mr, reader_cfg);

  // Write the result to a different set of files
  writer_cfg.replace_files(false);
  io::write_detector(det2, names2, writer_cfg);

  // Compare writing round-trip
  std::string geometry_file{names.get_detector_name() + "_geometry_2.json"};
  EXPECT_TRUE(compare_files(file_names["geometry"], geometry_file));
  std::filesystem::remove(geometry_file);
  std::filesystem::remove(file_names["geometry"]);

  // Check a homogeneous material description, if present
  if (auto search = file_names.find("homogeneous_material");
      search != file_names.end()) {
    std::string hom_mat_file{names.get_detector_name() +
                             "_homogeneous_material_2.json"};
    EXPECT_TRUE(
        compare_files(file_names["homogeneous_material"], hom_mat_file));
    std::filesystem::remove(hom_mat_file);
    std::filesystem::remove(file_names["homogeneous_material"]);
  }

  // Check a material map description, if present
  if (auto search = file_names.find("material_maps");
      search != file_names.end()) {
    std::string mat_map_file{names.get_detector_name() +
                             "_material_maps_2.json"};
    EXPECT_TRUE(compare_files(file_names["material_maps"], mat_map_file));
    std::filesystem::remove(mat_map_file);
    std::filesystem::remove(file_names["material_maps"]);
  }

  // Check a homogeneous material description, if present
  if (auto search = file_names.find("surface_grids");
      search != file_names.end()) {
    std::string grids_file{names.get_detector_name() + "_surface_grids_2.json"};
    EXPECT_TRUE(compare_files(file_names["surface_grids"], grids_file));
    std::filesystem::remove(grids_file);
    std::filesystem::remove(file_names["surface_grids"]);
  }

  return std::make_pair(std::move(det2), std::move(names2));
}

}  // anonymous namespace

/// Test the reading and writing of a telescope detector
GTEST_TEST(io, json_telescope_detector_reader) {
  using test_algebra = test::algebra;
  using scalar = test::scalar;

  mask<rectangle2D, test_algebra> rec2{0u, 100.f, 100.f};

  // Surface positions
  std::vector<scalar> positions = {1.f,   50.f,  100.f, 150.f, 200.f, 250.f,
                                   300.f, 350.f, 400.f, 450.f, 500.f};

  tel_det_config tel_cfg{rec2};
  tel_cfg.positions(positions);

  // Telescope detector
  vecmem::host_memory_resource host_mr;
  auto [tel_det, tel_names] =
      build_telescope_detector<test_algebra>(host_mr, tel_cfg);

  std::map<std::string, std::string, std::less<>> file_names;
  file_names["geometry"] = "telescope_detector_geometry.json";
  file_names["homogeneous_material"] =
      "telescope_detector_homogeneous_material.json";

  auto [det_io, names_io] =
      test_detector_json_io(tel_det, tel_names, file_names, host_mr);

  const auto& mat_store = det_io.material_store();
  const auto& slabs =
      mat_store.get<decltype(tel_det)::material::id::e_material_slab>();

  EXPECT_EQ(det_io.volumes().size(), 1u);
  EXPECT_EQ(slabs.size(), positions.size());
}

/// Test the reading and writing of a toy detector geometry
GTEST_TEST(io, json_toy_geometry) {
  using metadata_t = test::toy_metadata;
  using test_algebra = metadata_t::algebra_type;
  using detector_t = detector<metadata_t>;
  using scalar = test::scalar;

  // Toy detector
  vecmem::host_memory_resource host_mr;
  toy_det_config<scalar> toy_cfg{};
  toy_cfg.use_material_maps(false);
  auto [toy_det, names] = build_toy_detector<test_algebra>(host_mr, toy_cfg);

  // Write the detector
  io::json_converter<detector_t, io::geometry_writer> geo_writer;
  auto file_name = geo_writer.write(
      toy_det, names, std::ios::out | std::ios::binary | std::ios::trunc);

  // Empty volume name map to be filled
  typename detector_t::name_map volume_name_map{};
  volume_name_map.set_detector_name("toy_detector");

  // Read the detector back in
  detector_builder<metadata_t> toy_builder;
  io::json_converter<detector_t, io::geometry_reader> geo_reader;
  geo_reader.read(toy_builder, file_name);
  auto det = toy_builder.build(host_mr, volume_name_map);

  // @TODO: Will only work again after IO can perform data deduplication
  // EXPECT_TRUE(toy_detector_test(det, volume_name_map));

  // Read the toy detector into the default detector type
  using default_metadata_t = test::default_metadata;
  detector_builder<default_metadata_t> comp_builder;
  io::json_converter<detector<default_metadata_t>, io::geometry_reader>
      comp_geo_reader;
  comp_geo_reader.read(comp_builder, file_name);
  volume_name_map.clear_names();
  auto comp_det = comp_builder.build(host_mr, volume_name_map);

  using mask_id = detector<default_metadata_t>::masks::id;
  const auto& masks = comp_det.mask_store();

  EXPECT_EQ(comp_det.volumes().size(), 22u);
  EXPECT_EQ(comp_det.surfaces().size(), 3230);
  EXPECT_EQ(comp_det.transform_store().size(), 3252);
  EXPECT_EQ(masks.template size<mask_id::e_rectangle2D>(), 2492u);
  EXPECT_EQ(masks.template size<mask_id::e_trapezoid2D>(), 648u);
  EXPECT_EQ(masks.template size<mask_id::e_annulus2D>(), 0u);
  EXPECT_EQ(masks.template size<mask_id::e_cylinder2D>(), 0u);
  EXPECT_EQ(masks.template size<mask_id::e_concentric_cylinder2D>(), 56u);
  EXPECT_EQ(masks.template size<mask_id::e_ring2D>(), 60u);
  EXPECT_EQ(masks.template size<mask_id::e_ring2D>(), 60u);
  EXPECT_EQ(masks.template size<mask_id::e_straw_tube>(), 0u);
  EXPECT_EQ(masks.template size<mask_id::e_drift_cell>(), 0u);

  detail::check_consistency(comp_det);
}

/// Test the reading and writing of a toy detector geometry "light"
GTEST_TEST(io, json_toy_detector_roundtrip_homogeneous_material) {
  using test_algebra = test::algebra;
  using scalar = test::scalar;

  // Toy detector
  vecmem::host_memory_resource host_mr;
  toy_det_config<scalar> toy_cfg{};
  toy_cfg.use_material_maps(false);
  const auto [toy_det, toy_names] =
      build_toy_detector<test_algebra>(host_mr, toy_cfg);

  std::map<std::string, std::string, std::less<>> file_names;
  file_names["geometry"] = "toy_detector_geometry.json";
  file_names["homogeneous_material"] = "toy_detector_homogeneous_material.json";
  file_names["surface_grids"] = "toy_detector_surface_grids.json";

  auto [det_io, names_io] =
      test_detector_json_io<1u>(toy_det, toy_names, file_names, host_mr);

  // Remove empty files as there are not material maps
  std::filesystem::remove("toy_detector_material_maps.json");
  std::filesystem::remove("toy_detector_material_maps_2.json");

  // @TODO: Will only work again after IO can perform data deduplication
  // EXPECT_TRUE(toy_detector_test(det_io, names_io));
}

/// Test the reading and writing of a toy detector geometry
GTEST_TEST(io, json_toy_detector_roundtrip_material_maps) {
  using test_algebra = test::algebra;
  using scalar = test::scalar;

  // Toy detector
  vecmem::host_memory_resource host_mr;
  toy_det_config<scalar> toy_cfg{};
  toy_cfg.use_material_maps(true);
  const auto [toy_det, toy_names] =
      build_toy_detector<test_algebra>(host_mr, toy_cfg);

  std::map<std::string, std::string, std::less<>> file_names;
  file_names["geometry"] = "toy_detector_geometry.json";
  file_names["homogeneous_material"] = "toy_detector_homogeneous_material.json";
  file_names["material_maps"] = "toy_detector_material_maps.json";
  file_names["surface_grids"] = "toy_detector_surface_grids.json";

  auto [det_io, names_io] =
      test_detector_json_io<1u>(toy_det, toy_names, file_names, host_mr);

  // @TODO: Will only work again after IO can perform data deduplication
  // EXPECT_TRUE(toy_detector_test(det_io, names_io));
}

/// Test the reading and writing of a wire chamber
GTEST_TEST(io, json_wire_chamber_reader) {
  using test_algebra = test::algebra;
  using scalar = test::scalar;

  // Wire chamber
  vecmem::host_memory_resource host_mr;
  wire_chamber_config<scalar> wire_cfg{};
  auto [wire_det, wire_names] =
      build_wire_chamber<test_algebra>(host_mr, wire_cfg);

  std::map<std::string, std::string, std::less<>> file_names;
  file_names["geometry"] = "wire_chamber_geometry.json";
  file_names["homogeneous_material"] = "wire_chamber_homogeneous_material.json";
  file_names["surface_grids"] = "wire_chamber_surface_grids.json";

  auto [det_io, names_io] =
      test_detector_json_io(wire_det, wire_names, file_names, host_mr);

  // Remove empty files as there are material
  std::filesystem::remove("wire_chamber_material_maps.json");
  std::filesystem::remove("wire_chamber_material_maps_2.json");

  EXPECT_EQ(det_io.volumes().size(), 11u);
}
