// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/builders/volume_builder_interface.hpp"
#include "detray/core/concepts.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/utils/consistency_checker.hpp"
#include "detray/utils/logging.hpp"

// Detray IO include(s)
#include "detray/io/backend/geometry_reader.hpp"
#include "detray/io/backend/geometry_writer.hpp"
#include "detray/io/json/detector_reader.hpp"
#include "detray/io/json/detector_writer.hpp"
#include "detray/io/json/json_converter.hpp"

// Detray test include(s)
#include "detray/test/common/build_nested_detector.hpp"
#include "detray/test/common/build_telescope_detector.hpp"
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/common/build_wire_chamber.hpp"
#include "detray/test/cpu/toy_detector_test.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <algorithm>
#include <filesystem>
#include <ios>
#include <vector>

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
template <std::size_t CAP = 0u, concepts::detector detector_t>
auto test_detector_json_io(
    const detector_t& det, const typename detector_t::name_map& names,
    std::map<std::string, std::string, std::less<>>& file_names,
    vecmem::host_memory_resource& host_mr, const bool deduplicate = true) {
  auto writer_cfg = io::detector_writer_config{}
                        .format(io::format::json)
                        .source("detray integration test")
                        .replace_files(true)
                        .write_grids(true)
                        .write_material(true);
  io::write_detector_json(det, names, writer_cfg);

  // Read the detector back in
  io::detector_reader_config reader_cfg{};
  reader_cfg.verbose_check(true).deduplicate(deduplicate);
  for (auto& [_, name] : file_names) {
    reader_cfg.add_file(name);
  }

  auto [det2, names2] =
      io::read_detector_json<detector_t, CAP>(host_mr, reader_cfg);

  // Write the result to a different set of files
  writer_cfg.replace_files(false);
  io::write_detector_json(det2, names2, writer_cfg);

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

  constexpr bool deduplicate{false};
  mask<rectangle2D, test_algebra> rec2{0u, -100.f, -100.f, 100.f, 100.f};

  // Surface positions
  std::vector<scalar> positions = {1.f,   50.f,  100.f, 150.f, 200.f, 250.f,
                                   300.f, 350.f, 400.f, 450.f, 500.f};

  tel_det_config tel_cfg{rec2};
  tel_cfg.positions(positions);
  tel_cfg.builder_options().deduplicate(deduplicate);

  // Telescope detector
  vecmem::host_memory_resource host_mr;
  auto [tel_det, tel_names] =
      build_telescope_detector<test_algebra>(host_mr, tel_cfg);

  std::map<std::string, std::string, std::less<>> file_names;
  file_names["geometry"] = "telescope_detector_geometry.json";
  file_names["homogeneous_material"] =
      "telescope_detector_homogeneous_material.json";

  auto [det_io, names_io] = test_detector_json_io(
      tel_det, tel_names, file_names, host_mr, deduplicate);

  const auto& mat_store = det_io.material_store();
  const auto& slabs =
      mat_store.get<decltype(tel_det)::material::id::e_material_slab>();

  EXPECT_EQ(det_io.volumes().size(), 1u);
  EXPECT_EQ(slabs.size(), positions.size());
}

/// Test the reading and writing of a telescope detector - deduplicated
GTEST_TEST(io, json_telescope_detector_reader_deduplicated) {
  using test_algebra = test::algebra;
  using scalar = test::scalar;

  constexpr bool deduplicate{true};
  mask<rectangle2D, test_algebra> rec2{0u, -100.f, -100.f, 100.f, 100.f};

  // Surface positions
  std::vector<scalar> positions = {1.f,   50.f,  100.f, 150.f, 200.f, 250.f,
                                   300.f, 350.f, 400.f, 450.f, 500.f};

  tel_det_config tel_cfg{rec2};
  tel_cfg.builder_options().deduplicate(deduplicate);

  // Telescope detector
  vecmem::host_memory_resource host_mr;
  auto [tel_det, tel_names] =
      build_telescope_detector<test_algebra>(host_mr, tel_cfg);

  std::map<std::string, std::string, std::less<>> file_names;
  file_names["geometry"] = "telescope_detector_geometry.json";
  file_names["homogeneous_material"] =
      "telescope_detector_homogeneous_material.json";

  auto [det_io, names_io] = test_detector_json_io(
      tel_det, tel_names, file_names, host_mr, deduplicate);

  const auto& mat_store = det_io.material_store();
  const auto& slabs =
      mat_store.get<decltype(tel_det)::material::id::e_material_slab>();

  EXPECT_EQ(det_io.volumes().size(), 1u);
  EXPECT_EQ(slabs.size(), 1u);
}

/// Collect the daughter volume indices of a volume
struct daughter_collector {
  void operator()(const dindex daughter_idx,
                  std::vector<dindex>& daughters) const {
    daughters.push_back(daughter_idx);
  }
};

/// Test the reading and writing of a detector with daughter volumes
GTEST_TEST(io, json_nested_detector_roundtrip) {
  using test_algebra = test::algebra;
  using point3 = test::point3;

  vecmem::host_memory_resource host_mr;
  const auto [nested_det, nested_names] =
      build_nested_detector<test_algebra>(host_mr);

  std::map<std::string, std::string, std::less<>> file_names;
  file_names["geometry"] = "nested_detector_geometry.json";

  auto [det_io, names_io] =
      test_detector_json_io(nested_det, nested_names, file_names, host_mr);

  ASSERT_EQ(det_io.volumes().size(), nested_det.volumes().size());

  // Same mother and daughter volumes
  auto get_daughters = [](const auto& det, const dindex vol_idx) {
    std::vector<dindex> daughters{};
    const auto vol = tracking_volume{det, vol_idx};
    if (vol.has_daughters()) {
      vol.template visit_daughter_volumes<daughter_collector>(daughters);
    }
    std::ranges::sort(daughters);
    return daughters;
  };

  for (const auto& vol_desc : nested_det.volumes()) {
    const dindex vol_idx{vol_desc.index()};
    const auto& vol_desc_io = det_io.volume(vol_idx);

    EXPECT_EQ(vol_desc_io.mother(), vol_desc.mother())
        << nested_names.at(vol_idx);
    EXPECT_EQ(get_daughters(det_io, vol_idx),
              get_daughters(nested_det, vol_idx))
        << nested_names.at(vol_idx);
  }

  // Same volume lookup for points in all hierarchy levels
  const std::vector<point3> points{{0.f, 0.f, 0.f},     {0.f, 0.f, -100.f},
                                   {0.f, 0.f, -130.f},  {0.f, 0.f, 160.f},
                                   {0.f, 130.f, 0.f},   {0.f, 130.f, 60.f},
                                   {150.f, 0.f, 150.f}, {250.f, 0.f, 0.f}};

  for (const auto& [i, p] : detray::views::enumerate(points)) {
    EXPECT_EQ(det_io.volume_index(p), nested_det.volume_index(p))
        << "point #" << i;
  }
}

/// Test that broken mother links in the geometry payload are rejected
GTEST_TEST(io, json_nested_detector_mother_links) {
  using test_algebra = test::algebra;
  using toy_detector_t = host::detector<toy_metadata<test_algebra>>;

  vecmem::host_memory_resource host_mr;
  const auto [nested_det, nested_names] =
      build_nested_detector<test_algebra>(host_mr);
  using nested_detector_t = std::remove_cvref_t<decltype(nested_det)>;

  // Fill the payload from the nested detector
  io::detector_payload payload{};
  auto writer_cfg = io::detector_writer_config{}.source("detray test");
  io::write_detector(nested_det, nested_names, writer_cfg, payload);

  // The mother links restore the daughter volumes (default brute force)
  EXPECT_NO_THROW(io::read_detector<nested_detector_t>(
      host_mr, io::detector_reader_config{}, payload));

  // Mother volume does not exist
  io::detector_payload broken_payload{payload};
  broken_payload.geometry.volumes.at(2u).mother->link = 100u;
  EXPECT_THROW(io::read_detector<nested_detector_t>(
                   host_mr, io::detector_reader_config{}, broken_payload),
               std::invalid_argument);

  // Wrong mother volume: Box G (in box A) is not contained in box B
  broken_payload = payload;
  broken_payload.geometry.volumes.at(2u).mother->link = 3u;
  EXPECT_THROW(io::read_detector<nested_detector_t>(
                   host_mr, io::detector_reader_config{}, broken_payload),
               std::invalid_argument);

  // The toy detector type cannot hold daughter volumes (only a warning,
  // until volume grids exist)
  EXPECT_NO_THROW(io::read_detector<toy_detector_t>(
      host_mr, io::detector_reader_config{}, payload));
}

/// Test the reading and writing of a toy detector geometry
GTEST_TEST(io, json_toy_geometry) {
  using metadata_t = test::toy_metadata;
  using test_algebra = metadata_t::algebra_type;
  using detector_t = host::detector<metadata_t>;
  using scalar = test::scalar;

  // Toy detector
  vecmem::host_memory_resource host_mr;
  toy_det_config<scalar> toy_cfg{};
  toy_cfg.use_material_maps(false);
  auto [toy_det, names] = build_toy_detector<test_algebra>(host_mr, toy_cfg);

  // Write the detector
  io::detector_payload out_payload{};
  out_payload.header.detector = toy_det.name(names);
  out_payload.header.source = "detray integration test";
  io::geometry_writer<detector_t>{}.to_payload(toy_det, names, out_payload);
  out_payload.header.geometry = io::to_header_payload(out_payload.geometry);

  const auto file_name =
      io::json_output_converter{}
          .write_to_file(out_payload, "./",
                         std::ios::out | std::ios::binary | std::ios::trunc)
          .front();

  // Read the detector back in
  io::detector_payload payload{};
  io::json_input_converter{}.add_to_payload(file_name, payload);

  detector_builder<metadata_t> toy_builder;
  io::geometry_reader<detector_t>{}.from_payload(toy_builder, payload);

  names.clear_names();
  volume_builder_options builder_opts{};
  auto det = toy_builder.build(host_mr, builder_opts, names);

  // @TODO: Will only work again after IO can perform data deduplication
  // EXPECT_TRUE(toy_detector_test(det, volume_name_map));

  // Read the toy detector into the default detector type
  using default_metadata_t = test::default_metadata;
  detector_builder<default_metadata_t> comp_builder;
  io::geometry_reader<host::detector<default_metadata_t>>{}.from_payload(
      comp_builder, payload);
  names.clear_names();
  auto comp_det = comp_builder.build(host_mr, builder_opts, names);

  using mask_id = host::detector<default_metadata_t>::masks::id;
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
  EXPECT_FALSE(io::file_exists("toy_detector_material_maps.json"));
  EXPECT_FALSE(io::file_exists("toy_detector_material_maps_2.json"));

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
  file_names["material_maps"] = "toy_detector_material_maps.json";
  file_names["surface_grids"] = "toy_detector_surface_grids.json";

  auto [det_io, names_io] =
      test_detector_json_io<1u>(toy_det, toy_names, file_names, host_mr);

  // Remove empty files as there are not material maps
  EXPECT_FALSE(io::file_exists("toy_detector_homogeneous_material.json"));
  EXPECT_FALSE(io::file_exists("toy_detector_homogeneous_material_2.json"));

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
  EXPECT_FALSE(io::file_exists("wire_chamber_material_maps.json"));
  EXPECT_FALSE(io::file_exists("wire_chamber_material_maps_2.json"));

  EXPECT_EQ(det_io.volumes().size(), 11u);
}
