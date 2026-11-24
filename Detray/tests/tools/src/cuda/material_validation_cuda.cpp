// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"

// Detray IO include(s)
#include "detray/io/json/detector_reader.hpp"

// Detray test include(s)
#include "detray/options/detector_io_options.hpp"
#include "detray/options/parse_options.hpp"
#include "detray/options/propagation_options.hpp"
#include "detray/options/track_generator_options.hpp"
#include "detray/test/device/cuda/material_validation.hpp"
#include "detray/test/validation/material_validation.hpp"
#include "detray/test/validation/material_validation_config.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
// Boost
#include "detray/options/boost_program_options.hpp"

// System include(s)
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace po = boost::program_options;

using namespace detray;

int main(int argc, char **argv) {
  // Use the most general type to be able to read in all detector files
  using metadata_t = test::default_metadata;
  using detector_t = host::detector<metadata_t>;

  // Specific options for this test
  po::options_description desc("\ndetray CUDA material validation options");

  desc.add_options()(
      "data_dir",
      po::value<std::string>()->default_value("./validation_data/material"),
      "Directory that contains the data files")(
      "material_tol", po::value<float>()->default_value(1.f),
      "Tolerance for comparing the material traces [%]")(
      "overlaps_tol",
      po::value<float>()->default_value(stepping::config{}.min_stepsize),
      "Tolerance for considering surfaces to be overlapping [mm]")(
      "step_contingency",
      po::value<std::size_t>()->default_value(
          test::material_validation_config<test::algebra>{}.step_contingency()),
      "Additional material steps per track on top of the material scan")(
      "truth_from_file",
      "Read the truth material from file instead of running the ray scan");

  // Configs to be filled
  detray::io::detector_reader_config reader_cfg{};
  test::material_validation_config<test::algebra> mat_val_cfg{};
  test::material_scan_config<test::algebra> mat_scan_cfg{};

  po::variables_map vm = detray::options::parse_options(
      desc, argc, argv, reader_cfg, mat_scan_cfg.track_generator(),
      mat_val_cfg.propagation());

  // General options
  if (vm.count("material_tol") != 0u) {
    mat_val_cfg.relative_error(vm["material_tol"].as<float>() / 100.f);
  }
  if (vm.count("overlaps_tol") != 0u) {
    mat_scan_cfg.overlaps_tol(vm["overlaps_tol"].as<float>());
  }
  mat_val_cfg.step_contingency(vm["step_contingency"].as<std::size_t>());
  mat_scan_cfg.truth_from_file(vm.count("truth_from_file") != 0u);
  const auto data_dir{vm["data_dir"].as<std::string>()};

  /// Vecmem memory resource for the device allocations
  vecmem::cuda::device_memory_resource dev_mr{};
  vecmem::host_memory_resource host_mr;

  const auto [det, names] =
      detray::io::read_detector_json<detector_t>(host_mr, reader_cfg);

  detector_t::geometry_context ctx{};

  // Print the detector's material as recorded by a ray scan
  mat_scan_cfg.name("material_scan_for_cuda");
  mat_scan_cfg.track_generator().uniform_eta(true);
  mat_scan_cfg.material_file(data_dir + "/material_scan_for_cuda");
  const auto scan_data =
      material_validator::material_scan(ctx, det, names, mat_scan_cfg);

  // Now trace the material during navigation and compare
  mat_val_cfg.name("material_validation_cuda");
  mat_val_cfg.material_file(data_dir + "/navigation_material_trace");
  mat_val_cfg.device_mr(&dev_mr);

  const auto result =
      material_validator::material_validation<cuda::run_material_validation>(
          &host_mr, det, names, scan_data, mat_val_cfg);

  if (!result.passed()) {
    std::cerr << result.n_failed << " of " << result.n_tracks
              << " tracks have a material discrepancy" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
