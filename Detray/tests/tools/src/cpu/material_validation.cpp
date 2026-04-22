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
#include "detray/test/cpu/material_validation.hpp"

#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"

// Detray IO include(s)
#include "detray/io/frontend/detector_reader.hpp"

// Detray test include(s)
#include "detray/options/detector_io_options.hpp"
#include "detray/options/parse_options.hpp"
#include "detray/options/propagation_options.hpp"
#include "detray/options/track_generator_options.hpp"
#include "detray/test/cpu/material_scan.hpp"
#include "detray/test/framework/register_checks.hpp"
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// Boost
#include "detray/options/boost_program_options.hpp"

// System include(s)
#include <sstream>
#include <stdexcept>
#include <string>

namespace po = boost::program_options;

using namespace detray;

int main(int argc, char **argv) {
  // Use the most general type to be able to read in all detector files
  using metadata_t = test::default_metadata;
  using detector_t = detector<metadata_t>;

  // Filter out the google test flags
  ::testing::InitGoogleTest(&argc, argv);

  // Specific options for this test
  po::options_description desc("\ndetray material validation options");

  desc.add_options()(
      "data_dir",
      po::value<std::string>()->default_value("./validation_data/material"),
      "Directory that contains the data files")(
      "material_tol", po::value<float>()->default_value(1.f),
      "Tolerance for comparing the material traces [%]")(
      "overlaps_tol",
      po::value<float>()->default_value(stepping::config{}.min_stepsize),
      "Tolerance for considering surfaces to be overlapping [mm]");

  // Configs to be filled
  detray::io::detector_reader_config reader_cfg{};
  detray::test::material_validation<detector_t>::config mat_val_cfg{};
  test::material_scan<detector_t>::config mat_scan_cfg{};

  po::variables_map vm = detray::options::parse_options(
      desc, argc, argv, reader_cfg, mat_scan_cfg.track_generator(),
      mat_val_cfg.propagation());

  // General options
  if (vm.count("material_tol")) {
    mat_val_cfg.relative_error(vm["material_tol"].as<float>() / 100.f);
  }
  if (vm.count("overlaps_tol")) {
    mat_scan_cfg.overlaps_tol(vm["overlaps_tol"].as<float>());
  }
  const auto data_dir{vm["data_dir"].as<std::string>()};

  vecmem::host_memory_resource host_mr;

  const auto [det, names] =
      detray::io::read_detector<detector_t>(host_mr, reader_cfg);

  auto white_board = std::make_shared<test::whiteboard>();
  detector_t::geometry_context ctx{};

  // Print the detector's material as recorded by a ray scan
  mat_scan_cfg.track_generator().uniform_eta(true);
  mat_scan_cfg.material_file(data_dir + "/material_scan");

  detray::test::register_checks<test::material_scan>(det, names, mat_scan_cfg,
                                                     ctx, white_board);

  // Now trace the material during navigation and compare
  mat_val_cfg.material_file(data_dir + "/navigation_material_trace");

  test::register_checks<detray::test::material_validation>(
      det, names, mat_val_cfg, ctx, white_board);

  // Run the checks
  return RUN_ALL_TESTS();
}
