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
#include "detray/navigation/volume_graph.hpp"
#include "detray/utils/logging.hpp"

// Detray IO include(s)
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/io/utils/file_handle.hpp"

// Example linear algebra plugin: std::array
#include "detray/tutorial/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <ios>
#include <iostream>
#include <stdexcept>
#include <string>

/// Read a detector from file and transform its volumes and navigation links to
/// a graph in dot format
int main(int argc, char** argv) {
  std::clog << "Detector Graph Tutorial\n=======================\n";

  // Input data file
  auto reader_cfg = detray::io::detector_reader_config{};
  if (argc == 2) {
    reader_cfg.add_file(argv[1]);
  } else {
    throw std::runtime_error("Please specify an input file name!");
  }

  std::clog << reader_cfg << std::endl;

  // Read a toy detector
  using metadata_t = detray::tutorial::toy_metadata;
  using detector_t = detray::detector<metadata_t>;

  // Create an empty detector to be filled
  vecmem::host_memory_resource host_mr;

  // Read the detector in
  const auto [det, names] =
      detray::io::read_detector<detector_t>(host_mr, reader_cfg);

  // Display the detector volume graph
  detray::volume_graph graph(det);

  const std::string file_stem{det.name(names) + "_dot"};
  const std::ios_base::openmode io_mode{std::ios::out | std::ios::trunc};
  detray::io::file_handle out_file{file_stem, ".txt", io_mode};

  *out_file << graph.to_dot_string() << std::endl;

  DETRAY_INFO_HOST("Wrote file: " << file_stem + ".txt");
}
