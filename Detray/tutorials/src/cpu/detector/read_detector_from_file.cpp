// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/navigation/volume_graph.hpp"
#include "detray/utils/logging.hpp"

// Example linear algebra plugin: std::array
#include "detray/tutorial/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <iostream>
#include <stdexcept>
#include <string>

/// Read a detector from file. For now: Read in the geometry by calling the
/// json geometry reader directly.
int main(int argc, char** argv) {
  std::clog << "Detector IO Tutorial\n====================\n";

  // Input data file
  auto reader_cfg = detray::io::detector_reader_config{};
  if (argc == 2) {
    reader_cfg.add_file(argv[1]);
  } else {
    throw std::runtime_error("Please specify an input file name!");
  }

  std::clog << reader_cfg << std::endl;

  // Read a toy detector
  using metadata_t = detray::tutorial::default_metadata;
  using detector_t = detray::detector<metadata_t>;

  // Create an empty detector to be filled
  vecmem::host_memory_resource host_mr;

  // Read the json files
  const auto [det, names] =
      detray::io::read_detector<detector_t>(host_mr, reader_cfg);

  // Print the detector volume graph
  detray::volume_graph graph(det);
  DETRAY_INFO_HOST("Read " << det.volumes().size() << " volumes from file "
                           << reader_cfg.files()[0u] << ":\n\n"
                           << graph.to_string());
}
