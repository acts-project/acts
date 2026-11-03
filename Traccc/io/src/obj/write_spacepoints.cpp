// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "write_spacepoints.hpp"

// System include(s).
#include <fstream>

namespace traccc::io::obj {

void write_spacepoints(
    std::string_view filename,
    edm::spacepoint_collection::const_view spacepoints_view) {
  // Open the output file.
  std::ofstream file(filename.data());
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + std::string(filename));
  }

  // Create a device collection around the spacepoint view.
  const edm::spacepoint_collection::const_device spacepoints(spacepoints_view);

  // Write the spacepoints.
  for (edm::spacepoint_collection::const_device::size_type i = 0u;
       i < spacepoints.size(); ++i) {
    const edm::spacepoint sp = spacepoints.at(i);
    file << "v " << sp.x() << " " << sp.y() << " " << sp.z() << "\n";
  }
}

}  // namespace traccc::io::obj
