// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/navigation/volume_graph.hpp"
#include "detray/utils/logging.hpp"

// Detray test include(s)
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/utils/hash_tree.hpp"

// Example linear algebra plugin: std::array
#include "detray/tutorial/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <iostream>

// Hash of the "correct" geometry
constexpr std::size_t root_hash = 9563506807235398581ul;

///
/// Work in progress (!)
///
/// Get a graph that represents the detector volumes (nodes) and their links via
/// their adjacent boundary surfaces (edges) and hash it to detect linking
/// changes that might affect the navigation.
/// This could be useful in a CI job, but is poorly tested at this point (!).
int main() {
  std::clog << "Volume Graph Tutorial\n=====================\n\n";

  // Get an example detector
  vecmem::host_memory_resource host_mr;
  const auto [det, names] =
      detray::build_toy_detector<detray::tutorial::algebra_t>(host_mr);

  // Build the graph and get its adjacency matrix
  detray::volume_graph graph(det);
  const auto &adj_mat = graph.adjacency_matrix();

  // Construct a hash tree on the graph and compare against existing hash to
  // detect changes
  auto geo_checker = detray::hash_tree(adj_mat);

  if (geo_checker.root() == root_hash) {
    DETRAY_INFO_HOST("Geometry links are consistent");
  } else {
    DETRAY_ERROR_HOST("Geometry linking has changed! Root hash is "
                      << geo_checker.root() << ". Please check geometry:\n"
                      << graph.to_string());
  }
}
