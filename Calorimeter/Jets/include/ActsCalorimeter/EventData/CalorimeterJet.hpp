// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <vector>

namespace ActsCalorimeter {

/// A reconstructed calorimeter jet, independent of FastJet object lifetimes.
struct CalorimeterJet {
  /// E-scheme sum in ACTS native units, ordered (px, py, pz, E).
  /// Individual cluster inputs are massless; a recombined jet can have mass.
  Acts::Vector4 fourMomentum = Acts::Vector4::Zero();
  /// Sorted indices into this event's input CalorimeterCluster collection.
  /// Follow cluster.hitIndices and hit.sourceIndices to reach deposits.
  /// These transient links require the named input collections to be retained.
  std::vector<std::size_t> clusterIndices;
};

using CalorimeterJetContainer = std::vector<CalorimeterJet>;

}  // namespace ActsCalorimeter
