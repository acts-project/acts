// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <string>
#include <unordered_map>

#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "Acts/Geometry/GeometryID.hpp"

namespace FW {

/// Create planar clusters from simulation hits.
template <typename seedfinder_t>
class SeedingAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Create the configuration struct with the explicit seed finder
    Config(const seedfinder_t& sfinder) : seedFinder(sfinder) {}
    /// Disallow the default construct
    Config() = delete;
    /// Input collection of space points hits.
    std::string inputSpacePoints = "";
    /// Output collection of clusters.
    std::string outputSeeds = "";
    /// The actual seed finder
    const seedfinder_t& seedFinder;
  };

  /// Construct the digitization algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  SeedingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Build clusters from input simulation hits.
  ///
  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final override;

 private:
  Config m_cfg;
};

}  // namespace FW
