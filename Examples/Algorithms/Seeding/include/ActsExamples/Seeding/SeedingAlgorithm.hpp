// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "ActsExamples//Framework/BareAlgorithm.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/Seeding/SimSpacePoint.hpp"

#include <memory>
#include <set>
#include <string>
#include <unordered_map>

namespace Acts {
class PlanarModuleStepper;
}  // namespace Acts

namespace ActsExamples {

/// Create planar clusters from simulation hits.
class SeedingAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Output collection of clusters.
    std::string outputSeeds;
    /// Output prototracks
    std::string outputProtoTracks;
    // input Clusters from the event#-hits.csv file.
    std::string inputClusters;
    // input particles for creating proto seeds
    std::string inputParticles;
    // used to get truth information into seeds about what particles are in what
    // space point.
    std::string inputHitParticlesMap;
  };

  /// Construct the digitization algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  SeedingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// @param cluster The hit with local hit information
  /// Returns a space point with a particle barcode stored in .particles for
  /// each particle that made this space point.
  ActsExamples::SimSpacePoint* transformSP(
      std::size_t hit_id, const Acts::GeometryIdentifier geoId,
      const Acts::PlanarModuleCluster& cluster,
      const IndexMultimap<ActsFatras::Barcode>& hitParticlesMap,
      const AlgorithmContext& ctx) const;

  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final override;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
