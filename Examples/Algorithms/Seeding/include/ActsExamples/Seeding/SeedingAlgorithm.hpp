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
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"

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
    // input Clusters from the event#-hits.csv file.
    std::string inputClusters;
    // used to get truth information into seeds about what particles are in what
    // space point.
    std::string inputHitParticlesMap;

    float rMax = 200.;
    float deltaRMin = 1.;
    float deltaRMax = 60.;
    float collisionRegionMin = -250;
    float collisionRegionMax = 250.;
    float zMin = -2000.;
    float zMax = 2000.;
    float maxSeedsPerSpM = 1;
    float cotThetaMax = 7.40627;  // 2.7 eta
    float sigmaScattering = 2.25;
    float radLengthPerSeed = 0.1;
    float minPt = 500.;
    float bFieldInZ = 0.00199724;
    Acts::Vector2D beamPos = {0., 0.};
    float impactMax = 3.;
    std::vector<int> seedVolumes = {7, 8, 9};

    Acts::SeedfinderConfig<SimSpacePoint> finderConf;

    Acts::SeedFilterConfig sfconf;

    Acts::SpacePointGridConfig gridConf;
  };

  /// Construct the digitization algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  SeedingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// @param cluster The hit with local hit information
  /// Returns a space point with a particle barcode stored in .particles for
  /// each particle that made this space point.
  std::unique_ptr<ActsExamples::SimSpacePoint> transformSP(
      std::size_t hit_id, const Acts::PlanarModuleCluster& cluster,
      const AlgorithmContext& ctx) const;

  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final override;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
