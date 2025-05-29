// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding2/CylindricalSpacePointGrid2.hpp"
#include "Acts/Seeding2/SeedFilter2.hpp"
#include "Acts/Seeding2/SeedFinder2.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace ActsExamples {

/// Construct track seeds from space points.
class SeedingAlgorithm2 final : public IAlgorithm {
 public:
  struct Config {
    /// Input space point collections.
    std::string inputSpacePoints;
    /// Output track seed collection.
    std::string outputSeeds;

    Acts::CylindricalSpacePointGrid2::Config gridConfig;

    Acts::SeedFinder2::Config finderConfig;
    Acts::SeedFinder2::Options finderOptions;
    Acts::SeedFilter2::Config filterConfig;

    // allow for different values of rMax in gridConfig and seedFinderConfig
    bool allowSeparateRMax = false;

    // vector containing the map of z bins in the top and bottom layers
    std::vector<std::pair<int, int>> zBinNeighborsTop;
    std::vector<std::pair<int, int>> zBinNeighborsBottom;
    // number of phiBin neighbors at each side of the current bin that will be
    // used to search for SPs
    int numPhiNeighbors = 1;

    // Connect custom selections on the space points or to the doublet
    // compatibility
    bool useExtraCuts = false;
  };

  /// Construct the seeding algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  SeedingAlgorithm2(const Config& cfg, Acts::Logging::Level lvl);

  /// Run the seeding algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  std::optional<Acts::SeedFinder2> m_seedFinder;
  std::unique_ptr<const Acts::GridBinFinder<3ul>> m_bottomBinFinder{nullptr};
  std::unique_ptr<const Acts::GridBinFinder<3ul>> m_topBinFinder{nullptr};

  Acts::Delegate<bool(const SimSpacePoint&)> m_spacePointSelector{
      Acts::DelegateFuncTag<voidSpacePointSelector>{}};

  static bool voidSpacePointSelector(const SimSpacePoint& /*sp*/) {
    return true;
  }

  ReadDataHandle<SimSpacePointContainer> m_inputSpacePoints{this,
                                                            "InputSpacePoints"};

  WriteDataHandle<SimSeedContainer> m_outputSeeds{this, "OutputSeeds"};

  static inline bool itkFastTrackingCuts(float bottomRadius, float cotTheta) {
    static float rMin = 45.;
    static float cotThetaMax = 1.5;

    if (bottomRadius < rMin &&
        (cotTheta > cotThetaMax || cotTheta < -cotThetaMax)) {
      return false;
    }
    return true;
  }

  static inline bool itkFastTrackingSPselect(const SimSpacePoint& sp) {
    // At small r we remove points beyond |z| > 200.
    float r = sp.r();
    float zabs = std::abs(sp.z());
    if (zabs > 200. && r < 45.) {
      return false;
    }

    /// Remove space points beyond eta=4 if their z is
    /// larger than the max seed z0 (150.)
    float cotTheta = 27.2899;  // corresponds to eta=4
    if ((zabs - 150.) > cotTheta * r) {
      return false;
    }
    return true;
  }
};

}  // namespace ActsExamples
