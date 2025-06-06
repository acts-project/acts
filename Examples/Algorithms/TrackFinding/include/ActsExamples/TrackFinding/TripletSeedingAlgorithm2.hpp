// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding2/CylindricalSpacePointGrid2.hpp"
#include "Acts/Seeding2/TripletSeedFilter2.hpp"
#include "Acts/Seeding2/TripletSeedFinder2.hpp"
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
class TripletSeedingAlgorithm2 final : public IAlgorithm {
 public:
  struct Config {
    /// Input space point collections.
    std::string inputSpacePoints;
    /// Output track seed collection.
    std::string outputSeeds;

    Acts::CylindricalSpacePointGrid2::Config gridConfig;

    Acts::TripletSeedFinder2::Config finderConfig;
    Acts::TripletSeedFinder2::Options finderOptions;
    Acts::TripletSeedFilter2::Config filterConfig;

    // Seeding parameters used in the space-point grid creation and bin finding

    /// Vector containing the z-bin edges for non equidistant binning in z
    std::vector<float> zBinEdges;

    // Seeding parameters used to define the region of interest for middle
    // space-point

    /// Radial range for middle space-point
    /// The range can be defined manually with (rMinMiddle, rMaxMiddle). If
    /// useVariableMiddleSPRange is set to false and the vector rRangeMiddleSP
    /// is empty, we use (rMinMiddle, rMaxMiddle) to cut the middle space-points
    float rMinMiddle = 60.f * Acts::UnitConstants::mm;
    float rMaxMiddle = 120.f * Acts::UnitConstants::mm;
    /// If useVariableMiddleSPRange is set to false, the vector rRangeMiddleSP
    /// can be used to define a fixed r range for each z bin: {{rMin, rMax},
    /// ...}
    bool useVariableMiddleSPRange = false;
    /// Range defined in vector for each z bin
    std::vector<std::vector<float>> rRangeMiddleSP;

    float deltaRMiddleMinSPRange = 10. * Acts::UnitConstants::mm;
    float deltaRMiddleMaxSPRange = 10. * Acts::UnitConstants::mm;

    /// Order of z bins to loop over when searching for SPs
    std::vector<std::size_t> zBinsCustomLooping;

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
  TripletSeedingAlgorithm2(const Config& cfg, Acts::Logging::Level lvl);

  /// Run the seeding algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  std::optional<Acts::TripletSeedFinder2> m_seedFinder;
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

  /// Get the proper radius validity range given a middle space point candidate.
  /// In case the radius range changes according to the z-bin we need to
  /// retrieve the proper range. We can do this computation only once, since all
  /// the middle space point candidates belong to the same z-bin
  /// @param spM space point candidate to be used as middle SP in a seed
  /// @param rMiddleSPRange range object containing the minimum and maximum r for middle SP for a certain z bin
  std::pair<float, float> retrieveRadiusRangeForMiddle(
      const Acts::ConstSpacePointProxy2& spM,
      const Acts::Range1D<float>& rMiddleSPRange) const;

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
