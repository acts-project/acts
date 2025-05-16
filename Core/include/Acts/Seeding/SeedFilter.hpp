// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointMutableData.hpp"
#include "Acts/Seeding/CandidatesForMiddleSp.hpp"
#include "Acts/Seeding/IExperimentCuts.hpp"
#include "Acts/Seeding/InternalSpacePointContainer.hpp"
#include "Acts/Seeding/SeedFilterConfig.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <vector>

namespace Acts {

struct SeedFilterState {
  // longitudinal impact parameter as defined by bottom and middle space point
  float zOrigin = 0;
  // number of minimum top SPs in seed confirmation
  std::size_t nTopSeedConf = 0;
  // radius of bottom component of seed that is used to define the number of
  // compatible top required
  float rMaxSeedConf = std::numeric_limits<float>::max();  // UnitConstants::mm
};

/// Filter seeds at various stages with the currently
/// available information.
class SeedFilter final {
 public:
  explicit SeedFilter(const SeedFilterConfig& config,
                      IExperimentCuts* expCuts = nullptr);
  explicit SeedFilter(const SeedFilterConfig& config,
                      std::unique_ptr<const Logger> logger,
                      IExperimentCuts* expCuts = nullptr);

  /// Create Seeds for the all seeds with the same bottom and middle
  /// space point and discard all others.
  /// @param mutableData Container for mutable variables used in the seeding
  /// @param bottomSP fixed bottom space point
  /// @param middleSP fixed middle space point
  /// @param topSpVec vector containing all space points that may be compatible
  ///                 with both bottom and middle space point
  /// @param invHelixDiameterVec vector containing 1/(2*r) values where r is the helix radius
  /// @param impactParametersVec vector containing the impact parameters
  /// @param seedFilterState holds quantities used in seed filter
  /// @param candidates_collector container for the seed candidates
  void filterSeeds_2SpFixed(const InternalSpacePointContainer& spacePoints,
                            const SpacePointMutableData& mutableData,
                            ConstInternalSpacePointProxy bottomSP,
                            ConstInternalSpacePointProxy middleSP,
                            const std::vector<std::size_t>& topSpVec,
                            const std::vector<float>& invHelixDiameterVec,
                            const std::vector<float>& impactParametersVec,
                            SeedFilterState& seedFilterState,
                            CandidatesForMiddleSp& candidates_collector) const;

  /// Filter seeds once all seeds for one middle space point have been created
  /// @param mutableData Container for mutable variables used in the seeding
  /// @param candidates_collector collection of seed candidates
  /// @param outputCollection Output container for the seeds
  /// for all seeds with the same middle space point
  template <typename external_spacepoint_t, typename collection_t>
  void filterSeeds_1SpFixed(const InternalSpacePointContainer& spacePoints,
                            SpacePointMutableData& mutableData,
                            CandidatesForMiddleSp& candidates_collector,
                            collection_t& outputCollection) const;

  /// Filter seeds once all seeds for one middle space point have been created
  /// @param mutableData Container for mutable variables used in the seeding
  /// @param candidates collection of seed candidates
  /// @param numQualitySeeds number of high quality seeds in seed confirmation
  /// @param outputCollection Output container for the seeds
  /// for all seeds with the same middle space point
  template <typename external_spacepoint_t, typename collection_t>
  void filterSeeds_1SpFixed(const InternalSpacePointContainer& spacePoints,
                            SpacePointMutableData& mutableData,
                            std::vector<TripletCandidate>& candidates,
                            const std::size_t numQualitySeeds,
                            collection_t& outputCollection) const;

  const SeedFilterConfig& getConfig() const { return m_cfg; }
  const IExperimentCuts* getExperimentCuts() const { return m_experimentCuts; }

 private:
  const Logger& logger() const { return *m_logger; }

  SeedFilterConfig m_cfg;
  std::unique_ptr<const Logger> m_logger =
      getDefaultLogger("Filter", Logging::Level::INFO);
  const IExperimentCuts* m_experimentCuts;
};

}  // namespace Acts

#include "Acts/Seeding/SeedFilter.ipp"
