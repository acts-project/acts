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
template <typename external_spacepoint_t>
class SeedFilter final {
 public:
  explicit SeedFilter(
      const SeedFilterConfig& config,
      IExperimentCuts<external_spacepoint_t>* expCuts = nullptr);
  explicit SeedFilter(
      const SeedFilterConfig& config, std::unique_ptr<const Logger> logger,
      IExperimentCuts<external_spacepoint_t>* expCuts = nullptr);

  /// Create Seeds for the all seeds with the same bottom and middle
  /// space point and discard all others.
  /// @param mutableData Container for mutable variables used in the seeding
  /// @param bottomSp fixed bottom space point
  /// @param middleSp fixed middle space point
  /// @param topSpVec vector containing all space points that may be compatible
  ///                 with both bottom and middle space point
  /// @param invHelixDiameterVec vector containing 1/(2*r) values where r is the helix radius
  /// @param impactParametersVec vector containing the impact parameters
  /// @param seedFilterState holds quantities used in seed filter
  /// @param candidatesCollector container for the seed candidates
  void filterSeeds_2SpFixed(
      const SpacePointMutableData& mutableData,
      const external_spacepoint_t& bottomSp,
      const external_spacepoint_t& middleSp,
      const std::vector<const external_spacepoint_t*>& topSpVec,
      const std::vector<float>& invHelixDiameterVec,
      const std::vector<float>& impactParametersVec,
      SeedFilterState& seedFilterState,
      CandidatesForMiddleSp<const external_spacepoint_t>& candidatesCollector)
      const;

  /// Filter seeds once all seeds for one middle space point have been created
  /// @param mutableData Container for mutable variables used in the seeding
  /// @param candidatesCollector collection of seed candidates
  /// @param outputCollection Output container for the seeds
  /// for all seeds with the same middle space point
  template <typename collection_t>
  void filterSeeds_1SpFixed(
      SpacePointMutableData& mutableData,
      CandidatesForMiddleSp<const external_spacepoint_t>& candidatesCollector,
      collection_t& outputCollection) const;

  /// Filter seeds once all seeds for one middle space point have been created
  /// @param mutableData Container for mutable variables used in the seeding
  /// @param candidates collection of seed candidates
  /// @param numQualitySeeds number of high quality seeds in seed confirmation
  /// @param outputCollection Output container for the seeds
  /// for all seeds with the same middle space point
  template <typename collection_t>
  void filterSeeds_1SpFixed(
      SpacePointMutableData& mutableData,
      std::vector<typename CandidatesForMiddleSp<
          const external_spacepoint_t>::value_type>& candidates,
      const std::size_t numQualitySeeds, collection_t& outputCollection) const;

  const SeedFilterConfig getSeedFilterConfig() const { return m_cfg; }
  const IExperimentCuts<external_spacepoint_t>* getExperimentCuts() const {
    return m_experimentCuts;
  }

 private:
  const Logger& logger() const { return *m_logger; }

  const SeedFilterConfig m_cfg;
  std::unique_ptr<const Logger> m_logger =
      getDefaultLogger("Filter", Logging::Level::INFO);
  const IExperimentCuts<external_spacepoint_t>* m_experimentCuts;
};

}  // namespace Acts

#include "Acts/Seeding/SeedFilter.ipp"
