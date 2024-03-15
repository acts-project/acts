// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SpacePointData.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Seeding/BinnedGroup.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"

#include <cmath>
#include <csignal>
#include <cstddef>
#include <iterator>
#include <limits>
#include <ostream>
#include <stdexcept>

namespace ActsExamples {
struct AlgorithmContext;
}  // namespace ActsExamples

ActsExamples::SeedingAlgorithm::SeedingAlgorithm(
    ActsExamples::SeedingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("SeedingAlgorithm", lvl), m_cfg(std::move(cfg)) {
  // Seed Finder config requires Seed Filter object before conversion to
  // internal units
  m_cfg.seedFilterConfig = m_cfg.seedFilterConfig.toInternalUnits();
  m_cfg.seedFinderConfig.seedFilter =
      std::make_unique<Acts::SeedFilter<SimSpacePoint>>(m_cfg.seedFilterConfig);

  m_cfg.seedFinderConfig =
      m_cfg.seedFinderConfig.toInternalUnits().calculateDerivedQuantities();
  m_cfg.seedFinderOptions =
      m_cfg.seedFinderOptions.toInternalUnits().calculateDerivedQuantities(
          m_cfg.seedFinderConfig);
  m_cfg.gridConfig = m_cfg.gridConfig.toInternalUnits();
  m_cfg.gridOptions = m_cfg.gridOptions.toInternalUnits();
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point input collections");
  }

  for (const auto& spName : m_cfg.inputSpacePoints) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto& handle = m_inputSpacePoints.emplace_back(
        std::make_unique<ReadDataHandle<SimSpacePointContainer>>(
            this,
            "InputSpacePoints#" + std::to_string(m_inputSpacePoints.size())));
    handle->initialize(spName);
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds output collection");
  }

  m_outputSeeds.initialize(m_cfg.outputSeeds);

  if (m_cfg.gridConfig.rMax != m_cfg.seedFinderConfig.rMax &&
      m_cfg.allowSeparateRMax == false) {
    throw std::invalid_argument(
        "Inconsistent config rMax: using different values in gridConfig and "
        "seedFinderConfig. If values are intentional set allowSeparateRMax to "
        "true");
  }

  if (m_cfg.seedFilterConfig.deltaRMin != m_cfg.seedFinderConfig.deltaRMin) {
    throw std::invalid_argument("Inconsistent config deltaRMin");
  }

  if (m_cfg.gridConfig.deltaRMax != m_cfg.seedFinderConfig.deltaRMax) {
    throw std::invalid_argument("Inconsistent config deltaRMax");
  }

  static_assert(
      std::numeric_limits<
          decltype(m_cfg.seedFinderConfig.deltaRMaxTopSP)>::has_quiet_NaN,
      "Value of deltaRMaxTopSP must support NaN values");

  static_assert(
      std::numeric_limits<
          decltype(m_cfg.seedFinderConfig.deltaRMinTopSP)>::has_quiet_NaN,
      "Value of deltaRMinTopSP must support NaN values");

  static_assert(
      std::numeric_limits<
          decltype(m_cfg.seedFinderConfig.deltaRMaxBottomSP)>::has_quiet_NaN,
      "Value of deltaRMaxBottomSP must support NaN values");

  static_assert(
      std::numeric_limits<
          decltype(m_cfg.seedFinderConfig.deltaRMinBottomSP)>::has_quiet_NaN,
      "Value of deltaRMinBottomSP must support NaN values");

  if (std::isnan(m_cfg.seedFinderConfig.deltaRMaxTopSP)) {
    m_cfg.seedFinderConfig.deltaRMaxTopSP = m_cfg.seedFinderConfig.deltaRMax;
  }

  if (std::isnan(m_cfg.seedFinderConfig.deltaRMinTopSP)) {
    m_cfg.seedFinderConfig.deltaRMinTopSP = m_cfg.seedFinderConfig.deltaRMin;
  }

  if (std::isnan(m_cfg.seedFinderConfig.deltaRMaxBottomSP)) {
    m_cfg.seedFinderConfig.deltaRMaxBottomSP = m_cfg.seedFinderConfig.deltaRMax;
  }

  if (std::isnan(m_cfg.seedFinderConfig.deltaRMinBottomSP)) {
    m_cfg.seedFinderConfig.deltaRMinBottomSP = m_cfg.seedFinderConfig.deltaRMin;
  }

  if (m_cfg.gridConfig.zMin != m_cfg.seedFinderConfig.zMin) {
    throw std::invalid_argument("Inconsistent config zMin");
  }

  if (m_cfg.gridConfig.zMax != m_cfg.seedFinderConfig.zMax) {
    throw std::invalid_argument("Inconsistent config zMax");
  }

  if (m_cfg.seedFilterConfig.maxSeedsPerSpM !=
      m_cfg.seedFinderConfig.maxSeedsPerSpM) {
    throw std::invalid_argument("Inconsistent config maxSeedsPerSpM");
  }

  if (m_cfg.gridConfig.cotThetaMax != m_cfg.seedFinderConfig.cotThetaMax) {
    throw std::invalid_argument("Inconsistent config cotThetaMax");
  }

  if (m_cfg.gridConfig.minPt != m_cfg.seedFinderConfig.minPt) {
    throw std::invalid_argument("Inconsistent config minPt");
  }

  if (m_cfg.gridOptions.bFieldInZ != m_cfg.seedFinderOptions.bFieldInZ) {
    throw std::invalid_argument("Inconsistent config bFieldInZ");
  }

  if (m_cfg.gridConfig.zBinEdges.size() - 1 != m_cfg.zBinNeighborsTop.size() &&
      m_cfg.zBinNeighborsTop.empty() == false) {
    throw std::invalid_argument("Inconsistent config zBinNeighborsTop");
  }

  if (m_cfg.gridConfig.zBinEdges.size() - 1 !=
          m_cfg.zBinNeighborsBottom.size() &&
      m_cfg.zBinNeighborsBottom.empty() == false) {
    throw std::invalid_argument("Inconsistent config zBinNeighborsBottom");
  }

  if (!m_cfg.seedFinderConfig.zBinsCustomLooping.empty()) {
    // check if zBinsCustomLooping contains numbers from 1 to the total number
    // of bin in zBinEdges
    for (std::size_t i = 1; i != m_cfg.gridConfig.zBinEdges.size(); i++) {
      if (std::find(m_cfg.seedFinderConfig.zBinsCustomLooping.begin(),
                    m_cfg.seedFinderConfig.zBinsCustomLooping.end(),
                    i) == m_cfg.seedFinderConfig.zBinsCustomLooping.end()) {
        throw std::invalid_argument(
            "Inconsistent config zBinsCustomLooping does not contain the same "
            "bins as zBinEdges");
      }
    }
  }

  if (m_cfg.seedFinderConfig.useDetailedDoubleMeasurementInfo) {
    m_cfg.seedFinderConfig.getTopHalfStripLength.connect(
        [](const void*, const SimSpacePoint& sp) -> float {
          return sp.topHalfStripLength();
        });

    m_cfg.seedFinderConfig.getBottomHalfStripLength.connect(
        [](const void*, const SimSpacePoint& sp) -> float {
          return sp.bottomHalfStripLength();
        });

    m_cfg.seedFinderConfig.getTopStripDirection.connect(
        [](const void*, const SimSpacePoint& sp) -> Acts::Vector3 {
          return sp.topStripDirection();
        });

    m_cfg.seedFinderConfig.getBottomStripDirection.connect(
        [](const void*, const SimSpacePoint& sp) -> Acts::Vector3 {
          return sp.bottomStripDirection();
        });

    m_cfg.seedFinderConfig.getStripCenterDistance.connect(
        [](const void*, const SimSpacePoint& sp) -> Acts::Vector3 {
          return sp.stripCenterDistance();
        });

    m_cfg.seedFinderConfig.getTopStripCenterPosition.connect(
        [](const void*, const SimSpacePoint& sp) -> Acts::Vector3 {
          return sp.topStripCenterPosition();
        });
  }

  m_bottomBinFinder = std::make_unique<const Acts::GridBinFinder<2ul>>(
      m_cfg.numPhiNeighbors, m_cfg.zBinNeighborsBottom);
  m_topBinFinder = std::make_unique<const Acts::GridBinFinder<2ul>>(
      m_cfg.numPhiNeighbors, m_cfg.zBinNeighborsTop);

  m_cfg.seedFinderConfig.seedFilter =
      std::make_unique<Acts::SeedFilter<SimSpacePoint>>(m_cfg.seedFilterConfig);
  m_seedFinder =
      Acts::SeedFinder<SimSpacePoint,
                       Acts::CylindricalSpacePointGrid<SimSpacePoint>>(
          m_cfg.seedFinderConfig);
}

ActsExamples::ProcessCode ActsExamples::SeedingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // construct the combined input container of space point pointers from all
  // configured input sources.
  // pre-compute the total size required so we only need to allocate once
  std::size_t nSpacePoints = 0;
  for (const auto& isp : m_inputSpacePoints) {
    nSpacePoints += (*isp)(ctx).size();
  }

  std::vector<const SimSpacePoint*> spacePointPtrs;
  spacePointPtrs.reserve(nSpacePoints);
  for (const auto& isp : m_inputSpacePoints) {
    for (const auto& spacePoint : (*isp)(ctx)) {
      // since the event store owns the space
      // points, their pointers should be stable and
      // we do not need to create local copies.
      spacePointPtrs.push_back(&spacePoint);
    }
  }

  // construct the seeding tools
  // covariance tool, extracts covariances per spacepoint as required
  auto extractGlobalQuantities = [=](const SimSpacePoint& sp, float, float,
                                     float) {
    Acts::Vector3 position{sp.x(), sp.y(), sp.z()};
    Acts::Vector2 covariance{sp.varianceR(), sp.varianceZ()};
    return std::make_tuple(position, covariance, sp.t());
  };

  // extent used to store r range for middle spacepoint
  Acts::Extent rRangeSPExtent;

  Acts::CylindricalSpacePointGrid<SimSpacePoint> grid =
      Acts::CylindricalSpacePointGridCreator::createGrid<SimSpacePoint>(
          m_cfg.gridConfig, m_cfg.gridOptions);
  Acts::CylindricalSpacePointGridCreator::fillGrid(
      m_cfg.seedFinderConfig, m_cfg.seedFinderOptions, grid,
      spacePointPtrs.begin(), spacePointPtrs.end(), extractGlobalQuantities,
      rRangeSPExtent);

  std::array<std::vector<std::size_t>, 2ul> navigation;
  navigation[1ul] = m_cfg.seedFinderConfig.zBinsCustomLooping;

  auto spacePointsGrouping = Acts::CylindricalBinnedGroup<SimSpacePoint>(
      std::move(grid), *m_bottomBinFinder, *m_topBinFinder,
      std::move(navigation));

  // safely clamp double to float
  float up = Acts::clampValue<float>(
      std::floor(rRangeSPExtent.max(Acts::binR) / 2) * 2);

  /// variable middle SP radial region of interest
  const Acts::Range1D<float> rMiddleSPRange(
      std::floor(rRangeSPExtent.min(Acts::binR) / 2) * 2 +
          m_cfg.seedFinderConfig.deltaRMiddleMinSPRange,
      up - m_cfg.seedFinderConfig.deltaRMiddleMaxSPRange);

  // run the seeding
  static thread_local SimSeedContainer seeds;
  seeds.clear();
  static thread_local decltype(m_seedFinder)::SeedingState state;
  state.spacePointData.resize(
      spacePointPtrs.size(),
      m_cfg.seedFinderConfig.useDetailedDoubleMeasurementInfo);

  if (m_cfg.seedFinderConfig.useDetailedDoubleMeasurementInfo) {
    for (std::size_t grid_glob_bin(0);
         grid_glob_bin < spacePointsGrouping.grid().size(); ++grid_glob_bin) {
      const auto& collection = spacePointsGrouping.grid().at(grid_glob_bin);
      for (const auto& sp : collection) {
        std::size_t index = sp->index();

        const float topHalfStripLength =
            m_cfg.seedFinderConfig.getTopHalfStripLength(sp->sp());
        const float bottomHalfStripLength =
            m_cfg.seedFinderConfig.getBottomHalfStripLength(sp->sp());
        const Acts::Vector3 topStripDirection =
            m_cfg.seedFinderConfig.getTopStripDirection(sp->sp());
        const Acts::Vector3 bottomStripDirection =
            m_cfg.seedFinderConfig.getBottomStripDirection(sp->sp());

        state.spacePointData.setTopStripVector(
            index, topHalfStripLength * topStripDirection);
        state.spacePointData.setBottomStripVector(
            index, bottomHalfStripLength * bottomStripDirection);
        state.spacePointData.setStripCenterDistance(
            index, m_cfg.seedFinderConfig.getStripCenterDistance(sp->sp()));
        state.spacePointData.setTopStripCenterPosition(
            index, m_cfg.seedFinderConfig.getTopStripCenterPosition(sp->sp()));
      }
    }
  }

  for (const auto [bottom, middle, top] : spacePointsGrouping) {
    m_seedFinder.createSeedsForGroup(
        m_cfg.seedFinderOptions, state, spacePointsGrouping.grid(),
        std::back_inserter(seeds), bottom, middle, top, rMiddleSPRange);
  }

  ACTS_DEBUG("Created " << seeds.size() << " track seeds from "
                        << spacePointPtrs.size() << " space points");

  m_outputSeeds(ctx, SimSeedContainer{seeds});
  return ActsExamples::ProcessCode::SUCCESS;
}
