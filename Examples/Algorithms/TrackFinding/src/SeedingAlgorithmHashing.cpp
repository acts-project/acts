// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedingAlgorithmHashing.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SpacePointData.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Seeding/BinnedGroup.hpp"
#include "Acts/Seeding/Hashing/HashingAlgorithm.hpp"
#include "Acts/Seeding/Hashing/HashingTraining.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedFinder.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

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

ActsExamples::SeedingAlgorithmHashing::SeedingAlgorithmHashing(
    ActsExamples::SeedingAlgorithmHashing::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("SeedingAlgorithmHashing", lvl),
      m_cfg(std::move(cfg)) {
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
  m_outputBuckets.initialize(m_cfg.outputBuckets);

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
  m_Hashing = Acts::HashingAlgorithm<const SimSpacePoint*,
                                     std::vector<const SimSpacePoint*>>(
      m_cfg.hashingConfig);
  m_HashingTraining =
      Acts::HashingTrainingAlgorithm<std::vector<const SimSpacePoint*>>(
          m_cfg.hashingTrainingConfig);
}

ActsExamples::ProcessCode ActsExamples::SeedingAlgorithmHashing::execute(
    const AlgorithmContext& ctx) const {
  ACTS_DEBUG("Start of SeedingAlgorithmHashing execute");
  using SpacePointPtrVector = std::vector<const SimSpacePoint*>;

  // construct the combined input container of space point pointers from all
  // configured input sources.
  // pre-compute the total size required so we only need to allocate once
  std::size_t nSpacePoints = 0;
  for (const auto& isp : m_inputSpacePoints) {
    nSpacePoints += (*isp)(ctx).size();
  }

  SpacePointPtrVector spacePointPtrs;
  spacePointPtrs.reserve(nSpacePoints);
  for (const auto& isp : m_inputSpacePoints) {
    for (const SimSpacePoint& spacePoint : (*isp)(ctx)) {
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

  // Hashing Training
  Acts::AnnoyModel annoyModel = m_HashingTraining.execute(spacePointPtrs);
  // Hashing
  static thread_local std::vector<SpacePointPtrVector> bucketsPtrs;
  bucketsPtrs.clear();
  VectorPolicy bucketsPolicy(bucketsPtrs);
  GenericBackInserter bucketsBackInserter(bucketsPolicy);
  m_Hashing.execute(spacePointPtrs, &annoyModel, bucketsBackInserter);

  // pre-compute the maximum size required so we only need to allocate once
  // doesn't combine the input containers of space point pointers
  std::size_t maxNSpacePoints = 0, inSpacePoints = 0;
  for (const SpacePointPtrVector& bucket : bucketsPtrs) {
    inSpacePoints = bucket.size();
    if (inSpacePoints > maxNSpacePoints) {
      maxNSpacePoints = inSpacePoints;
    }
  }

  static thread_local std::set<ActsExamples::SimSeed> seedsSet;
  seedsSet.clear();
  static thread_local decltype(m_seedFinder)::SeedingState state;
  state.spacePointData.resize(
      maxNSpacePoints, m_cfg.seedFinderConfig.useDetailedDoubleMeasurementInfo);

  for (const SpacePointPtrVector& bucket : bucketsPtrs) {
    // extent used to store r range for middle spacepoint
    Acts::Extent rRangeSPExtent;
    // construct the seeding tools
    Acts::CylindricalSpacePointGrid<SimSpacePoint> grid =
        Acts::CylindricalSpacePointGridCreator::createGrid<SimSpacePoint>(
            m_cfg.gridConfig, m_cfg.gridOptions);
    Acts::CylindricalSpacePointGridCreator::fillGrid(
        m_cfg.seedFinderConfig, m_cfg.seedFinderOptions, grid, bucket.begin(),
        bucket.end(), extractGlobalQuantities, rRangeSPExtent);

    std::array<std::vector<std::size_t>, 2ul> navigation;
    navigation[1ul] = m_cfg.seedFinderConfig.zBinsCustomLooping;

    state.spacePointData.clear();

    // groups spacepoints
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
              index,
              m_cfg.seedFinderConfig.getTopStripCenterPosition(sp->sp()));
        }
      }
    }

    for (const auto [bottom, middle, top] : spacePointsGrouping) {
      SetPolicy setPolicy(seedsSet);
      GenericBackInserter backInserter(setPolicy);
      m_seedFinder.createSeedsForGroup(m_cfg.seedFinderOptions, state,
                                       spacePointsGrouping.grid(), backInserter,
                                       bottom, middle, top, rMiddleSPRange);
    }
  }

  static thread_local SimSeedContainer seeds;
  seeds.clear();
  seeds.reserve(seedsSet.size());
  for (const ActsExamples::SimSeed& seed : seedsSet) {
    seeds.push_back(seed);
  }

  ACTS_INFO("Created " << seeds.size() << " track seeds from "
                       << spacePointPtrs.size() << " space points");

  m_outputSeeds(ctx, SimSeedContainer{seeds});
  std::vector<SimSpacePointContainer> buckets;
  for (const SpacePointPtrVector& bucket : bucketsPtrs) {
    SimSpacePointContainer bucketSP;
    for (const SimSpacePoint* spacePoint : bucket) {
      bucketSP.push_back(*spacePoint);
    }
    buckets.push_back(bucketSP);
  }
  m_outputBuckets(ctx, std::vector<SimSpacePointContainer>{buckets});

  ACTS_DEBUG("End of SeedingAlgorithmHashing execute");
  return ActsExamples::ProcessCode::SUCCESS;
}
