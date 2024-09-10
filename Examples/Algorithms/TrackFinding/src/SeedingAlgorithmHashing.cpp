// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedingAlgorithmHashing.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/Seed.hpp"
#include "Acts/EventData/SpacePointData.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Plugins/Hashing/HashingAlgorithm.hpp"
#include "Acts/Plugins/Hashing/HashingTraining.hpp"
#include "Acts/Seeding/BinnedGroup.hpp"
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

namespace ActsExamples {
struct AlgorithmContext;
}  // namespace ActsExamples

ActsExamples::SeedingAlgorithmHashing::SeedingAlgorithmHashing(
    ActsExamples::SeedingAlgorithmHashing::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("SeedingAlgorithmHashing", lvl),
      m_cfg(std::move(cfg)) {
  using SpacePointProxy_type = typename Acts::SpacePointContainer<
      ActsExamples::SpacePointContainer<std::vector<const SimSpacePoint*>>,
      Acts::detail::RefHolder>::SpacePointProxyType;

  // Seed Finder config requires Seed Filter object before conversion to
  // internal units
  m_cfg.seedFilterConfig = m_cfg.seedFilterConfig.toInternalUnits();
  m_cfg.seedFinderConfig.seedFilter =
      std::make_unique<Acts::SeedFilter<SpacePointProxy_type>>(
          m_cfg.seedFilterConfig);

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

  if (m_cfg.useExtraCuts) {
    // This function will be applied to select space points during grid filling
    m_cfg.seedFinderConfig.spacePointSelector
        .connect<itkFastTrackingSPselect>();

    // This function will be applied to the doublet compatibility selection
    m_cfg.seedFinderConfig.experimentCuts.connect<itkFastTrackingCuts>();
  }

  m_bottomBinFinder = std::make_unique<const Acts::GridBinFinder<2ul>>(
      m_cfg.numPhiNeighbors, m_cfg.zBinNeighborsBottom);
  m_topBinFinder = std::make_unique<const Acts::GridBinFinder<2ul>>(
      m_cfg.numPhiNeighbors, m_cfg.zBinNeighborsTop);

  m_cfg.seedFinderConfig.seedFilter =
      std::make_unique<Acts::SeedFilter<SpacePointProxy_type>>(
          m_cfg.seedFilterConfig);
  m_seedFinder =
      Acts::SeedFinder<SpacePointProxy_type,
                       Acts::CylindricalSpacePointGrid<SpacePointProxy_type>>(
          m_cfg.seedFinderConfig);
  m_hashing = Acts::HashingAlgorithm<const SpacePointProxy_type*,
                                     std::vector<const SpacePointProxy_type*>>(
      m_cfg.hashingConfig);
  m_hashingTraining =
      Acts::HashingTrainingAlgorithm<std::vector<const SpacePointProxy_type*>>(
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

  // Config
  Acts::SpacePointContainerConfig spConfig;
  spConfig.useDetailedDoubleMeasurementInfo =
      m_cfg.seedFinderConfig.useDetailedDoubleMeasurementInfo;
  // Options
  Acts::SpacePointContainerOptions spOptions;
  spOptions.beamPos = {0., 0.};

  // Prepare interface SpacePoint backend-ACTS
  ActsExamples::SpacePointContainer container(spacePointPtrs);
  // Prepare Acts API
  Acts::SpacePointContainer<decltype(container), Acts::detail::RefHolder>
      spContainer(spConfig, spOptions, container);

  using value_type = typename Acts::SpacePointContainer<
      ActsExamples::SpacePointContainer<std::vector<const SimSpacePoint*>>,
      Acts::detail::RefHolder>::SpacePointProxyType;
  using seed_type = Acts::Seed<value_type>;

  std::vector<const value_type*> spacePointProxyPtr;
  for (const value_type& proxy : spContainer) {
    spacePointProxyPtr.push_back(&proxy);
  }

  // Hashing Training
  Acts::AnnoyModel annoyModel = m_hashingTraining.execute(spacePointProxyPtr);
  // Hashing
  static thread_local std::vector<std::vector<const value_type*>> bucketsPtrs;
  bucketsPtrs.clear();
  m_hashing.execute(spacePointProxyPtr, &annoyModel, bucketsPtrs);

  // pre-compute the maximum size required so we only need to allocate once
  // doesn't combine the input containers of space point pointers
  std::size_t maxNSpacePoints = 0;
  std::size_t inSpacePoints = 0;
  for (const std::vector<const value_type*>& bucket : bucketsPtrs) {
    inSpacePoints = bucket.size();
    if (inSpacePoints > maxNSpacePoints) {
      maxNSpacePoints = inSpacePoints;
    }
  }

  // Create the set with custom comparison function
  static thread_local std::set<seed_type, SeedComparison<value_type>> seedsSet;
  seedsSet.clear();
  static thread_local decltype(m_seedFinder)::SeedingState state;
  state.spacePointMutableData.resize(spContainer.size());

  for (const std::vector<const value_type*>& bucket : bucketsPtrs) {
    std::vector<value_type> buck;
    buck.reserve(bucket.size());
    for (const value_type* val : bucket) {
      buck.push_back(*val);
    }

    // extent used to store r range for middle spacepoint
    Acts::Extent rRangeSPExtent;
    // construct the seeding tools
    Acts::CylindricalSpacePointGrid<value_type> grid =
        Acts::CylindricalSpacePointGridCreator::createGrid<value_type>(
            m_cfg.gridConfig, m_cfg.gridOptions);
    Acts::CylindricalSpacePointGridCreator::fillGrid(
        m_cfg.seedFinderConfig, m_cfg.seedFinderOptions, grid, buck.begin(),
        buck.end(), rRangeSPExtent);

    std::array<std::vector<std::size_t>, 2ul> navigation;
    navigation[1ul] = m_cfg.seedFinderConfig.zBinsCustomLooping;

    // groups spacepoints
    auto spacePointsGrouping = Acts::CylindricalBinnedGroup<value_type>(
        std::move(grid), *m_bottomBinFinder, *m_topBinFinder,
        std::move(navigation));

    // safely clamp double to float
    float up = Acts::clampValue<float>(
        std::floor(rRangeSPExtent.max(Acts::BinningValue::binR) / 2) * 2);

    /// variable middle SP radial region of interest
    const Acts::Range1D<float> rMiddleSPRange(
        std::floor(rRangeSPExtent.min(Acts::BinningValue::binR) / 2) * 2 +
            m_cfg.seedFinderConfig.deltaRMiddleMinSPRange,
        up - m_cfg.seedFinderConfig.deltaRMiddleMaxSPRange);

    for (const auto [bottom, middle, top] : spacePointsGrouping) {
      m_seedFinder.createSeedsForGroup(m_cfg.seedFinderOptions, state,
                                       spacePointsGrouping.grid(), seedsSet,
                                       bottom, middle, top, rMiddleSPRange);
    }
  }
  /*
  // convert seed of proxies to seed of simseeds
  ActsExamples::SimSeedContainer seeds;
  seeds.clear();
  seeds.reserve(seedsSet.size());
  for (const seed_type& seed : seedsSet) {
    const std::array<const value_type*, 3>& sps = seed.sp();
    const SimSpacePoint* bottom = sps[0]->externalSpacePoint();
    const SimSpacePoint* middle = sps[1]->externalSpacePoint();
    const SimSpacePoint* top = sps[2]->externalSpacePoint();

    seeds.emplace_back(*bottom, *middle, *top);
    seeds.back().setZvertex(seed.z());
    seeds.back().setQuality(seed.seedQuality());
  }

  ACTS_INFO("Created " << seeds.size() << " track seeds from "
                       << spacePointPtrs.size() << " space points");

  m_outputSeeds(ctx, SimSeedContainer{seeds});
  std::vector<SimSpacePointContainer> buckets;
  for (const std::vector<value_type>& bucket : bucketsPtrs) {
    SimSpacePointContainer bucketSP;
    for (const value_type& spacePoint : bucket) {
      bucketSP.push_back(*spacePoint.externalSpacePoint());
    }
    buckets.push_back(bucketSP);
  }
  m_outputBuckets(ctx, std::vector<SimSpacePointContainer>{buckets});
  */
  ACTS_DEBUG("End of SeedingAlgorithmHashing execute");
  return ActsExamples::ProcessCode::SUCCESS;
}
