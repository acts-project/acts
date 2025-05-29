// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedingAlgorithm2.hpp"

#include "Acts/EventData2/SeedContainer2.hpp"
#include "Acts/EventData2/SpacePointContainer2.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"

#include <cmath>
#include <csignal>
#include <cstddef>
#include <limits>
#include <stdexcept>

using namespace Acts::HashedStringLiteral;

namespace ActsExamples {

namespace {

class SpacePointContainerProxy {
 public:
  explicit SpacePointContainerProxy(const SimSpacePointContainer& container)
      : m_container(container) {}

  std::size_t size() const { return m_container.size(); }

 private:
  const SimSpacePointContainer m_container;
};

}  // namespace

SeedingAlgorithm2::SeedingAlgorithm2(const Config& cfg,
                                     Acts::Logging::Level lvl)
    : IAlgorithm("SeedingAlgorithm2", lvl), m_cfg(cfg) {
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_outputSeeds.initialize(m_cfg.outputSeeds);

  if (m_cfg.gridConfig.rMax != m_cfg.finderConfig.rMax &&
      m_cfg.allowSeparateRMax == false) {
    throw std::invalid_argument(
        "Inconsistent config rMax: using different values in gridConfig and "
        "finderConfig. If values are intentional set allowSeparateRMax to "
        "true");
  }

  if (m_cfg.filterConfig.deltaRMin != m_cfg.finderConfig.deltaRMin) {
    throw std::invalid_argument("Inconsistent config deltaRMin");
  }

  if (m_cfg.gridConfig.deltaRMax != m_cfg.finderConfig.deltaRMax) {
    throw std::invalid_argument("Inconsistent config deltaRMax");
  }

  static_assert(std::numeric_limits<
                    decltype(m_cfg.finderConfig.deltaRMaxTopSP)>::has_quiet_NaN,
                "Value of deltaRMaxTopSP must support NaN values");

  static_assert(std::numeric_limits<
                    decltype(m_cfg.finderConfig.deltaRMinTopSP)>::has_quiet_NaN,
                "Value of deltaRMinTopSP must support NaN values");

  static_assert(
      std::numeric_limits<
          decltype(m_cfg.finderConfig.deltaRMaxBottomSP)>::has_quiet_NaN,
      "Value of deltaRMaxBottomSP must support NaN values");

  static_assert(
      std::numeric_limits<
          decltype(m_cfg.finderConfig.deltaRMinBottomSP)>::has_quiet_NaN,
      "Value of deltaRMinBottomSP must support NaN values");

  if (std::isnan(m_cfg.finderConfig.deltaRMaxTopSP)) {
    m_cfg.finderConfig.deltaRMaxTopSP = m_cfg.finderConfig.deltaRMax;
  }

  if (std::isnan(m_cfg.finderConfig.deltaRMinTopSP)) {
    m_cfg.finderConfig.deltaRMinTopSP = m_cfg.finderConfig.deltaRMin;
  }

  if (std::isnan(m_cfg.finderConfig.deltaRMaxBottomSP)) {
    m_cfg.finderConfig.deltaRMaxBottomSP = m_cfg.finderConfig.deltaRMax;
  }

  if (std::isnan(m_cfg.finderConfig.deltaRMinBottomSP)) {
    m_cfg.finderConfig.deltaRMinBottomSP = m_cfg.finderConfig.deltaRMin;
  }

  if (m_cfg.gridConfig.zMin != m_cfg.finderConfig.zMin) {
    throw std::invalid_argument("Inconsistent config zMin");
  }

  if (m_cfg.gridConfig.zMax != m_cfg.finderConfig.zMax) {
    throw std::invalid_argument("Inconsistent config zMax");
  }

  if (m_cfg.filterConfig.maxSeedsPerSpM != m_cfg.finderConfig.maxSeedsPerSpM) {
    throw std::invalid_argument("Inconsistent config maxSeedsPerSpM");
  }

  if (m_cfg.gridConfig.cotThetaMax != m_cfg.finderConfig.cotThetaMax) {
    throw std::invalid_argument("Inconsistent config cotThetaMax");
  }

  if (m_cfg.gridConfig.minPt != m_cfg.finderConfig.minPt) {
    throw std::invalid_argument("Inconsistent config minPt");
  }

  if (m_cfg.gridConfig.bFieldInZ != m_cfg.finderOptions.bFieldInZ) {
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

  if (!m_cfg.finderConfig.zBinsCustomLooping.empty()) {
    // check that the bins required in the custom bin looping
    // are contained in the bins defined by the total number of edges

    for (std::size_t i : m_cfg.finderConfig.zBinsCustomLooping) {
      if (i >= m_cfg.gridConfig.zBinEdges.size()) {
        throw std::invalid_argument(
            "Inconsistent config zBinsCustomLooping does not contain a subset "
            "of bins defined by zBinEdges");
      }
    }
  }

  if (m_cfg.useExtraCuts) {
    // This function will be applied to select space points during grid filling
    m_spacePointSelector.connect<itkFastTrackingSPselect>();

    // This function will be applied to the doublet compatibility selection
    m_cfg.finderConfig.experimentCuts.connect<itkFastTrackingCuts>();
  }

  m_bottomBinFinder = std::make_unique<const Acts::GridBinFinder<3ul>>(
      m_cfg.numPhiNeighbors, cfg.zBinNeighborsBottom, 0);
  m_topBinFinder = std::make_unique<const Acts::GridBinFinder<3ul>>(
      m_cfg.numPhiNeighbors, m_cfg.zBinNeighborsTop, 0);

  m_cfg.finderConfig.seedFilter = std::make_unique<Acts::SeedFilter2>(
      m_cfg.filterConfig.derive(), logger().cloneWithSuffix("SeedFilter2"));

  m_seedFinder = Acts::SeedFinder2(m_cfg.finderConfig.derive(),
                                   logger().cloneWithSuffix("SeedFinder2"));
}

ProcessCode SeedingAlgorithm2::execute(const AlgorithmContext& ctx) const {
  SimSpacePointContainer spacePoints = m_inputSpacePoints(ctx);

  Acts::SpacePointContainer2 coreSpacePoints;
  coreSpacePoints.reserve(spacePoints.size());
  for (const auto& sp : spacePoints) {
    // check if the space point passes the selection
    if (m_spacePointSelector(sp)) {
      auto newSp = coreSpacePoints.makeSpacePoint(sp.sourceLinks()[0], sp.x(),
                                                  sp.y(), sp.z());
      newSp.varianceR() = sp.varianceR();
      newSp.varianceZ() = sp.varianceZ();
    }
  }

  Acts::CylindricalSpacePointGrid2 grid(m_cfg.gridConfig.derive(),
                                        logger().cloneWithSuffix("Grid"));

  grid.fill(coreSpacePoints);

  // Compute radius Range
  // we rely on the fact the grid is storing the proxies
  // with a sorting in the radius
  const Acts::Range1D<float> rRange = grid.computeRadiusRange(coreSpacePoints);

  std::array<std::vector<std::size_t>, 3ul> navigation;
  navigation[1ul] = m_cfg.finderConfig.zBinsCustomLooping;

  auto spacePointsGrouping = std::move(grid).binnedGround(
      *m_bottomBinFinder, *m_topBinFinder, navigation);

  Acts::SeedFinder2::Options finderOptions = m_cfg.finderOptions;

  /// variable middle SP radial region of interest
  finderOptions.rMiddleSpRange = {
      std::floor(rRange.min() / 2) * 2 +
          m_cfg.finderConfig.deltaRMiddleMinSPRange,
      std::floor(rRange.max() / 2) * 2 -
          m_cfg.finderConfig.deltaRMiddleMaxSPRange};

  // run the seeding
  static thread_local Acts::SeedContainer2 seeds;
  static thread_local Acts::SeedFinder2::State state;

  auto derivedOptions = finderOptions.derive(m_seedFinder->config());

  std::vector<std::vector<Acts::SpacePointIndex2>> bottomSpGroups;
  std::vector<Acts::SpacePointIndex2> middleSpGroup;
  std::vector<std::vector<Acts::SpacePointIndex2>> topSpGroups;

  for (const auto [bottom, middle, top] : spacePointsGrouping) {
    ACTS_VERBOSE("Process middle " << middle);

    bottomSpGroups.clear();
    for (const auto b : bottom) {
      bottomSpGroups.push_back(spacePointsGrouping.grid().at(b));
    }
    middleSpGroup.clear();
    middleSpGroup = spacePointsGrouping.grid().at(middle);
    topSpGroups.clear();
    for (const auto t : top) {
      topSpGroups.push_back(spacePointsGrouping.grid().at(t));
    }

    m_seedFinder->createSeeds(derivedOptions, state, coreSpacePoints,
                              bottomSpGroups, middleSpGroup, topSpGroups,
                              seeds);
  }

  ACTS_DEBUG("Created " << seeds.size() << " track seeds from "
                        << spacePoints.size() << " space points");

  SimSeedContainer seedContainerForStorage;
  m_outputSeeds(ctx, std::move(seedContainerForStorage));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
