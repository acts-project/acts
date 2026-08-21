// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/GridTripletSeedingAlgorithm.hpp"

#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding/DoubletSeedFinder.hpp"
#include "Acts/Seeding/TripletSeedFinder.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"

#include <cmath>
#include <csignal>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ActsExamples {

namespace {

static inline bool itkFastTrackingCuts(
    const Acts::ConstSpacePointProxy& /*middle*/,
    const Acts::ConstSpacePointProxy& other, float cotTheta,
    bool isBottomCandidate) {
  static float rMin = 45;
  static float cotThetaMax = 1.5;

  if (isBottomCandidate && other.zr()[1] < rMin &&
      (cotTheta > cotThetaMax || cotTheta < -cotThetaMax)) {
    return false;
  }
  return true;
}

static inline bool itkFastTrackingSPselect(const ConstSpacePointProxy& sp) {
  // At small r we remove points beyond |z| > 200.
  float r = sp.r();
  float zabs = std::abs(sp.z());
  if (zabs > 200. && r < 45.) {
    return false;
  }

  // Remove space points beyond eta=4 if their z is larger than the max seed
  // z0 (150.)
  float cotTheta = 27.2899;  // corresponds to eta=4
  if ((zabs - 150.) > cotTheta * r) {
    return false;
  }
  return true;
}

/// Stateful experiment-cut functor that rejects a doublet whose z-origin does
/// not fall inside any of the provided per-event z-windows (built from the
/// reconstructed vertices). Bound to the DoubletSeedFinder's `experimentCuts`
/// delegate for the duration of one event. An empty window list keeps every
/// doublet (no constraint).
struct VertexZCuts {
  std::vector<std::pair<float, float>> windows;  // [zmin, zmax] in mm
  // per-event diagnostics (mutable: operator() is const, called during finding)
  mutable std::size_t nTested{0};
  mutable std::size_t nRejected{0};

  bool operator()(const Acts::ConstSpacePointProxy& middle,
                  const Acts::ConstSpacePointProxy& /*other*/, float cotTheta,
                  bool /*isBottomCandidate*/) const {
    if (windows.empty()) {
      return true;
    }
    ++nTested;
    // z-origin of the doublet extrapolated to the beam line. NOTE: the core
    // seeding SpacePointContainer fills the packed zr column (z = zr[0],
    // r = zr[1]), not the separate z()/r() columns -- calling middle.z()/r()
    // here reads unallocated data and segfaults.
    const float zM = middle.zr()[0];
    const float rM = middle.zr()[1];
    const float zOrigin = zM - rM * cotTheta;
    for (const auto& [lo, hi] : windows) {
      if (zOrigin >= lo && zOrigin <= hi) {
        return true;
      }
    }
    ++nRejected;
    return false;  // outside every window -> reject
  }
};

}  // namespace

GridTripletSeedingAlgorithm::GridTripletSeedingAlgorithm(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("GridTripletSeedingAlgorithm", std::move(logger)), m_cfg(cfg) {
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_outputSeeds.initialize(m_cfg.outputSeeds);
  if (!m_cfg.inputVertices.empty()) {
    m_inputVertices.initialize(m_cfg.inputVertices);
    // Note: the ctor parameter is named `logger`, which shadows the logger()
    // method the ACTS_* macros use, so log via this->logger() explicitly.
    ACTS_LOG_WITH_LOGGER(this->logger(), Acts::Logging::INFO,
                         "Vertex-z seeding constraint ENABLED: vertices='"
                             << m_cfg.inputVertices
                             << "', nSigma=" << m_cfg.vertexZNSigma
                             << ", margin=" << m_cfg.vertexZMargin << " mm");
  }

  // check that the bins required in the custom bin looping
  // are contained in the bins defined by the total number of edges
  for (std::size_t i : m_cfg.zBinsCustomLooping) {
    if (i >= m_cfg.zBinEdges.size()) {
      throw std::invalid_argument(
          "Inconsistent config zBinsCustomLooping does not contain a subset "
          "of bins defined by zBinEdges");
    }
  }

  if (m_cfg.useExtraCuts) {
    // This function will be applied to select space points during grid filling
    m_spacePointSelector.connect<itkFastTrackingSPselect>();
  }

  m_gridConfig.minPt = m_cfg.minPt;
  // TODO switch to `m_cfg.rMin`
  m_gridConfig.rMin = 0;
  m_gridConfig.rMax = m_cfg.rMax;
  m_gridConfig.zMin = m_cfg.zMin;
  m_gridConfig.zMax = m_cfg.zMax;
  m_gridConfig.deltaRMax = m_cfg.deltaRMax;
  m_gridConfig.cotThetaMax = m_cfg.cotThetaMax;
  m_gridConfig.impactMax = m_cfg.impactMax;
  m_gridConfig.phiMin = m_cfg.phiMin;
  m_gridConfig.phiMax = m_cfg.phiMax;
  m_gridConfig.phiBinDeflectionCoverage = m_cfg.phiBinDeflectionCoverage;
  m_gridConfig.maxPhiBins = m_cfg.maxPhiBins;
  m_gridConfig.rBinEdges = {};
  m_gridConfig.zBinEdges = m_cfg.zBinEdges;
  m_gridConfig.bFieldInZ = m_cfg.bFieldInZ;
  m_gridConfig.bottomBinFinder.emplace(m_cfg.numPhiNeighbors,
                                       m_cfg.zBinNeighborsBottom, 0);
  m_gridConfig.topBinFinder.emplace(m_cfg.numPhiNeighbors,
                                    m_cfg.zBinNeighborsTop, 0);
  m_gridConfig.navigation[0ul] = {};
  m_gridConfig.navigation[1ul] = m_cfg.zBinsCustomLooping;
  m_gridConfig.navigation[2ul] = {};

  m_filterConfig.deltaInvHelixDiameter = m_cfg.deltaInvHelixDiameter;
  m_filterConfig.deltaRMin = m_cfg.deltaRMin;
  m_filterConfig.compatSeedWeight = m_cfg.compatSeedWeight;
  m_filterConfig.impactWeightFactor = m_cfg.impactWeightFactor;
  m_filterConfig.zOriginWeightFactor = m_cfg.zOriginWeightFactor;
  m_filterConfig.maxSeedsPerSpM = m_cfg.maxSeedsPerSpM;
  m_filterConfig.compatSeedLimit = m_cfg.compatSeedLimit;
  m_filterConfig.seedWeightIncrement = m_cfg.seedWeightIncrement;
  m_filterConfig.numSeedIncrement = m_cfg.numSeedIncrement;
  m_filterConfig.seedConfirmation = m_cfg.seedConfirmation;
  m_filterConfig.centralSeedConfirmationRange =
      m_cfg.centralSeedConfirmationRange;
  m_filterConfig.forwardSeedConfirmationRange =
      m_cfg.forwardSeedConfirmationRange;
  m_filterConfig.maxSeedsPerSpMConf = m_cfg.maxSeedsPerSpMConf;
  m_filterConfig.maxQualitySeedsPerSpMConf = m_cfg.maxQualitySeedsPerSpMConf;
  m_filterConfig.useDeltaRinsteadOfTopRadius =
      m_cfg.useDeltaRinsteadOfTopRadius;

  m_filterLogger = this->logger().cloneWithSuffix("Filter");

  m_seedFinder = Acts::TripletSeeder(this->logger().cloneWithSuffix("Finder"));
}

ProcessCode GridTripletSeedingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const SpacePointContainer& spacePoints = m_inputSpacePoints(ctx);

  // Build the per-event vertex z-windows used to constrain the doublet
  // z-origin. Lives on the stack for the whole event; the DoubletSeedFinders
  // (below) hold a delegate pointing at it and are only used within this call.
  VertexZCuts vertexZCuts;
  if (!m_cfg.inputVertices.empty()) {
    const VertexContainer& vertices = m_inputVertices(ctx);
    vertexZCuts.windows.reserve(vertices.size());
    for (const auto& v : vertices) {
      const double z = v.position().z();
      const double sigmaZ = std::sqrt(v.covariance()(2, 2));
      const double half = m_cfg.vertexZNSigma * sigmaZ + m_cfg.vertexZMargin;
      vertexZCuts.windows.emplace_back(static_cast<float>(z - half),
                                       static_cast<float>(z + half));
    }
    std::ostringstream oss;
    for (const auto& [lo, hi] : vertexZCuts.windows) {
      oss << " [" << lo << ", " << hi << "]";
    }
    ACTS_VERBOSE("Vertex-z seeding constraint: read "
                 << vertices.size() << " vertex(es) from '"
                 << m_cfg.inputVertices << "' -> " << vertexZCuts.windows.size()
                 << " z-window(s) [mm]:" << oss.str());
  }

  Acts::CylindricalSpacePointGrid grid(m_gridConfig,
                                       logger().cloneWithSuffix("Grid"));

  for (std::size_t i = 0; i < spacePoints.size(); ++i) {
    const auto& sp = spacePoints[i];

    // check if the space point passes the selection
    if (m_spacePointSelector.connected() && !m_spacePointSelector(sp)) {
      continue;
    }

    float phi = std::atan2(sp.y(), sp.x());
    grid.insert(i, phi, sp.z(), sp.r());
  }

  for (std::size_t i = 0; i < grid.numberOfBins(); ++i) {
    std::ranges::sort(grid.at(i), [&](const Acts::SpacePointIndex& a,
                                      const Acts::SpacePointIndex& b) {
      return spacePoints[a].r() < spacePoints[b].r();
    });
  }

  Acts::SpacePointContainer coreSpacePoints(
      Acts::SpacePointColumns::CopiedFromIndex |
      Acts::SpacePointColumns::PackedXY | Acts::SpacePointColumns::PackedZR |
      Acts::SpacePointColumns::VarianceZ | Acts::SpacePointColumns::VarianceR);
  coreSpacePoints.reserve(grid.numberOfSpacePoints());
  std::vector<Acts::SpacePointIndexRange> gridSpacePointRanges;
  gridSpacePointRanges.reserve(grid.numberOfBins());
  for (std::size_t i = 0; i < grid.numberOfBins(); ++i) {
    std::uint32_t begin = coreSpacePoints.size();
    for (Acts::SpacePointIndex spIndex : grid.at(i)) {
      const ConstSpacePointProxy& sp = spacePoints[spIndex];

      auto newSp = coreSpacePoints.createSpacePoint();
      newSp.copiedFromIndex() = sp.index();
      newSp.xy() = std::array<float, 2>{static_cast<float>(sp.x()),
                                        static_cast<float>(sp.y())};
      newSp.zr() = std::array<float, 2>{static_cast<float>(sp.z()),
                                        static_cast<float>(sp.r())};
      newSp.varianceZ() = static_cast<float>(sp.varianceZ());
      newSp.varianceR() = static_cast<float>(sp.varianceR());
    }
    std::uint32_t end = coreSpacePoints.size();
    gridSpacePointRanges.emplace_back(begin, end);
  }

  // Compute radius range. We rely on the fact the grid is storing the proxies
  // with a sorting in the radius
  const Acts::Range1D<float> rRange = [&]() -> Acts::Range1D<float> {
    float minRange = std::numeric_limits<float>::max();
    float maxRange = std::numeric_limits<float>::lowest();
    for (const Acts::SpacePointIndexRange& range : gridSpacePointRanges) {
      if (range.first == range.second) {
        continue;
      }
      auto first = coreSpacePoints[range.first];
      auto last = coreSpacePoints[range.second - 1];
      minRange = std::min(first.zr()[1], minRange);
      maxRange = std::max(last.zr()[1], maxRange);
    }
    return {minRange, maxRange};
  }();

  Acts::DoubletSeedFinder::Config bottomDoubletFinderConfig;
  bottomDoubletFinderConfig.spacePointsSortedByRadius = true;
  bottomDoubletFinderConfig.candidateDirection = Acts::Direction::Backward();
  bottomDoubletFinderConfig.deltaRMin = std::isnan(m_cfg.deltaRMinBottom)
                                            ? m_cfg.deltaRMin
                                            : m_cfg.deltaRMinBottom;
  bottomDoubletFinderConfig.deltaRMax = std::isnan(m_cfg.deltaRMaxBottom)
                                            ? m_cfg.deltaRMax
                                            : m_cfg.deltaRMaxBottom;
  bottomDoubletFinderConfig.deltaZMin = m_cfg.deltaZMin;
  bottomDoubletFinderConfig.deltaZMax = m_cfg.deltaZMax;
  bottomDoubletFinderConfig.impactMax = m_cfg.impactMax;
  bottomDoubletFinderConfig.interactionPointCut = m_cfg.interactionPointCut;
  bottomDoubletFinderConfig.collisionRegionMin = m_cfg.collisionRegionMin;
  bottomDoubletFinderConfig.collisionRegionMax = m_cfg.collisionRegionMax;
  bottomDoubletFinderConfig.cotThetaMax = m_cfg.cotThetaMax;
  bottomDoubletFinderConfig.minPt = m_cfg.minPt;
  bottomDoubletFinderConfig.helixCutTolerance = m_cfg.helixCutTolerance;
  // Vertex-z constraint takes the single experimentCuts slot when enabled;
  // otherwise fall back to the (optional) ITk fast-tracking cuts. The top
  // doublet config below is copied from this one, so the delegate is shared.
  if (!m_cfg.inputVertices.empty()) {
    bottomDoubletFinderConfig.experimentCuts.connect<&VertexZCuts::operator()>(
        &vertexZCuts);
  } else if (m_cfg.useExtraCuts) {
    bottomDoubletFinderConfig.experimentCuts.connect<itkFastTrackingCuts>();
  }
  auto bottomDoubletFinder =
      Acts::DoubletSeedFinder::create(Acts::DoubletSeedFinder::DerivedConfig(
          bottomDoubletFinderConfig, m_cfg.bFieldInZ));

  Acts::DoubletSeedFinder::Config topDoubletFinderConfig =
      bottomDoubletFinderConfig;
  topDoubletFinderConfig.candidateDirection = Acts::Direction::Forward();
  topDoubletFinderConfig.deltaRMin =
      std::isnan(m_cfg.deltaRMinTop) ? m_cfg.deltaRMin : m_cfg.deltaRMinTop;
  topDoubletFinderConfig.deltaRMax =
      std::isnan(m_cfg.deltaRMaxTop) ? m_cfg.deltaRMax : m_cfg.deltaRMaxTop;
  auto topDoubletFinder =
      Acts::DoubletSeedFinder::create(Acts::DoubletSeedFinder::DerivedConfig(
          topDoubletFinderConfig, m_cfg.bFieldInZ));

  Acts::TripletSeedFinder::Config tripletFinderConfig;
  tripletFinderConfig.useStripInfo = false;
  tripletFinderConfig.sortedByCotTheta = true;
  tripletFinderConfig.minPt = m_cfg.minPt;
  tripletFinderConfig.sigmaScattering = m_cfg.sigmaScattering;
  tripletFinderConfig.radLengthPerSeed = m_cfg.radLengthPerSeed;
  tripletFinderConfig.impactMax = m_cfg.impactMax;
  tripletFinderConfig.helixCutTolerance = m_cfg.helixCutTolerance;
  tripletFinderConfig.toleranceParam = m_cfg.toleranceParam;
  auto tripletFinder =
      Acts::TripletSeedFinder::create(Acts::TripletSeedFinder::DerivedConfig(
          tripletFinderConfig, m_cfg.bFieldInZ));

  // variable middle SP radial region of interest
  Acts::Range1D<float> rMiddleSpRange = {
      std::floor(rRange.min() / 2) * 2 + m_cfg.deltaRMiddleMinSPRange,
      std::floor(rRange.max() / 2) * 2 - m_cfg.deltaRMiddleMaxSPRange};

  // run the seeding
  Acts::BroadTripletSeedFilter::State filterState;
  Acts::BroadTripletSeedFilter::Cache filterCache;
  Acts::BroadTripletSeedFilter seedFilter(m_filterConfig, filterState,
                                          filterCache, *m_filterLogger);
  static thread_local Acts::TripletSeeder::Cache cache;

  std::vector<Acts::SpacePointContainer::ConstRange> bottomSpRanges;
  std::optional<Acts::SpacePointContainer::ConstRange> middleSpRange;
  std::vector<Acts::SpacePointContainer::ConstRange> topSpRanges;

  Acts::SeedContainer seeds;
  seeds.assignSpacePointContainer(spacePoints);

  for (const auto [bottom, middle, top] : grid.binnedGroup()) {
    ACTS_VERBOSE("Process middle " << middle);

    bottomSpRanges.clear();
    for (const auto b : bottom) {
      bottomSpRanges.push_back(
          coreSpacePoints.range(gridSpacePointRanges.at(b)).asConst());
    }
    middleSpRange =
        coreSpacePoints.range(gridSpacePointRanges.at(middle)).asConst();
    topSpRanges.clear();
    for (const auto t : top) {
      topSpRanges.push_back(
          coreSpacePoints.range(gridSpacePointRanges.at(t)).asConst());
    }

    if (middleSpRange->empty()) {
      ACTS_DEBUG("No middle space points in this group, skipping");
      continue;
    }

    // we compute this here since all middle space point candidates belong to
    // the same z-bin
    Acts::ConstSpacePointProxy firstMiddleSp = middleSpRange->front();
    std::pair<float, float> radiusRangeForMiddle =
        retrieveRadiusRangeForMiddle(firstMiddleSp, rMiddleSpRange);
    ACTS_VERBOSE("Validity range (radius) for the middle space point is ["
                 << radiusRangeForMiddle.first << ", "
                 << radiusRangeForMiddle.second << "]");

    m_seedFinder->createSeedsFromGroups(
        cache, *bottomDoubletFinder, *topDoubletFinder, *tripletFinder,
        seedFilter, coreSpacePoints, bottomSpRanges, *middleSpRange,
        topSpRanges, radiusRangeForMiddle, seeds);
  }

  ACTS_DEBUG("Created " << seeds.size() << " track seeds from "
                        << spacePoints.size() << " space points");

  // Report how strongly the vertex-z constraint fired this event.
  if (!vertexZCuts.windows.empty()) {
    ACTS_VERBOSE("Vertex-z seeding constraint: rejected "
                 << vertexZCuts.nRejected << " / " << vertexZCuts.nTested
                 << " doublet candidates; " << seeds.size()
                 << " seeds produced");
  }

  // update seed space point indices to original space point container
  for (auto seed : seeds) {
    for (auto& spIndex : seed.spacePointIndices()) {
      spIndex = coreSpacePoints.at(spIndex).copiedFromIndex();
    }
  }

  m_outputSeeds(ctx, std::move(seeds));
  return ProcessCode::SUCCESS;
}

std::pair<float, float>
GridTripletSeedingAlgorithm::retrieveRadiusRangeForMiddle(
    const Acts::ConstSpacePointProxy& spM,
    const Acts::Range1D<float>& rMiddleSpRange) const {
  if (m_cfg.useVariableMiddleSPRange) {
    return {rMiddleSpRange.min(), rMiddleSpRange.max()};
  }
  if (m_cfg.rRangeMiddleSP.empty()) {
    return {m_cfg.rMinMiddle, m_cfg.rMaxMiddle};
  }

  // get zBin position of the middle SP
  auto pVal = std::ranges::lower_bound(m_cfg.zBinEdges, spM.zr()[0]);
  std::size_t zBin = std::distance(m_cfg.zBinEdges.begin(), pVal);
  // protects against zM at the limit of zBinEdges
  zBin == 0 ? zBin : --zBin;
  return {m_cfg.rRangeMiddleSP[zBin][0], m_cfg.rRangeMiddleSP[zBin][1]};
}

}  // namespace ActsExamples
