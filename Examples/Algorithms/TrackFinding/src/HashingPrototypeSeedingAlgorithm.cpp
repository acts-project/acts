// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/HashingPrototypeSeedingAlgorithm.hpp"

#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding2/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding2/DoubletSeedFinder.hpp"
#include "Acts/Seeding2/TripletSeedFinder.hpp"
#include "Acts/Utilities/AngleHelpers.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"

#include <cmath>
#include <csignal>
#include <cstddef>
#include <stdexcept>

#include <annoy/annoylib.h>
#include <annoy/kissrandom.h>

namespace ActsExamples {

namespace {

static inline bool itkFastTrackingCuts(
    const Acts::ConstSpacePointProxy2& /*middle*/,
    const Acts::ConstSpacePointProxy2& other, float cotTheta,
    bool isBottomCandidate) {
  static float rMin = 45;
  static float cotThetaMax = 1.5;

  if (isBottomCandidate && other.zr()[1] < rMin &&
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

  // Remove space points beyond eta=4 if their z is larger than the max seed
  // z0 (150.)
  float cotTheta = 27.2899;  // corresponds to eta=4
  if ((zabs - 150.) > cotTheta * r) {
    return false;
  }
  return true;
}

using AnnoyMetric = Annoy::AngularEuclidean;
using AnnoyModel =
    Annoy::AnnoyIndex<std::uint32_t, float, AnnoyMetric, Annoy::Kiss32Random,
                      Annoy::AnnoyIndexSingleThreadedBuildPolicy>;

AnnoyModel createModel(std::uint32_t f, std::uint32_t annoySeed) {
  AnnoyModel model(f);

  model.set_seed(annoySeed);

  return model;
}

void trainModel(AnnoyModel& model, const std::uint32_t index, const float phi,
                const float eta) {
  std::array<float, 2> vec{};
  if (model.get_f() >= 1) {
    vec[0] = phi;
  }
  if (model.get_f() >= 2) {
    vec[1] = eta;
  }

  model.add_item(index, vec.data());
}

void buildModel(AnnoyModel& model) {
  const std::uint32_t nTrees = 2 * model.get_f();
  model.build(nTrees);
}

std::vector<std::vector<Acts::SpacePointIndex2>> computeSpacePointsBuckets(
    const AnnoyModel& annoyModel, const Acts::SpacePointContainer2& spacePoints,
    const std::size_t bucketSize, const std::size_t zBins,
    const std::size_t phiBins, const double layerRMin, const double layerRMax,
    const double layerZMin, const double layerZMax) {
  std::vector<std::set<Acts::SpacePointIndex2>> resultSets;

  std::size_t nBins = 0;
  if (zBins > 0) {
    nBins = zBins;
  } else if (phiBins > 0) {
    nBins = phiBins;
  } else {
    throw std::runtime_error("No bins defined");
  }
  resultSets.resize(nBins);

  const auto layerSelection = [&layerRMin, &layerRMax, &layerZMin, &layerZMax](
                                  const float r2, const float z) {
    const bool isInside =
        (r2 >= layerRMin * layerRMin && r2 <= layerRMax * layerRMax) &&
        (z >= layerZMin && z <= layerZMax);
    return isInside;
  };

  const auto getBinIndexZ = [&zBins, &layerZMin, &layerZMax](const float z) {
    const float binSize = (layerZMax - layerZMin) / zBins;
    return static_cast<int>((z - layerZMin + 0.5f * binSize) / binSize);
  };

  const auto getBinIndexPhi = [&phiBins](const float phi) {
    const float binSize = 2 * std::numbers::pi / phiBins;
    return static_cast<int>((phi + std::numbers::pi) / binSize);
  };

  const auto getBinIndex = [&zBins, &phiBins, &getBinIndexZ, &getBinIndexPhi](
                               const float z, const float phi) -> int {
    if (zBins > 0) {
      return getBinIndexZ(z);
    } else if (phiBins > 0) {
      return getBinIndexPhi(phi);
    } else {
      throw std::runtime_error("No bins defined");
    }
  };

  for (const auto spacePoint : spacePoints) {
    const float z = spacePoint.zr()[0];
    const float r = spacePoint.zr()[1];
    const float r2 = r * r;
    const float phi = spacePoint.phi();

    if (!layerSelection(r2, z)) {
      continue;
    }

    const int binIndex = getBinIndex(z, phi);
    if (binIndex < 0 || static_cast<std::uint32_t>(binIndex) >= nBins) {
      throw std::runtime_error("binIndex outside of bins covering");
    }

    /// Get the `bucketSize` closest spacePoints
    std::vector<std::uint32_t> neighborSpacePointIndices;
    annoyModel.get_nns_by_item(spacePoint.index(), bucketSize, -1,
                               &neighborSpacePointIndices, nullptr);
    resultSets[binIndex].insert(spacePoint.index());
    for (const auto& neighborSpacePointIndex : neighborSpacePointIndices) {
      resultSets[binIndex].insert(neighborSpacePointIndex);
    }
  }

  std::vector<std::vector<Acts::SpacePointIndex2>> result;
  result.reserve(resultSets.size());
  for (const auto& spSet : resultSets) {
    if (!spSet.empty()) {
      result.emplace_back(spSet.begin(), spSet.end());
    }
  }

  const std::size_t nBuckets = result.size();
  const std::size_t nSpacePoints = spacePoints.size();

  // Check if the number of buckets is greater than the number of space points
  if (nBuckets > nSpacePoints) {
    throw std::runtime_error("More buckets than the number of space points");
  }

  return result;
}

struct SeedComparison {
  bool operator()(const Acts::Seed<SimSpacePoint>& seed1,
                  const Acts::Seed<SimSpacePoint>& seed2) const {
    const auto& sp1 = seed1.sp();
    const auto& sp2 = seed2.sp();

    for (std::size_t i = 0; i < sp1.size(); ++i) {
      if (sp1[i]->z() != sp2[i]->z()) {
        return sp1[i]->z() < sp2[i]->z();
      }
    }

    for (std::size_t i = 0; i < sp1.size(); ++i) {
      if (sp1[i]->x() != sp2[i]->x()) {
        return sp1[i]->x() < sp2[i]->x();
      }
    }

    for (std::size_t i = 0; i < sp1.size(); ++i) {
      if (sp1[i]->y() != sp2[i]->y()) {
        return sp1[i]->y() < sp2[i]->y();
      }
    }

    return false;
  }
};

}  // namespace

HashingPrototypeSeedingAlgorithm::HashingPrototypeSeedingAlgorithm(
    Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("HashingPrototypeSeedingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_outputSeeds.initialize(m_cfg.outputSeeds);
  m_outputBuckets.initialize(m_cfg.outputBuckets);

  if (m_cfg.useExtraCuts) {
    // This function will be applied to select space points during grid filling
    m_spacePointSelector.connect<itkFastTrackingSPselect>();
  }

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

  m_filterLogger = logger().cloneWithSuffix("Filter");

  m_seedFinder = Acts::TripletSeeder(logger().cloneWithSuffix("Finder"));
}

ProcessCode HashingPrototypeSeedingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  Acts::DoubletSeedFinder::Config bottomDoubletFinderConfig;
  bottomDoubletFinderConfig.spacePointsSortedByRadius = true;
  bottomDoubletFinderConfig.candidateDirection = Acts::Direction::Backward();
  bottomDoubletFinderConfig.deltaRMin = std::isnan(m_cfg.deltaRMaxBottom)
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
  if (m_cfg.useExtraCuts) {
    bottomDoubletFinderConfig.experimentCuts.connect<itkFastTrackingCuts>();
  }
  auto bottomDoubletFinder =
      Acts::DoubletSeedFinder::create(Acts::DoubletSeedFinder::DerivedConfig(
          bottomDoubletFinderConfig, m_cfg.bFieldInZ));

  Acts::DoubletSeedFinder::Config topDoubletFinderConfig =
      bottomDoubletFinderConfig;
  topDoubletFinderConfig.candidateDirection = Acts::Direction::Forward();
  topDoubletFinderConfig.deltaRMin =
      std::isnan(m_cfg.deltaRMaxTop) ? m_cfg.deltaRMin : m_cfg.deltaRMinTop;
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

  // run the seeding
  Acts::BroadTripletSeedFilter::State filterState;
  Acts::BroadTripletSeedFilter::Cache filterCache;
  Acts::BroadTripletSeedFilter seedFilter(m_filterConfig, filterState,
                                          filterCache, *m_filterLogger);
  static thread_local Acts::TripletSeeder::Cache cache;

  const SimSpacePointContainer& spacePoints = m_inputSpacePoints(ctx);

  Acts::SpacePointContainer2 coreSpacePoints(
      Acts::SpacePointColumns::SourceLinks | Acts::SpacePointColumns::PackedXY |
      Acts::SpacePointColumns::PackedZR | Acts::SpacePointColumns::VarianceZ |
      Acts::SpacePointColumns::VarianceR | Acts::SpacePointColumns::Phi);

  // create and train the hashing model
  AnnoyModel hashingModel = createModel(m_cfg.f, m_cfg.annoySeed);
  for (const auto& sp : spacePoints) {
    // check if the space point passes the selection
    if (!m_spacePointSelector(sp)) {
      continue;
    }

    auto newSp = coreSpacePoints.createSpacePoint();
    newSp.assignSourceLinks(
        std::array<Acts::SourceLink, 1>{Acts::SourceLink(&sp)});
    newSp.xy() = std::array<float, 2>{static_cast<float>(sp.x()),
                                      static_cast<float>(sp.y())};
    newSp.zr() = std::array<float, 2>{static_cast<float>(sp.z()),
                                      static_cast<float>(sp.r())};
    newSp.varianceZ() = static_cast<float>(sp.varianceZ());
    newSp.varianceR() = static_cast<float>(sp.varianceR());

    const float phi = std::atan2(newSp.xy()[1], newSp.xy()[0]);
    const float eta = Acts::AngleHelpers::etaFromTheta(
        std::atan2(newSp.zr()[1], newSp.zr()[0]));
    trainModel(hashingModel, newSp.index(), phi, eta);

    newSp.phi() = phi;
  }
  buildModel(hashingModel);

  // create buckets based on hashing model
  std::vector<std::vector<Acts::SpacePointIndex2>> buckets =
      computeSpacePointsBuckets(hashingModel, coreSpacePoints, m_cfg.bucketSize,
                                m_cfg.zBins, m_cfg.phiBins, m_cfg.layerRMin,
                                m_cfg.layerRMax, m_cfg.layerZMin,
                                m_cfg.layerZMax);

  ACTS_DEBUG("Created " << buckets.size() << " buckets  from "
                        << coreSpacePoints.size() << " space points");

  // sort buckets by radius
  for (auto& bucket : buckets) {
    std::ranges::sort(bucket, {}, [&](const Acts::SpacePointIndex2& spIndex) {
      return coreSpacePoints.at(spIndex).zr()[1];
    });
  }

  std::set<Acts::Seed<SimSpacePoint>, SeedComparison> uniqueSeeds;

  for (const auto& bucket : buckets) {
    auto bucketSpacePoints = coreSpacePoints.subset(bucket).asConst();

    Acts::SeedContainer2 seeds;

    for (auto middleSp : bucketSpacePoints) {
      m_seedFinder->createSeedsFromGroup(
          cache, *bottomDoubletFinder, *topDoubletFinder, *tripletFinder,
          seedFilter, coreSpacePoints, bucketSpacePoints, middleSp,
          bucketSpacePoints, seeds);
    }

    ACTS_DEBUG("Created " << seeds.size() << " track seeds from "
                          << bucketSpacePoints.size() << " space points");

    for (const auto& seed : seeds) {
      auto sps = seed.spacePointIndices();

      Acts::Seed<SimSpacePoint> newSeed(*coreSpacePoints.at(sps[0])
                                             .sourceLinks()[0]
                                             .get<const SimSpacePoint*>(),
                                        *coreSpacePoints.at(sps[1])
                                             .sourceLinks()[0]
                                             .get<const SimSpacePoint*>(),
                                        *coreSpacePoints.at(sps[2])
                                             .sourceLinks()[0]
                                             .get<const SimSpacePoint*>());
      newSeed.setVertexZ(seed.vertexZ());
      newSeed.setQuality(seed.quality());

      uniqueSeeds.insert(newSeed);
    }
  }

  ACTS_DEBUG("Total unique seeds created: " << uniqueSeeds.size());

  std::vector<SimSpacePointContainer> outputBuckets;
  for (const auto& bucket : buckets) {
    SimSpacePointContainer spContainer;
    for (const auto& spIndex : bucket) {
      const auto& sp = *coreSpacePoints.at(spIndex)
                            .sourceLinks()[0]
                            .get<const SimSpacePoint*>();
      spContainer.push_back(sp);
    }
    outputBuckets.push_back(std::move(spContainer));
  }

  m_outputSeeds(ctx, SimSeedContainer(uniqueSeeds.begin(), uniqueSeeds.end()));
  m_outputBuckets(ctx, std::move(outputBuckets));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
