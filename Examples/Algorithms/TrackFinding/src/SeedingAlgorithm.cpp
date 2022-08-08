// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SpacePoint.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>

ActsExamples::SeedingAlgorithm::SeedingAlgorithm(
    ActsExamples::SeedingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("SeedingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point input collections");
  }
  for (const auto& i : m_cfg.inputSpacePoints) {
    if (i.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing proto tracks output collection");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds output collection");
  }

  if (m_cfg.gridConfig.rMax != m_cfg.seedFinderConfig.rMax and
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

  if (m_cfg.gridConfig.bFieldInZ != m_cfg.seedFinderConfig.bFieldInZ) {
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

  if (m_cfg.seedFinderConfig.zBinsCustomLooping.size() != 0) {
    // check if zBinsCustomLooping contains numbers from 1 to the total number
    // of bin in zBinEdges
    for (size_t i = 1; i != m_cfg.gridConfig.zBinEdges.size(); i++) {
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
        [](const void*, const Acts::SpacePoint& sp) -> float {
          const auto* simSPSL =
              static_cast<const ActsExamples::SimSPSourceLink*>(
                  sp.getSourceLinks()[0]);
          return simSPSL->getSimSP()->topHalfStripLength();
        });

    m_cfg.seedFinderConfig.getBottomHalfStripLength.connect(
        [](const void*, const Acts::SpacePoint& sp) -> float {
          const auto* simSPSL =
              static_cast<const ActsExamples::SimSPSourceLink*>(
                  sp.getSourceLinks()[0]);
          return simSPSL->getSimSP()->bottomHalfStripLength();
        });

    m_cfg.seedFinderConfig.getTopStripDirection.connect(
        [](const void*, const Acts::SpacePoint& sp) -> Acts::Vector3 {
          const auto* simSPSL =
              static_cast<const ActsExamples::SimSPSourceLink*>(
                  sp.getSourceLinks()[0]);
          return simSPSL->getSimSP()->topStripDirection();
        });

    m_cfg.seedFinderConfig.getBottomStripDirection.connect(
        [](const void*, const Acts::SpacePoint& sp) -> Acts::Vector3 {
          const auto* simSPSL =
              static_cast<const ActsExamples::SimSPSourceLink*>(
                  sp.getSourceLinks()[0]);
          return simSPSL->getSimSP()->bottomStripDirection();
        });

    m_cfg.seedFinderConfig.getStripCenterDistance.connect(
        [](const void*, const Acts::SpacePoint& sp) -> Acts::Vector3 {
          const auto* simSPSL =
              static_cast<const ActsExamples::SimSPSourceLink*>(
                  sp.getSourceLinks()[0]);
          return simSPSL->getSimSP()->stripCenterDistance();
        });

    m_cfg.seedFinderConfig.getTopStripCenterPosition.connect(
        [](const void*, const Acts::SpacePoint& sp) -> Acts::Vector3 {
          const auto* simSPSL =
              static_cast<const ActsExamples::SimSPSourceLink*>(
                  sp.getSourceLinks()[0]);
          return simSPSL->getSimSP()->topStripCenterPosition();
        });
  }

  m_cfg.seedFinderConfig.seedFilter =
      std::make_unique<Acts::SeedFilter>(m_cfg.seedFilterConfig);
}

ActsExamples::ProcessCode ActsExamples::SeedingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // extent used to store r range for middle spacepoint
  Acts::Extent rRangeSPExtent;

  // create iterator and wrapper container to loop over all space points
  // without copying them
  struct pointer_iterator : std::iterator_traits<Acts::SpacePoint*> {
    using value_type = std::add_pointer_t<Acts::SpacePoint>;

    value_type ptr;

    pointer_iterator(value_type i) : ptr(i){};

    void operator++(int) { ptr++; }

    void operator++() { ++ptr; }

    bool operator!=(const pointer_iterator& o) { return ptr != o.ptr; }

    value_type operator*() { return ptr; }
  };

  struct SPContainer {
    SPContainer(std::vector<Acts::SpacePoint>& vec) : spacePoints(vec) {}

    std::vector<Acts::SpacePoint>& spacePoints;
    size_t index = 0;

    void operator++() { ++index; }
    void operator++(int) { index++; }

    Acts::SpacePoint* operator*() { return &(spacePoints.at(index)); }

    pointer_iterator begin() const {
      return pointer_iterator(&*spacePoints.begin());
    }

    pointer_iterator end() const {
      return pointer_iterator(&*spacePoints.end());
    }
  };

  std::vector<Acts::SpacePoint> actsSpacePoints;
  for (const auto& isp : m_cfg.inputSpacePoints) {
    for (const auto& spacePoint :
         ctx.eventStore.get<SimSpacePointContainer>(isp)) {
      // since the event store owns the space points, their pointers should be
      // stable and we do not need to create local copies.
      Acts::SpacePoint sp(
          {spacePoint.x(), spacePoint.y(), spacePoint.z()},
          {0., 0.},  // beamspot position in x,y of the Seedfinder coordinate
                     // system, used as offset to position
          {spacePoint.varianceR(), spacePoint.varianceZ()},
          {0., 0.},  // alignment uncertainty, added to measurement uncertainty
                     // (varianceR, varianceZ)
          m_cfg.seedFinderConfig.sigmaError,
          {static_cast<const Acts::SourceLink*>(
              new SimSPSourceLink(&spacePoint))});
      actsSpacePoints.push_back(sp);

      // store x,y,z values in extent
      rRangeSPExtent.check({spacePoint.x(), spacePoint.y(), spacePoint.z()});
    }
  }
  SPContainer actsSpPtrs(actsSpacePoints);

  auto bottomBinFinder = std::make_shared<Acts::BinFinder>(
      Acts::BinFinder(m_cfg.zBinNeighborsBottom, m_cfg.numPhiNeighbors));
  auto topBinFinder = std::make_shared<Acts::BinFinder>(
      Acts::BinFinder(m_cfg.zBinNeighborsTop, m_cfg.numPhiNeighbors));
  auto grid = Acts::SpacePointGridCreator::createGrid(m_cfg.gridConfig);
  auto spacePointsGrouping = Acts::BinnedSPGroup(
      actsSpPtrs.begin(), actsSpPtrs.end(), bottomBinFinder, topBinFinder,
      std::move(grid), m_cfg.seedFinderConfig);
  auto finder = Acts::Seedfinder(m_cfg.seedFinderConfig);

  // run the seeding
  static thread_local Acts::SeedContainer seeds;
  seeds.clear();
  static thread_local decltype(finder)::State state;

  auto group = spacePointsGrouping.begin();
  auto groupEnd = spacePointsGrouping.end();
  for (; !(group == groupEnd); ++group) {
    finder.createSeedsForGroup(state, std::back_inserter(seeds), group.bottom(),
                               group.middle(), group.top(), rRangeSPExtent);
  }

  // extract proto tracks, i.e. groups of measurement indices, from tracks seeds
  size_t nSeeds = seeds.size();
  static thread_local ProtoTrackContainer protoTracks;
  protoTracks.clear();

  protoTracks.reserve(nSeeds);
  for (const auto& seed : seeds) {
    ProtoTrack& protoTrack = protoTracks.emplace_back();
    protoTrack.reserve(seed.sp.size());
    for (auto spacePointPtr : seed.sp) {
      protoTrack.push_back(static_cast<const SimSPSourceLink*>(
                               spacePointPtr->getSourceLinks()[0])
                               ->getSimSP()
                               ->measurementIndex());
    }
  }

  ACTS_DEBUG("Created " << seeds.size() << " track seeds from "
                        << actsSpacePoints.size() << " space points");

  ctx.eventStore.add(m_cfg.outputSeeds, Acts::SeedContainer{seeds});
  ctx.eventStore.add(m_cfg.outputProtoTracks, ProtoTrackContainer{protoTracks});
  return ActsExamples::ProcessCode::SUCCESS;
}
