// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/Seeding/GbtsDataStorage.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsLutParser.hpp"
#include "Acts/Seeding/SeedFinderGbtsConfig.hpp"
#include "Acts/TrackFinding/RoiDescriptor.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace Acts::Experimental {

// defined the tuple template used to carry the spacepoint components
using SPContainerComponentsType =
    std::tuple<SpacePointContainer2, SpacePointColumnProxy<int, true>,
               SpacePointColumnProxy<float, true>,
               SpacePointColumnProxy<float, true>>;

class SeedFinderGbts {
 public:
  SeedFinderGbts(const SeedFinderGbtsConfig config, const GbtsGeometry* gbtsGeo,
                 const std::vector<TrigInDetSiLayer>* layerGeometry,
                 const GbtsLutParser* gbtsLutParser,
                 std::unique_ptr<const Acts::Logger> logger =
                     Acts::getDefaultLogger("Finder",
                                            Acts::Logging::Level::INFO));

  struct seedProperties {
    seedProperties(float quality, int clone, std::vector<unsigned int> sps)
        : seedQuality(quality), isClone(clone), spacepoints(std::move(sps)) {}

    float seedQuality{};
    int isClone{};
    std::vector<unsigned int> spacepoints{};

    bool operator<(seedProperties const& o) const {
      return std::tie(seedQuality, isClone, spacepoints) <
             std::tie(o.seedQuality, o.isClone, o.spacepoints);
    }
  };

  SeedContainer2 CreateSeeds(
      const RoiDescriptor& roi,
      const SPContainerComponentsType& SpContainerComponents,
      int max_layers) const;

  std::vector<std::vector<GbtsNode>> CreateNodes(const auto& container,
                                                 int MaxLayers) const;

  std::pair<int, int> buildTheGraph(
      const RoiDescriptor& roi, const std::unique_ptr<GbtsDataStorage>& storage,
      std::vector<GbtsEdge>& edgeStorage) const;

  int runCCA(int nEdges, std::vector<GbtsEdge>& edgeStorage) const;

  void extractSeedsFromTheGraph(
      int maxLevel, int nEdges, int nHits, std::vector<GbtsEdge>& edgeStorage,
      std::vector<seedProperties>& vSeedCandidates) const;

 private:
  SeedFinderGbtsConfig m_config;

  const GbtsGeometry* m_geo;

  const std::vector<TrigInDetSiLayer>* m_layerGeometry;

  const GbtsLutParser* m_lutParser;

  std::unique_ptr<const Acts::Logger> m_logger =
      Acts::getDefaultLogger("Finder", Acts::Logging::Level::INFO);

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace Acts::Experimental
