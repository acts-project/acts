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
    std::tuple<SpacePointContainer2, SpacePointColumnProxy<std::uint32_t, true>,
               SpacePointColumnProxy<float, true>,
               SpacePointColumnProxy<float, true>>;

class SeedFinderGbts {
 public:
  SeedFinderGbts(const SeedFinderGbtsConfig config,
                 std::unique_ptr<GbtsGeometry> gbtsGeo,
                 const std::vector<TrigInDetSiLayer>* layerGeometry,
                 std::unique_ptr<const Acts::Logger> logger =
                     Acts::getDefaultLogger("Finder",
                                            Acts::Logging::Level::INFO));

  struct SeedProperties {
    SeedProperties(float quality, std::int32_t clone,
                   std::vector<std::uint32_t> sps)
        : seedQuality(quality), isClone(clone), spacepoints(std::move(sps)) {}

    float seedQuality{};
    std::int32_t isClone{};
    std::vector<std::uint32_t> spacepoints{};

    bool operator<(SeedProperties const& o) const {
      return std::tie(seedQuality, isClone, spacepoints) <
             std::tie(o.seedQuality, o.isClone, o.spacepoints);
    }
  };

  SeedContainer2 createSeeds(
      const RoiDescriptor& roi,
      const SPContainerComponentsType& SpContainerComponents,
      std::int32_t maxLayers) const;

  std::vector<std::vector<GbtsNode>> createNodes(
      const SPContainerComponentsType& container, std::int32_t maxLayers) const;

  GbtsMLLookupTable parseGbtsMLLookupTable(const std::string& lutInputFile);

  std::pair<std::int32_t, std::int32_t> buildTheGraph(
      const RoiDescriptor& roi, const std::unique_ptr<GbtsDataStorage>& storage,
      std::vector<GbtsEdge>& edgeStorage) const;

  std::int32_t runCCA(std::int32_t nEdges,
                      std::vector<GbtsEdge>& edgeStorage) const;

  void extractSeedsFromTheGraph(
      std::int32_t maxLevel, std::int32_t nEdges, std::int32_t nHits,
      std::vector<GbtsEdge>& edgeStorage,
      std::vector<SeedProperties>& vSeedCandidates) const;

 private:
  SeedFinderGbtsConfig m_config;

  const std::shared_ptr<const GbtsGeometry> m_geo;

  const std::vector<TrigInDetSiLayer>* m_layerGeometry;

  GbtsMLLookupTable m_mlLut{};

  std::unique_ptr<const Acts::Logger> m_logger =
      Acts::getDefaultLogger("Finder", Acts::Logging::Level::INFO);

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace Acts::Experimental
