// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/Seeding/GbtsConfig.hpp"
#include "Acts/Seeding/GbtsDataStorage.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/TrackFinding/RoiDescriptor.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace Acts::Experimental {

/// Tuple template used to carry the spacepoint components
using SPContainerComponentsType =
    std::tuple<SpacePointContainer2, SpacePointColumnProxy<std::uint32_t, true>,
               SpacePointColumnProxy<float, true>,
               SpacePointColumnProxy<float, true>>;

/// Seed finder implementing the GBTs seeding workflow.
class SeedFinderGbts {
 public:
  /// Seed metadata produced by the GBTs algorithm.
  struct SeedProperties {
    /// Constructor.
    /// @param quality Seed quality score
    /// @param clone Clone flag
    /// @param sps Spacepoint indices
    SeedProperties(float quality, std::int32_t clone,
                   std::vector<std::uint32_t> sps)
        : seedQuality(quality), isClone(clone), spacePoints(std::move(sps)) {}

    /// Seed quality score.
    float seedQuality{};
    /// Clone flag.
    std::int32_t isClone{};
    /// Spacepoint indices.
    std::vector<std::uint32_t> spacePoints;

    /// Comparison operator.
    /// @param o Other seed properties to compare
    /// @return True if this is less than other
    auto operator<=>(const SeedProperties& o) const = default;
  };

  /// Constructor.
  /// @param config Configuration for the seed finder
  /// @param gbtsGeo GBTs geometry
  /// @param layerGeometry Layer geometry information
  /// @param logger Logging instance
  SeedFinderGbts(const GbtsConfig& config,
                 std::unique_ptr<GbtsGeometry> gbtsGeo,
                 const std::vector<TrigInDetSiLayer>& layerGeometry,
                 std::unique_ptr<const Acts::Logger> logger =
                     Acts::getDefaultLogger("Finder",
                                            Acts::Logging::Level::INFO));

  /// Create seeds from spacepoints in a region of interest.
  /// @param roi Region of interest descriptor
  /// @param SpContainerComponents Spacepoint container components
  /// @param maxLayers Maximum number of layers
  /// @return Container with generated seeds
  SeedContainer2 createSeeds(
      const RoiDescriptor& roi,
      const SPContainerComponentsType& SpContainerComponents,
      std::uint32_t maxLayers) const;

  /// Create graph nodes from spacepoints.
  /// @param container Spacepoint container components
  /// @param maxLayers Maximum number of layers
  /// @return Vector of node vectors organized by layer
  std::vector<std::vector<GbtsNode>> createNodes(
      const SPContainerComponentsType& container,
      std::uint32_t maxLayers) const;

  /// Parse machine learning lookup table from file.
  /// @param lutInputFile Path to the lookup table input file
  /// @return Parsed machine learning lookup table
  GbtsMLLookupTable parseGbtsMLLookupTable(const std::string& lutInputFile);

  /// Build doublet graph from nodes.
  /// @param roi Region of interest descriptor
  /// @param storage Data storage containing nodes
  /// @param edgeStorage Storage for generated edges
  /// @return Pair of edge count and maximum level
  std::pair<std::int32_t, std::int32_t> buildTheGraph(
      const RoiDescriptor& roi, const std::unique_ptr<GbtsDataStorage>& storage,
      std::vector<GbtsEdge>& edgeStorage) const;

  /// Run connected component analysis on the graph.
  /// @param nEdges Number of edges in the graph
  /// @param edgeStorage Storage containing graph edges
  /// @return Number of connected components found
  std::int32_t runCCA(std::uint32_t nEdges,
                      std::vector<GbtsEdge>& edgeStorage) const;

  /// Extract seed candidates from the graph.
  /// @param maxLevel Maximum level in the graph
  /// @param nEdges Number of edges
  /// @param nHits Number of hits
  /// @param edgeStorage Storage containing edges
  /// @param vSeedCandidates Output vector for seed candidates
  void extractSeedsFromTheGraph(
      std::uint32_t maxLevel, std::uint32_t nEdges, std::int32_t nHits,
      std::vector<GbtsEdge>& edgeStorage,
      std::vector<SeedProperties>& vSeedCandidates) const;

 private:
  GbtsConfig m_cfg{};

  const std::shared_ptr<const GbtsGeometry> m_geo;

  const std::vector<TrigInDetSiLayer>* m_layerGeometry{};

  GbtsMLLookupTable m_mlLut;

  std::unique_ptr<const Acts::Logger> m_logger =
      Acts::getDefaultLogger("Finder", Acts::Logging::Level::INFO);

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace Acts::Experimental
