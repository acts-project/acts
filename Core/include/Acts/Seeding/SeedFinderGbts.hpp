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

/// Tuple template used to carry the spacepoint components
using SPContainerComponentsType =
    std::tuple<SpacePointContainer2, SpacePointColumnProxy<int, true>,
               SpacePointColumnProxy<float, true>,
               SpacePointColumnProxy<float, true>>;

/// Seed finder implementing the GBTs seeding workflow.
class SeedFinderGbts {
 public:
  /// Constructor.
  /// @param config Configuration for the seed finder
  /// @param gbtsGeo GBTs geometry
  /// @param layerGeometry Layer geometry information
  /// @param logger Logging instance
  SeedFinderGbts(const SeedFinderGbtsConfig config,
                 std::unique_ptr<GbtsGeometry> gbtsGeo,
                 const std::vector<TrigInDetSiLayer>* layerGeometry,
                 std::unique_ptr<const Acts::Logger> logger =
                     Acts::getDefaultLogger("Finder",
                                            Acts::Logging::Level::INFO));

  /// Seed metadata produced by the GBTs algorithm.
  struct seedProperties {
    /// Constructor.
    /// @param quality Seed quality score
    /// @param clone Clone flag
    /// @param sps Spacepoint indices
    seedProperties(float quality, int clone, std::vector<unsigned int> sps)
        : seedQuality(quality), isClone(clone), spacepoints(std::move(sps)) {}

    /// Seed quality score.
    float seedQuality{};
    /// Clone flag.
    int isClone{};
    /// Spacepoint indices.
    std::vector<unsigned int> spacepoints{};

    /// Comparison operator.
    /// @param o Other seed properties to compare
    /// @return True if this is less than other
    bool operator<(seedProperties const& o) const {
      return std::tie(seedQuality, isClone, spacepoints) <
             std::tie(o.seedQuality, o.isClone, o.spacepoints);
    }
  };

  /// Create seeds from spacepoints in a region of interest.
  /// @param roi Region of interest descriptor
  /// @param SpContainerComponents Spacepoint container components
  /// @param max_layers Maximum number of layers
  /// @return Container with generated seeds
  SeedContainer2 createSeeds(
      const RoiDescriptor& roi,
      const SPContainerComponentsType& SpContainerComponents,
      int max_layers) const;

  /// Create graph nodes from spacepoints.
  /// @param container Spacepoint container components
  /// @param MaxLayers Maximum number of layers
  /// @return Vector of node vectors organized by layer
  std::vector<std::vector<GbtsNode>> createNodes(
      const SPContainerComponentsType& container, int MaxLayers) const;

  /// Parse machine learning lookup table from file.
  /// @param lutInputFile Path to the lookup table input file
  /// @return Parsed machine learning lookup table
  GbtsMLLookupTable parseGbtsMLLookupTable(const std::string& lutInputFile);

  /// Build doublet graph from nodes.
  /// @param roi Region of interest descriptor
  /// @param storage Data storage containing nodes
  /// @param edgeStorage Storage for generated edges
  /// @return Pair of edge count and maximum level
  std::pair<int, int> buildTheGraph(
      const RoiDescriptor& roi, const std::unique_ptr<GbtsDataStorage>& storage,
      std::vector<GbtsEdge>& edgeStorage) const;

  /// Run connected component analysis on the graph.
  /// @param nEdges Number of edges in the graph
  /// @param edgeStorage Storage containing graph edges
  /// @return Number of connected components found
  int runCCA(int nEdges, std::vector<GbtsEdge>& edgeStorage) const;

  /// Extract seed candidates from the graph.
  /// @param maxLevel Maximum level in the graph
  /// @param nEdges Number of edges
  /// @param nHits Number of hits
  /// @param edgeStorage Storage containing edges
  /// @param vSeedCandidates Output vector for seed candidates
  void extractSeedsFromTheGraph(
      int maxLevel, int nEdges, int nHits, std::vector<GbtsEdge>& edgeStorage,
      std::vector<seedProperties>& vSeedCandidates) const;

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
