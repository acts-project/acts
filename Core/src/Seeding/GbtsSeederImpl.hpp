// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"
#include "Acts/Seeding/GbtsRoiDescriptor.hpp"
#include "Acts/Seeding/GbtsSeedFit.hpp"
#include "Acts/Seeding/GbtsTauLookupTable.hpp"
#include "Acts/Seeding/GraphBasedTrackSeeder.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "GbtsGeometryImpl.hpp"
#include "GbtsGraphBuilder.hpp"
#include "GbtsGraphTypes.hpp"
#include "GbtsNodeStorageImpl.hpp"
#include "GbtsTrackingFilter.hpp"

namespace Acts::detail {

struct GbtsGraph;
class GbtsGraphBuilder;

/// Seed finder implementing the GBTS seeding workflow.
class GbtsSeederImpl {
 public:
  /// Configuration struct for the GBTS seeding algorithm.
  using Config = Acts::GraphBasedTrackSeeder::Config;

  /// Derived configuration struct that contains calculated parameters based on
  /// the configuration.
  struct DerivedConfig : public Config {
    /// Construct derived configuration from base configuration.
    /// @param config Base configuration to derive from
    explicit DerivedConfig(const Config& config);

    /// Phi slice width
    float phiSliceWidth = std::numeric_limits<float>::quiet_NaN();
  };

  /// Optional inputs for variables passed in
  /// or derived during runtime.
  using Options = Acts::GraphBasedTrackSeeder::Options;

  /// @param config Configuration for the seed finder
  /// @param geometry GBTS geometry
  /// @param logger Logging instance
  GbtsSeederImpl(const DerivedConfig& config,
                 std::shared_ptr<const GbtsGeometryImpl> geometry,
                 std::unique_ptr<const Acts::Logger> logger =
                     Acts::getDefaultLogger("Finder",
                                            Acts::Logging::Level::INFO));

  /// Create an empty node storage matching this seeder's configuration.
  ///
  /// Fill it through GbtsNodeStorageImpl::insert, then call
  /// GbtsNodeStorageImpl::finalize before handing it to createSeeds.
  /// @return An empty node storage
  std::unique_ptr<GbtsNodeStorageImpl> makeNodeStorage() const;

  /// Create seeds from an ACTS space point container in a region of interest.
  ///
  /// Convenience wrapper that builds and finalizes the node storage itself. The
  /// container must carry the `gbtsLayerIndex`, `clusterWidth` and
  /// `localPositionY` columns.
  /// @param spacePoints Space point container
  /// @param roi Region of interest descriptor
  /// @param options Event based options such as magnetic field strength
  /// @param outputSeeds Container with generated seeds
  void createSeeds(const SpacePointContainer& spacePoints,
                   const GbtsRoiDescriptor& roi, const Options& options,
                   SeedContainer& outputSeeds) const;

  /// Create seeds from a finalized node storage in a region of interest.
  /// @param nodeStorage Finalized graph node storage
  /// @param roi Region of interest descriptor
  /// @param options Event based options such as magnetic field strength
  /// @param outputSeeds Container with generated seeds
  /// @param outputFits Fit of the tracking filter to each generated seed
  void createSeeds(
      GbtsNodeStorageImpl& nodeStorage, const GbtsRoiDescriptor& roi,
      const Options& options, SeedContainer& outputSeeds,
      std::vector<std::optional<GbtsSeedFit>>* outputFits = nullptr) const;

 private:
  /// candidate seed metadata produced by the GBTS algorithm.
  struct SeedCandidateProperties {
    /// @param quality Seed quality score
    /// @param clone Whether the candidate was rejected as a clone
    /// @param sps Vector of graph node indices
    /// @param splitFlag used to flag if seed needs to be split in two
    SeedCandidateProperties(float quality, bool clone,
                            std::vector<SpacePointIndex> sps, bool splitFlag)
        : seedQuality(quality),
          isClone(clone),
          nodes(std::move(sps)),
          forSeedSplitting(splitFlag) {}

    /// Seed quality score.
    float seedQuality{};
    /// Clone flag.
    bool isClone{};
    /// Graph node indices.
    std::vector<SpacePointIndex> nodes;
    /// Flag for seed splitting.
    bool forSeedSplitting{};
    /// Fit of the tracking filter.
    std::optional<GbtsSeedFit> fit;
  };

  /// Output seed metadata
  struct OutputSeedProperties {
    /// @param quality Seed quality score
    /// @param sps Vector of space point indices in the seed
    OutputSeedProperties(float quality, std::vector<std::uint32_t> sps)
        : seedQuality(quality), spacePoints(std::move(sps)) {}

    /// Quality of seed.
    float seedQuality{};
    /// Index of spacepoints in seed.
    std::vector<std::uint32_t> spacePoints;
    /// Fit of the tracking filter.
    std::optional<GbtsSeedFit> fit;
  };

  GbtsGraphBuilder m_graphBuilder;
  GbtsTrackingFilter m_filter;
  DerivedConfig m_cfg;

  std::shared_ptr<const GbtsGeometryImpl> m_geometry;

  std::unique_ptr<const Acts::Logger> m_logger =
      Acts::getDefaultLogger("Finder", Acts::Logging::Level::INFO);

  const Acts::Logger& logger() const { return *m_logger; }

  /// Extract seed candidates from the graph.
  /// @param nodeStorage Storage containing the graph nodes
  /// @param graph The graph
  /// @param vOutputSeeds Output vector for seed candidates
  /// @param filter Tracking filter to be applied
  /// @param withFits Whether the fit of the filter to each seed is kept
  void extractSeedsFromTheGraph(const GbtsNodeStorageImpl& nodeStorage,
                                const GbtsGraph& graph,
                                std::vector<OutputSeedProperties>& vOutputSeeds,
                                const GbtsTrackingFilter& filter,
                                bool withFits) const;

  /// Extract the best fitted path from every edge inwards to a root edge.
  /// @param nodeView View of the node positions and layers
  /// @param graph The graph
  /// @param filter Tracking filter fitting the paths
  /// @return The fitted paths, each with its head edge first
  std::vector<detail::GbtsEdgeState> extractBestPathPerEdge(
      const detail::GbtsNodeView& nodeView, const GbtsGraph& graph,
      const GbtsTrackingFilter& filter) const;

  /// Estimate the inverse radius of the circle through three nodes.
  /// @param nodeView View of the node positions and layers
  /// @param nodes The three graph nodes, innermost first
  /// @return The estimated inverse radius
  float estimateCurvature(const detail::GbtsNodeView& nodeView,
                          const std::array<SpacePointIndex, 3>& nodes) const;
};

}  // namespace Acts::detail
