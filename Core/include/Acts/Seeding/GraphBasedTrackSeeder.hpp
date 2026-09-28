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
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"
#include "Acts/Seeding/GbtsNodeStorage.hpp"
#include "Acts/Seeding/GbtsRoiDescriptor.hpp"
#include "Acts/Seeding/GbtsTrackingFilter.hpp"
#include "Acts/Seeding/detail/GbtsGraphTypes.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace Acts::Experimental {

class GbtsGraphBuilder;

/// Seed finder implementing the GBTS seeding workflow.
class GraphBasedTrackSeeder {
 public:
  /// Configuration struct for the GBTS seeding algorithm.
  struct Config {
    /// Accepted tau range per cluster width bin, needed by the cluster width
    /// cuts and ignored without them.
    detail::GbtsTauLookupTable tauLookupTable;
    /// Enable the cluster width cuts: wide endcap rejection and tau narrowing.
    bool useClusterWidthCuts = false;
    /// Maximum number of phi slices.
    float nMaxPhiSlice = 53;  // used to calculate phi slices

    // Seed extraction options
    //
    // @note The chain length a candidate must reach, and the `addTriplets`
    //       relaxation of it, are `GbtsGraphBuilder::Config`: the graph applies
    //       them when it picks the chain heads, and extraction here has to
    //       agree.

    /// Smallest seed size that is split into drop-out candidates.
    std::uint32_t minSplitSeedSize = 4;
    /// Largest seed size that is split.
    std::uint32_t maxSplitSeedSize = 5;
    /// Minimum eta for edge masking.
    float edgeMaskMinEta = 1.5;
    /// Threshold for hit sharing between seeds.
    float hitShareThreshold = 0.49;
    /// Max seed eta value considered for splitting.
    float maxSeedSplitEta = 0.6;
    /// Max allowed curvature for seed self consistency check.
    float maxInvRadDiff = 0.7e-2 / UnitConstants::m;
    // GbtsNodeStorage options
    /// Maximum endcap cluster width.
    float maxEndcapClusterWidth = 0.35f * Acts::UnitConstants::mm;
    /// Half-length in local y of a pixel module, against which the distance of
    /// a cluster to the module edge is measured.
    float moduleHalfLengthY = 10.0f * Acts::UnitConstants::mm;
    /// Distance to the module edge below which a cluster may be shortened,
    /// which switches to the tau lookup table's near-edge bounds.
    float moduleEdgeTolerance = 0.3f * Acts::UnitConstants::mm;
    /// Cluster width covered by one bin of the tau lookup table.
    float tauLutBinWidth = 0.05f * Acts::UnitConstants::mm;
    /// Multiples of the phi slice width duplicated either side of the
    /// wrap-around, so a sliding window never has to wrap.
    float phiIndexMargin = 1.5f;
    /// Buckets used to sort an eta bin by phi, at most
    /// `GbtsNodeStorage::kMaxPhiSortBuckets`.
    std::uint32_t phiSortBuckets = 31;
  };

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
  struct Options {
    /// Magnetic field in z
    /// units of GeV/(e*mm).
    float bFieldInZ{};
  };

  /// @param config Configuration for the seed finder
  /// @param geometry GBTS geometry
  /// @param logger Logging instance
  GraphBasedTrackSeeder(const DerivedConfig& config,
                        std::shared_ptr<GbtsGeometry> geometry,
                        std::unique_ptr<const Acts::Logger> logger =
                            Acts::getDefaultLogger("Finder",
                                                   Acts::Logging::Level::INFO));

  /// Create an empty node storage matching this seeder's configuration.
  ///
  /// Fill it through GbtsNodeStorage::insert, then call
  /// GbtsNodeStorage::finalize before handing it to createSeeds.
  /// @return An empty node storage
  GbtsNodeStorage makeNodeStorage() const;

  /// Create seeds from an ACTS space point container in a region of interest.
  ///
  /// Convenience wrapper that builds and finalizes the node storage itself. The
  /// container must carry the `gbtsLayerIndex`, `clusterWidth` and
  /// `localPositionY` columns.
  /// @param spacePoints Space point container
  /// @param roi Region of interest descriptor
  /// @param graphBuilder Doublet graph builder, which also carries the chain
  ///              selection that this seeder's extraction agrees with
  /// @param filter Tracking filter to be applied
  /// @param options Event based options such as magnetic field strength
  /// @param outputSeeds Container with generated seeds
  void createSeeds(const SpacePointContainer& spacePoints,
                   const GbtsRoiDescriptor& roi,
                   const GbtsGraphBuilder& graphBuilder,
                   const GbtsTrackingFilter& filter, const Options& options,
                   SeedContainer& outputSeeds) const;

  /// Create seeds from a finalized node storage in a region of interest.
  /// @param nodeStorage Finalized graph node storage
  /// @param roi Region of interest descriptor
  /// @param graphBuilder Doublet graph builder, which also carries the chain
  ///              selection that this seeder's extraction agrees with
  /// @param filter Tracking filter to be applied
  /// @param options Event based options such as magnetic field strength
  /// @param outputSeeds Container with generated seeds
  void createSeeds(GbtsNodeStorage& nodeStorage, const GbtsRoiDescriptor& roi,
                   const GbtsGraphBuilder& graphBuilder,
                   const GbtsTrackingFilter& filter, const Options& options,
                   SeedContainer& outputSeeds) const;

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
  };

  DerivedConfig m_cfg;

  std::shared_ptr<const GbtsGeometry> m_geometry;

  std::unique_ptr<const Acts::Logger> m_logger =
      Acts::getDefaultLogger("Finder", Acts::Logging::Level::INFO);

  const Acts::Logger& logger() const { return *m_logger; }

  /// Extract seed candidates from the graph.
  /// @param nodeStorage Storage containing the graph nodes
  /// @param edgeStorage Storage containing edges
  /// @param vOutputSeeds Output vector for seed candidates
  /// @param filter Tracking filter to be applied
  /// @param vChainHeads Chain heads the graph selected
  /// @param graphBuilder Doublet graph builder, read for its chain selection
  void extractSeedsFromTheGraph(const GbtsNodeStorage& nodeStorage,
                                std::vector<detail::GbtsEdge>& edgeStorage,
                                std::vector<OutputSeedProperties>& vOutputSeeds,
                                const GbtsTrackingFilter& filter,
                                std::vector<detail::GbtsEdge*>& vChainHeads,
                                const GbtsGraphBuilder& graphBuilder) const;

  /// Estimate the inverse radius of the circle through three nodes.
  /// @param nodeView View of the node positions and layers
  /// @param nodes The three graph nodes, innermost first
  /// @return The estimated inverse radius
  float estimateCurvature(const detail::GbtsNodeView& nodeView,
                          const std::array<SpacePointIndex, 3>& nodes) const;
};

}  // namespace Acts::Experimental
