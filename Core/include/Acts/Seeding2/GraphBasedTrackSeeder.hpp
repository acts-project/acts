// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding2/GbtsDataStorage.hpp"
#include "Acts/Seeding2/GbtsGeometry.hpp"
#include "Acts/Seeding2/GbtsRoiDescriptor.hpp"
#include "Acts/Seeding2/GbtsTrackingFilter.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace Acts::Experimental {

/// Seed finder implementing the GBTS seeding workflow.
class GraphBasedTrackSeeder {
 public:
  /// Configuration struct for the GBTS seeding algorithm.
  struct Config {
    // GbtsSeedingAlgorithm options
    /// Enable beam spot correction.
    bool beamSpotCorrection = false;

    // Path to the connector configuration file that defines the layer
    // connections
    /// Connector configuration file path.
    std::string connectorInputFile;

    /// Look-up table input file path.
    std::string lutInputFile;

    // SeedFinderGbts option
    /// Enable Large Radius Tracking mode.
    bool lrtMode = false;
    /// Use machine learning features (e.g., cluster width).
    bool useMl = false;
    /// Match seeds before creating them.
    bool matchBeforeCreate = false;
    /// Use legacy tuning parameters.
    bool useOldTunings = false;
    /// description here
    bool validateTriplets = true;
    /// description here
    bool useAdaptiveCuts = true;
    /// Tau ratio cut threshold.
    float tauRatioCut = 0.007;
    /// Tau ratio precut threshold.
    float tauRatioPrecut = 0.009f;
    /// description here
    float tauRatioCorr = 0.006;
    /// Eta bin width override (0 uses default from connection file).
    // specify non-zero to override eta bin width from connection file (default
    // 0.2 in createLinkingScheme.py)
    float etaBinWidthOverride = 0.0f;

    /// Maximum number of phi slices.
    float nMaxPhiSlice = 53;  // used to calculate phi slices
    /// Minimum transverse momentum.
    float minPt = 1.0f * UnitConstants::GeV;

    // graph building options
    /// Use eta binning from geometry structure.
    bool useEtaBinning = true;
    /// Apply RZ cuts on doublets.
    bool doubletFilterRZ = true;
    /// Maximum number of Gbts edges/doublets.
    std::uint32_t nMaxEdges = 2000000;
    /// Minimum delta radius between layers.
    float minDeltaRadius = 2.0 * Acts::UnitConstants::mm;
    /// Maximum d0 impact parameter when validating edge connection triplet
    float d0Max = 3.0 * UnitConstants::mm;

    // Seed extraction options
    /// Minimum eta for edge masking.
    float edgeMaskMinEta = 1.5;
    /// Threshold for hit sharing between seeds.
    float hitShareThreshold = 0.49;
    /// Max seed eta value considered for splitting.
    float maxSeedSplitEta = 0.6;
    /// Max allowed curvature for seed self consistency check.
    float maxInvRadDiff = 0.7e-2 / UnitConstants::m;
    // GbtsDataStorage options
    /// Maximum endcap cluster width.
    float maxEndcapClusterWidth = 0.35 * Acts::UnitConstants::mm;
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
    /// @param bFieldInZ_ the magnetic field in z
    explicit Options(float bFieldInZ_);

    /// Magnetic field in z
    /// units of GeV/(e*mm).
    float bFieldInZ{};

    /// Transverse momentum coefficient (~0.3*B/2 - assumes nominal field of
    /// 2*T).
    double ptCoeff{};
  };

  /// candidate seed metadata produced by the GBTS algorithm.
  struct SeedCandidateProperties {
    /// @param quality Seed quality score
    /// @param clone Clone flag
    /// @param sps Vector of pointers to actual space points
    /// @param splitFlag used to flag if seed needs to be split in two
    SeedCandidateProperties(float quality, std::int32_t clone,
                            std::vector<const GbtsNode*> sps,
                            std::uint32_t splitFlag)
        : seedQuality(quality),
          isClone(clone),
          spacePoints(std::move(sps)),
          forSeedSplitting(splitFlag) {}

    /// Seed quality score.
    float seedQuality{};
    /// Clone flag.
    std::int32_t isClone{};
    /// Space point indices.
    std::vector<const GbtsNode*> spacePoints;
    /// Flag for seed splitting.
    std::uint32_t forSeedSplitting{};
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

  /// Sliding window in phi used to define range used for edge creation
  struct SlidingWindow {
    /// sliding window position
    std::uint32_t firstIt{};
    /// window half-width;
    float deltaPhi{};
    /// active or not
    bool hasNodes{};
    /// associated eta bin
    const GbtsEtaBin* bin{};
  };

  /// @param config Configuration for the seed finder
  /// @param geometry GBTS geometry
  /// @param logger Logging instance
  GraphBasedTrackSeeder(const DerivedConfig& config,
                        std::shared_ptr<GbtsGeometry> geometry,
                        std::unique_ptr<const Acts::Logger> logger =
                            Acts::getDefaultLogger("Finder",
                                                   Acts::Logging::Level::INFO));

  /// Create seeds from space points in a region of interest.
  /// @param spacePoints Space point container
  /// @param roi Region of interest descriptor
  /// @param maxLayers Maximum number of layers
  /// @param filter Tracking filter to be applied
  /// @param options Event based options such as magnetic field strength
  /// @param outputSeeds Container with generated seeds
  void createSeeds(const SpacePointContainer2& spacePoints,
                   const GbtsRoiDescriptor& roi, std::uint32_t maxLayers,
                   const GbtsTrackingFilter& filter, const Options& options,
                   SeedContainer2& outputSeeds) const;

 private:
  DerivedConfig m_cfg;

  std::shared_ptr<const GbtsGeometry> m_geometry;

  GbtsMlLookupTable m_mlLut;

  std::unique_ptr<const Acts::Logger> m_logger =
      Acts::getDefaultLogger("Finder", Acts::Logging::Level::INFO);

  const Acts::Logger& logger() const { return *m_logger; }

  /// Create graph nodes from space points.
  /// @param spacePoints Space point container
  /// @param maxLayers Maximum number of layers
  /// @return Vector of node vectors organized by layer
  std::vector<std::vector<GbtsNode>> createNodes(
      const SpacePointContainer2& spacePoints, std::uint32_t maxLayers) const;

  /// Parse machine learning lookup table from file.
  /// @param lutInputFile Path to the lookup table input file
  /// @return Parsed machine learning lookup table
  GbtsMlLookupTable parseGbtsMlLookupTable(const std::string& lutInputFile);

  /// Build doublet graph from nodes.
  /// @param roi Region of interest descriptor
  /// @param nodeStorage Data storage containing nodes
  /// @param edgeStorage Storage for generated edges
  /// @param options Event based options such as magnetic field strength
  /// @return Pair of edge count and maximum level
  std::pair<std::int32_t, std::int32_t> buildTheGraph(
      const GbtsRoiDescriptor& roi, GbtsNodeStorage& nodeStorage,
      std::vector<GbtsEdge>& edgeStorage, const Options& options) const;

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
  /// @param vOutputSeeds Output vector for seed candidates
  /// @param filter Tracking filter to be applied
  void extractSeedsFromTheGraph(std::uint32_t maxLevel, std::uint32_t nEdges,
                                std::int32_t nHits,
                                std::vector<GbtsEdge>& edgeStorage,
                                std::vector<OutputSeedProperties>& vOutputSeeds,
                                const GbtsTrackingFilter& filter) const;

  /// Check to see if z0 of segment is within the expected z range of the
  /// beamspot
  /// @param z0BitMask Sets allowed bins of allowed z value
  /// @param z0 Estimated z0 of segments z value at beamspot
  /// @param minZ0 Minimum value of beam spot z coordinate
  /// @param z0HistoCoeff Scalfactor that converts z coodindate into bin index
  /// @return Whether segment is within beamspot range
  bool checkZ0BitMask(std::uint16_t z0BitMask, float z0, float minZ0,
                      float z0HistoCoeff) const;

  float estimateCurvature(const std::array<const GbtsNode*, 3>& nodes) const;

  bool validateTriplet(const std::array<const GbtsNode*, 3> candidateTriplet,
                       float tripletMinPt, float tauRatio, float tauRatioCut,
                       const Options& options) const;
};

}  // namespace Acts::Experimental
