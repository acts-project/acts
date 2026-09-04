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
#include <string>
#include <utility>
#include <vector>

namespace Acts::Experimental {

/// Seed finder implementing the GBTS seeding workflow.
class GraphBasedTrackSeeder {
 public:
  /// Configuration struct for the GBTS seeding algorithm.
  struct Config {
    /// Enable beam spot correction.
    bool beamSpotCorrection = false;

    /// Look-up table input file path.
    std::string lutInputFile;

    /// Take the strip-to-strip layer connections from the connector file
    /// instead of the pixel-to-pixel ones. Read where the file is loaded, not
    /// by the seeder itself.
    bool useStripConnections = false;
    /// Enable the cluster width cuts: wide endcap rejection and tau narrowing.
    bool useClusterWidthCuts = false;
    /// Match seeds before creating them.
    bool matchBeforeCreate = false;
    /// optional validation for barrel triplets
    bool validateTriplets = true;
    /// widens allowed variation in tau ratio
    /// if layer is missed in edge connecting
    bool useAdaptiveCuts = true;
    /// optionally add 3 sp seeds within a certain eta range
    ///
    /// @note Worth little until `maxAbsEtaAddTriplets` is opened past
    ///       `edgeMaskMinEta`; matters most where there are few layers.
    bool addTriplets = false;
    /// Tau ratio cut threshold.
    float tauRatioCut = 0.007;
    /// Tau ratio precut threshold.
    float tauRatioPrecut = 0.009f;
    /// correction applied to tau acceptance
    /// if a layer is missed during edge connecting
    float tauRatioCorr = 0.006;
    /// The same for a triplet any of whose three nodes a strip module made,
    /// whose two doublets resolved the shared node's along-strip coordinate
    /// separately. Reaches nothing without a strip in the triplet.
    float tauRatioCorrStrip = 0.03f;
    /// the maximum allowed eta value in which
    /// three spacepoint seeds are passed through
    float maxAbsEtaAddTriplets = 1.5;
    /// Maximum number of phi slices.
    float nMaxPhiSlice = 53;  // used to calculate phi slices
    /// Minimum transverse momentum.
    float minPt = 1.0f * UnitConstants::GeV;
    /// Fraction of `minPt` a triplet may fall to, allowing for three-point
    /// pT resolution.
    float tripletPtFraction = 0.8f;

    // graph building options
    /// Use eta binning from geometry structure.
    bool useEtaBinning = true;
    /// Apply RZ cuts on doublets.
    bool doubletFilterRZ = true;
    /// Maximum number of Gbts edges/doublets.
    std::uint32_t nMaxEdges = 2000000;
    /// Minimum delta radius between layers.
    float minDeltaRadius = 2.0 * Acts::UnitConstants::mm;
    /// Largest |cot(theta)| accepted for a doublet. The default corresponds to
    /// |eta| of about 4.3, beyond the acceptance of any current tracker.
    float maxAbsTau = 36.0f;
    /// Maximum d0 impact parameter when validating edge connection triplet
    float d0Max = 3.0 * UnitConstants::mm;
    /// Maximum difference in allowed tangent between candidate edge connection
    float cutDPhiMax = 0.012f;
    /// Maximum allowed curvature tolerance for candidate edge connections
    float cutDCurvMax = 0.001f;
    /// Minimum z0 value. In pixel mode the value is picked from the RoI.
    float minZ0 = -600;
    /// Maximum z0 value. In pixel mode the value is picked from the RoI.
    float maxZ0 = 600;

    /// pT the default cut coefficients were tuned at; they scale by
    /// `tuningPt / minPt`.
    float tuningPt = 0.9f * UnitConstants::GeV;
    /// Maximum |curvature| above `curvatureSplitAbsTau`, before that scaling.
    float maxCurvatureHighEta = 4.75e-4f / UnitConstants::mm;
    /// Maximum |curvature| below `curvatureSplitAbsTau`, before that scaling.
    float maxCurvatureLowEta = 3.75e-4f / UnitConstants::mm;
    /// |cot(theta)| separating the two curvature cuts, |eta| of about 2.1.
    float curvatureSplitAbsTau = 4.0f;
    /// Radial separation splitting the two phi window slopes below.
    float phiWindowSplitDeltaRadius = 60.0f * UnitConstants::mm;
    /// Phi window below `phiWindowSplitDeltaRadius`, as offset plus slope times
    /// the radial separation.
    float phiWindowNearOffset = 0.002f;
    /// Slope of the near phi window, per unit radial separation. Scaled by
    /// `tuningPt / minPt`.
    float phiWindowNearSlope = 4.33e-4f / UnitConstants::mm;
    /// Phi window above `phiWindowSplitDeltaRadius`, in the same form.
    float phiWindowFarOffset = 0.015f;
    /// Slope of the far phi window, per unit radial separation. Scaled by
    /// `tuningPt / minPt`.
    float phiWindowFarSlope = 2.2e-4f / UnitConstants::mm;
    /// Incoming edge count below which a node is accepted without a tau match.
    std::uint32_t matchBeforeCreateMaxEdges = 2;
    /// Layers whose nodes are cut against the z0 histogram of their outer
    /// neighbourhood, and whose isolated nodes are skipped.
    std::vector<std::uint32_t> z0HistogramLayerIds{80000};
    /// Layers `matchBeforeCreate` applies to, when it is enabled.
    std::vector<std::uint32_t> matchBeforeCreateLayerIds{80000, 81000};
    /// Half-width of the z0 window a node is matched against in the histogram.
    float z0Resolution = 2.5f * UnitConstants::mm;
    /// Maximum radius of pixel detector
    float maxOuterRadius = 550.0f;

    /// Resolve a doublet's strip endpoints against its own direction before
    /// cutting on them. Nothing is written back; the correction is the pair's.
    bool calibrateStrips = true;
    /// How far along a strip a crossing may land and still be recovered, as a
    /// multiple of the strip half-length, so 1 is the strip itself. The same
    /// quantity as `TripletSeedFinder::Config::toleranceParam`.
    float maxStripLengthFraction = 1.1f;

    // Seed extraction options
    /// Maximum number of connected-component iterations.
    std::uint32_t ccaMaxIterations = 15;
    /// Chain length a seed candidate must reach: a triplet plus one
    /// confirmation.
    std::uint32_t minSeedLevel = 3;
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
    float maxEndcapClusterWidth = 0.35 * Acts::UnitConstants::mm;
    /// Half-length in local y of a pixel module, against which the distance of
    /// a cluster to the module edge is measured.
    float moduleHalfLengthY = 10.0 * Acts::UnitConstants::mm;
    /// Distance to the module edge below which a cluster may be shortened,
    /// which switches to the tau lookup table's near-edge bounds.
    float moduleEdgeTolerance = 0.3 * Acts::UnitConstants::mm;
    /// Cluster width covered by one bin of the tau lookup table.
    float tauLutBinWidth = 0.05 * Acts::UnitConstants::mm;
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
    /// @param bFieldInZ_ the magnetic field in z
    explicit Options(float bFieldInZ_);

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
  /// container must carry the `layerId`, `clusterWidth` and `localPositionY`
  /// columns.
  /// @param spacePoints Space point container
  /// @param roi Region of interest descriptor
  /// @param filter Tracking filter to be applied
  /// @param options Event based options such as magnetic field strength
  /// @param outputSeeds Container with generated seeds
  void createSeeds(const SpacePointContainer& spacePoints,
                   const GbtsRoiDescriptor& roi,
                   const GbtsTrackingFilter& filter, const Options& options,
                   SeedContainer& outputSeeds) const;

  /// Create seeds from a finalized node storage in a region of interest.
  /// @param nodeStorage Finalized graph node storage
  /// @param roi Region of interest descriptor
  /// @param filter Tracking filter to be applied
  /// @param options Event based options such as magnetic field strength
  /// @param outputSeeds Container with generated seeds
  void createSeeds(GbtsNodeStorage& nodeStorage, const GbtsRoiDescriptor& roi,
                   const GbtsTrackingFilter& filter, const Options& options,
                   SeedContainer& outputSeeds) const;

 private:
  /// candidate seed metadata produced by the GBTS algorithm.
  struct SeedCandidateProperties {
    /// @param quality Seed quality score
    /// @param clone Clone flag
    /// @param sps Vector of graph node indices
    /// @param splitFlag used to flag if seed needs to be split in two
    SeedCandidateProperties(float quality, std::int32_t clone,
                            std::vector<SpacePointIndex> sps,
                            std::uint32_t splitFlag)
        : seedQuality(quality),
          isClone(clone),
          nodes(std::move(sps)),
          forSeedSplitting(splitFlag) {}

    /// Seed quality score.
    float seedQuality{};
    /// Clone flag.
    std::int32_t isClone{};
    /// Graph node indices.
    std::vector<SpacePointIndex> nodes;
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

  /// Sliding window in phi used to define range used for edge creation.
  ///
  /// Covers one non-empty source eta bin, whose phi-ordered node list it holds
  /// directly so that the innermost loop does not reach through the bin.
  struct SlidingWindow {
    /// phi-ordered nodes of the bin, including the wrap-around duplicates
    const std::pair<float, SpacePointIndex>* phiNodes{};
    /// number of entries in @c phiNodes
    std::uint32_t numPhiNodes{};
    /// sliding window position
    std::uint32_t firstIt{};
    /// window half-width;
    float deltaPhi{};
    /// GBTS layer ID of the bin
    std::uint32_t layerId{};
    /// Type of the bin's layer.
    GbtsLayerType type{};
    /// Technology of the bin's layer.
    GbtsLayerTechnology technology{};
  };
  DerivedConfig m_cfg;

  std::shared_ptr<const GbtsGeometry> m_geometry;

  detail::GbtsTauLookupTable m_tauLut;

  std::unique_ptr<const Acts::Logger> m_logger =
      Acts::getDefaultLogger("Finder", Acts::Logging::Level::INFO);

  const Acts::Logger& logger() const { return *m_logger; }

  /// Parse the tau lookup table from file.
  /// @param lutInputFile Path to the lookup table input file
  /// @return Parsed tau lookup table
  detail::GbtsTauLookupTable parseTauLookupTable(
      const std::string& lutInputFile) const;

  /// Build doublet graph from nodes.
  /// @param roi Region of interest descriptor
  /// @param nodeStorage Data storage containing nodes
  /// @param edgeStorage Storage for generated edges
  /// @param options Event based options such as magnetic field strength
  /// @return Pair of edge count and maximum level
  std::pair<std::int32_t, std::int32_t> buildTheGraph(
      const GbtsRoiDescriptor& roi, GbtsNodeStorage& nodeStorage,
      std::vector<detail::GbtsEdge>& edgeStorage, const Options& options) const;

  /// Run connected component analysis on the graph.
  /// @param nEdges Number of edges in the graph
  /// @param edgeStorage Storage containing graph edges
  /// @return Number of connected components found
  std::int32_t runCCA(std::uint32_t nEdges,
                      std::vector<detail::GbtsEdge>& edgeStorage) const;

  /// Extract seed candidates from the graph.
  /// @param maxLevel Maximum level in the graph
  /// @param nEdges Number of edges
  /// @param nodeStorage Storage containing the graph nodes
  /// @param edgeStorage Storage containing edges
  /// @param vOutputSeeds Output vector for seed candidates
  /// @param filter Tracking filter to be applied
  void extractSeedsFromTheGraph(std::uint32_t maxLevel, std::uint32_t nEdges,
                                const GbtsNodeStorage& nodeStorage,
                                std::vector<detail::GbtsEdge>& edgeStorage,
                                std::vector<OutputSeedProperties>& vOutputSeeds,
                                const GbtsTrackingFilter& filter) const;

  /// Check to see if z0 of segment is within the expected z range of the
  /// beamspot
  /// @param z0BitMask Sets allowed bins of allowed z value
  /// @param z0 Estimated z0 of segments z value at beamspot
  /// @param z0HistoCoeff Scalfactor that converts z coodindate into bin index
  /// @return Whether segment is within beamspot range
  bool checkZ0BitMask(std::uint16_t z0BitMask, float z0,
                      float z0HistoCoeff) const;

  /// Estimate the inverse radius of the circle through three nodes.
  /// @param nodeView View of the node positions and layers
  /// @param nodes The three graph nodes, innermost first
  /// @return The estimated inverse radius
  float estimateCurvature(const detail::GbtsNodeView& nodeView,
                          const std::array<SpacePointIndex, 3>& nodes) const;

  /// Check a triplet against the pT and d0 cuts.
  /// @param nodeView View of the node positions and layers
  /// @param candidateTriplet The three graph nodes
  /// @param tripletMinPt Minimum transverse momentum
  /// @param tauRatio Tau ratio of the triplet
  /// @param tauRatioCut Tau ratio cut threshold
  /// @param options Event based options such as magnetic field strength
  /// @return Whether the triplet is accepted
  bool validateTriplet(const detail::GbtsNodeView& nodeView,
                       const std::array<SpacePointIndex, 3>& candidateTriplet,
                       float tripletMinPt, float tauRatio, float tauRatioCut,
                       const Options& options) const;
};

}  // namespace Acts::Experimental
