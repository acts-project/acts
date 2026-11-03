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
#include "Acts/Seeding/GbtsNodeStorage.hpp"
#include "Acts/Seeding/GbtsRoiDescriptor.hpp"
#include "Acts/Seeding/GbtsSeedFit.hpp"
#include "Acts/Seeding/GbtsTauLookupTable.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

namespace Acts {

namespace detail {
class GbtsSeederImpl;
}
class GbtsGeometry;

/// Seed finder implementing the GBTS seeding workflow.
/// Concurrent calls may share a seeder, using separate storage and outputs.
class GraphBasedTrackSeeder {
 public:
  /// Doublet and triplet acceptance configuration.
  struct GraphConfig {
    /// Match seeds before creating them.
    bool matchBeforeCreate = false;

    /// Optional validation for barrel triplets.
    bool validateTriplets = true;

    /// Tau ratio scale of the `validateTriplets` high-pT check. The linking cut
    /// comes from the layer pairs, see `GbtsLayerConnection`.
    float tauRatioCut = 0.007f;

    /// Tau ratio precut threshold.
    float tauRatioPrecut = 0.009f;

    /// Minimum transverse momentum.
    float minPt = 1.0f * Acts::UnitConstants::GeV;

    /// Fraction of `minPt` a triplet may fall to, allowing for three-point
    /// pT resolution.
    float tripletPtFraction = 0.8f;

    // Graph-building options

    /// Use eta binning from geometry structure.
    bool useEtaBinning = true;

    /// Apply RZ cuts on doublets.
    bool doubletFilterRZ = true;

    /// Maximum number of GBTS edges/doublets per space point.
    float maxEdgesPerSP = 30.0f;

    /// Minimum delta radius between layers.
    float minDeltaRadius = 2.0f * Acts::UnitConstants::mm;

    /// Largest |cot(theta)| accepted for a doublet. The default corresponds to
    /// |eta| of about 4.3, beyond the acceptance of any current tracker.
    float maxAbsTau = 36.0f;

    /// Maximum d0 impact parameter when validating an edge-connection triplet.
    float d0Max = 3.0f * Acts::UnitConstants::mm;

    /// Maximum difference in allowed tangent between candidate edge
    /// connections.
    float cutDPhiMax = 0.012f;

    /// Maximum allowed curvature tolerance for candidate edge connections.
    float cutDCurvMax = 0.001f;

    /// Reject a candidate edge connection whose tau ratio, phi and curvature
    /// differences, each divided by its own cut, add up to more than
    /// `sumRatioCutMax`: a pair that is marginal in all three at once.
    bool useSumRatioCut = false;

    /// Threshold of the summed ratios, between 1 (only one difference
    /// can use its full cut) and 3 (no effect beyond the individual cuts).
    float sumRatioCutMax = 1.3f;

    /// Minimum z0 value. In pixel mode the value is picked from the RoI.
    float minZ0 = -600.0f;

    /// Maximum z0 value. In pixel mode the value is picked from the RoI.
    float maxZ0 = 600.0f;

    /// pT at which the default cut coefficients were tuned; they scale by
    /// `tuningPt / minPt`.
    float tuningPt = 0.9f * Acts::UnitConstants::GeV;

    /// Maximum |curvature| above `curvatureSplitAbsTau`, before that scaling.
    float maxCurvatureHighEta = 4.75e-4f / Acts::UnitConstants::mm;

    /// Maximum |curvature| below `curvatureSplitAbsTau`, before that scaling.
    float maxCurvatureLowEta = 3.75e-4f / Acts::UnitConstants::mm;

    /// |cot(theta)| separating the two curvature cuts, corresponding to |eta|
    /// of about 2.1.
    float curvatureSplitAbsTau = 4.0f;

    /// Radial separation splitting the two phi-window slopes below.
    float phiWindowSplitDeltaRadius = 60.0f * Acts::UnitConstants::mm;

    /// Phi window below `phiWindowSplitDeltaRadius`, as an offset plus a slope
    /// times the radial separation.
    float phiWindowNearOffset = 0.002f;

    /// Slope of the near phi window per unit radial separation. Scaled by
    /// `tuningPt / minPt`.
    float phiWindowNearSlope = 4.33e-4f / Acts::UnitConstants::mm;

    /// Phi window above `phiWindowSplitDeltaRadius`, in the same form.
    float phiWindowFarOffset = 0.015f;

    /// Slope of the far phi window per unit radial separation. Scaled by
    /// `tuningPt / minPt`.
    float phiWindowFarSlope = 2.2e-4f / Acts::UnitConstants::mm;

    /// Incoming edge count below which a node is accepted without a tau match.
    std::uint32_t matchBeforeCreateMaxEdges = 2;

    /// Highest pixel barrel layer, counted inside out, whose nodes are cut
    /// against the z0 histogram of their outer neighbourhood and whose isolated
    /// nodes are skipped. A negative value disables the cut.
    std::int32_t z0HistogramMaxBarrelOrder = 0;

    /// Highest pixel barrel layer, counted inside out, to which
    /// `matchBeforeCreate` applies when it is enabled. A negative value
    /// disables it.
    std::int32_t matchBeforeCreateMaxBarrelOrder = 1;

    /// Half-width of the z0 window against which a node is matched in the
    /// histogram.
    float z0Resolution = 2.5f * Acts::UnitConstants::mm;

    /// Maximum radius of the pixel detector.
    float maxOuterRadius = 550.0f;

    /// Resolve a doublet's strip endpoints against its own direction before
    /// cutting on them. Nothing is written back; the correction belongs to the
    /// pair.
    bool calibrateStrips = true;

    /// How far along a strip a crossing may land and still be recovered, as a
    /// multiple of the strip half-length, so 1 is the strip itself. This is the
    /// same quantity as `TripletSeedFinder::Config::toleranceParam`.
    float maxStripLengthFraction = 1.1f;
  };

  /// Path tracking filter configuration.
  struct FilterConfig {
    /// Multiple scattering sigma (for 900 MeV track at eta=0).
    float sigmaMS = 0.016;
    /// Radiation length fraction per layer (2.5% per layer).
    float radLen = 0.025;

    /// Measurement uncertainty in x direction.
    float sigmaX = 0.08 * Acts::UnitConstants::mm;
    /// Measurement uncertainty in y direction.
    float sigmaY = 0.25 * Acts::UnitConstants::mm;

    /// Across the strips only. The along-strip resolution is not configured:
    /// it is derived per node from the node's own strip. See `update`.
    float sigmaXStrip = 0.03 * Acts::UnitConstants::mm;

    /// Measurement weight in x direction.
    float weightX = 0.5;
    /// Measurement weight in y direction.
    float weightY = 0.5;

    /// Maximum delta chi2 in x direction.
    float maxDChi2X = 5.0;
    /// Maximum delta chi2 in y direction.
    float maxDChi2Y = 6.0;

    /// Chi2 penalty for adding a hit.
    float addHit = 14.0;

    /// Maximum track curvature.
    float maxCurvature = 1e-3f / Acts::UnitConstants::mm;
    /// Maximum longitudinal impact parameter.
    float maxZ0 = 170.0 * Acts::UnitConstants::mm;

    /// Initial variance of the transverse state (y, dy/dx, d2y/dx2).
    std::array<float, 3> initialVarianceX = {0.25f, 0.001f, 0.001f};
    /// Initial variance of the longitudinal state (z, dz/dr).
    std::array<float, 2> initialVarianceY = {1.5f, 0.001f};
  };

  /// Configuration struct for the GBTS seeding algorithm.
  struct Config {
    /// Graph construction cuts.
    GraphConfig graph;
    /// Path tracking filter cuts.
    FilterConfig filter;
    /// Accepted tau range per cluster width bin, needed by the cluster width
    /// cuts and ignored without them.
    GbtsTauLookupTable tauLookupTable;
    /// Enable the cluster width cuts: wide endcap rejection and tau narrowing.
    bool useClusterWidthCuts = false;
    /// Maximum number of phi slices.
    float nMaxPhiSlice = 53;  // used to calculate phi slices

    // Seed extraction options

    /// Longest chain of edges taken into account. Edges in longer chains take
    /// no part in the seed extraction.
    std::uint32_t maxChainLength = 15;
    /// Chain length a seed candidate must reach: a triplet plus one
    /// confirmation.
    std::uint8_t minSeedLevel = 3;
    /// optionally add 3 sp seeds within a certain eta range
    ///
    /// @note Matters most where there are few layers.
    bool addTriplets = false;
    /// the maximum allowed eta value in which
    /// three spacepoint seeds are passed through
    float maxAbsEtaAddTriplets = 1.5;
    /// Smallest seed size that is split into drop-out candidates.
    std::uint32_t minSplitSeedSize = 4;
    /// Largest seed size that is split.
    std::uint32_t maxSplitSeedSize = 5;
    /// Threshold for hit sharing between seeds.
    float hitShareThreshold = 0.49;
    /// Max seed eta value considered for splitting.
    float maxSeedSplitEta = 0.6;
    /// Max allowed curvature for seed self consistency check.
    float maxInvRadDiff = 0.7e-2 / UnitConstants::m;
    /// Reject a complete chain whose fitted |curvature| times this factor
    /// exceeds `FilterConfig::maxCurvature`: a tighter minimum
    /// pT once the whole chain has been fitted. Values below 1 does nothing.
    float postFitCurvatureFactor = 0.0f;
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
    /// Buckets used to sort an eta bin by phi, at most 31.
    std::uint32_t phiSortBuckets = 31;
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
  GraphBasedTrackSeeder(const Config& config,
                        std::shared_ptr<const GbtsGeometry> geometry,
                        std::unique_ptr<const Acts::Logger> logger =
                            Acts::getDefaultLogger("Finder",
                                                   Acts::Logging::Level::INFO));

  /// Destroy the seeder and its implementation.
  ~GraphBasedTrackSeeder();
  /// Move the seeder without moving its implementation.
  GraphBasedTrackSeeder(GraphBasedTrackSeeder&&) noexcept;
  /// Move-assign the seeder.
  GraphBasedTrackSeeder& operator=(GraphBasedTrackSeeder&&) noexcept;
  GraphBasedTrackSeeder(const GraphBasedTrackSeeder&) = delete;
  GraphBasedTrackSeeder& operator=(const GraphBasedTrackSeeder&) = delete;

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
  /// @param options Event based options such as magnetic field strength
  /// @param outputSeeds Container with generated seeds
  void createSeeds(const SpacePointContainer& spacePoints,
                   const GbtsRoiDescriptor& roi, const Options& options,
                   SeedContainer& outputSeeds) const;

  /// Create seeds from a finalized node storage in a region of interest.
  /// @param nodeStorage Finalized storage from the same configuration and geometry
  /// @param roi Region of interest descriptor
  /// @param options Event based options such as magnetic field strength
  /// @param outputSeeds Container with generated seeds
  /// @param outputFits Fit of the tracking filter to each generated seed
  void createSeeds(
      GbtsNodeStorage& nodeStorage, const GbtsRoiDescriptor& roi,
      const Options& options, SeedContainer& outputSeeds,
      std::vector<std::optional<GbtsSeedFit>>* outputFits = nullptr) const;

 private:
  std::unique_ptr<detail::GbtsSeederImpl> m_impl;
};

}  // namespace Acts
