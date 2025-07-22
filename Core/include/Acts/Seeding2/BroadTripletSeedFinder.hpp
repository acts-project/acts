// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding2/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding2/DoubletSeedFinder.hpp"
#include "Acts/Seeding2/detail/CandidatesForMiddleSp2.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstdint>
#include <vector>

namespace Acts::Experimental {

/// @brief Triplet seeding algorithm front-end
///
/// This class implements the triplet seeding algorithm, which is typical
/// procedure to find track seeds using space points in a cylindrical detector.
/// It is designed to be fast, flexible, and configurable.
///
/// The algorithm works by first finding compatible doublets of space points,
/// two space points that can be connected by a track coming from the
/// interaction region, and then forming triplets by combinding these
/// doublets at a common middle space point. The triplets are then filtered
/// using a seed filter to produce a set of track seeds.
///
/// Note that this algorithm is designed and tuned for cylindrical detectors and
/// uses R-Z coordinates for the space points.
class BroadTripletSeedFinder {
 public:
  struct Options {
    float bFieldInZ = 2 * UnitConstants::T;

    /// Delegates for accessors to detailed information on double measurement
    /// that produced the space point. This is mainly referring to space points
    /// produced when combining measurement from strips on back-to-back modules.
    /// Enables setting of the following delegates.
    bool useStripMeasurementInfo = false;

    /// Whether the input space points are sorted by radius
    bool spacePointsSortedByRadius = false;
  };

  struct TripletCuts {
    /// Minimum transverse momentum (pT) used to check the r-z slope
    /// compatibility of triplets with maximum multiple scattering effect
    /// (produced by the minimum allowed pT particle) + a certain uncertainty
    /// term. Check the documentation for more information
    /// https://acts.readthedocs.io/en/latest/core/reconstruction/pattern_recognition/seeding.html
    float minPt = 400 * UnitConstants::MeV;
    /// Number of sigmas of scattering angle to be considered in the minimum pT
    /// scattering term
    float sigmaScattering = 5;
    /// Term that accounts for the thickness of scattering medium in radiation
    /// lengths in the Lynch & Dahl correction to the Highland equation default
    /// is 5%
    float radLengthPerSeed = 0.05;
    /// Maximum transverse momentum for scattering calculation
    float maxPtScattering = 10 * UnitConstants::GeV;
    /// Maximum value of impact parameter estimation of the seed candidates
    float impactMax = 20 * UnitConstants::mm;
    /// Parameter which can loosen the tolerance of the track seed to form a
    /// helix. This is useful for e.g. misaligned seeding.
    float helixCutTolerance = 1;

    /// Tolerance parameter used to check the compatibility of space-point
    /// coordinates in xyz. This is only used in a detector specific check for
    /// strip modules
    float toleranceParam = 1.1 * UnitConstants::mm;
  };

  struct DerivedTripletCuts : public TripletCuts {
    DerivedTripletCuts(const TripletCuts& cuts, float bFieldInZ);

    float bFieldInZ = 0;
    float highland = 0;
    float pTPerHelixRadius = std::numeric_limits<float>::quiet_NaN();
    float minHelixDiameter2 = std::numeric_limits<float>::quiet_NaN();
    float sigmapT2perRadius = std::numeric_limits<float>::quiet_NaN();
    float multipleScattering2 = std::numeric_limits<float>::quiet_NaN();
  };

  struct TripletCache {
    std::vector<std::uint32_t> sortedBottoms;
    std::vector<std::uint32_t> sortedTops;
  };

  struct TripletTopCandidates {
    std::vector<SpacePointIndex2> topSpacePoints;
    std::vector<float> curvatures;
    std::vector<float> impactParameters;

    std::size_t size() const { return topSpacePoints.size(); }

    void resize(std::size_t size) {
      topSpacePoints.resize(size);
      curvatures.resize(size);
      impactParameters.resize(size);
    }

    void clear() {
      topSpacePoints.clear();
      curvatures.clear();
      impactParameters.clear();
    }

    void emplace_back(SpacePointIndex2 spT, float curvature,
                      float impactParameter) {
      topSpacePoints.emplace_back(spT);
      curvatures.emplace_back(curvature);
      impactParameters.emplace_back(impactParameter);
    }
  };

  struct Cache {
    DoubletSeedFinder::DoubletsForMiddleSp bottomDoublets;
    DoubletSeedFinder::DoubletsForMiddleSp topDoublets;

    TripletCache tripletCache;
    TripletTopCandidates tripletTopCandidates;

    BroadTripletSeedFilter::Cache filter;

    CandidatesForMiddleSp2 candidatesCollector;
    std::vector<TripletCandidate2> sortedCandidates;
  };

  struct State {
    BroadTripletSeedFilter::State filter;
  };

  explicit BroadTripletSeedFinder(std::unique_ptr<const Logger> logger =
                                      getDefaultLogger("BroadTripletSeedFinder",
                                                       Logging::Level::INFO));

  /// Create all possible seeds from bottom, middle, and top space points.
  ///
  /// @param options Configuration options for the seed finder
  /// @param state State of the seed finder
  /// @param cache Cache object to store intermediate results
  /// @param bottomFinder Finder for bottom doublets
  /// @param topFinder Finder for top doublets
  /// @param tripletCuts Derived cuts for the triplet space points
  /// @param filter Triplet seed filter that defines the filtering criteria
  /// @param spacePoints Space point container
  /// @param bottomSps Subset of space points to be used as innermost SP in a seed
  /// @param middleSp Space point candidate to be used as middle SP in a seed
  /// @param topSps Subset of space points to be used as outermost SP in a seed
  /// @param outputSeeds Output container for the seeds
  void createSeedsFromGroup(const Options& options, State& state, Cache& cache,
                            const DoubletSeedFinder& bottomFinder,
                            const DoubletSeedFinder& topFinder,
                            const DerivedTripletCuts& tripletCuts,
                            const BroadTripletSeedFilter& filter,
                            const SpacePointContainer2& spacePoints,
                            SpacePointContainer2::ConstSubset& bottomSps,
                            const ConstSpacePointProxy2& middleSp,
                            SpacePointContainer2::ConstSubset& topSps,
                            SeedContainer2& outputSeeds) const;

  /// Create all possible seeds from bottom, middle, and top space points.
  ///
  /// @param options Configuration options for the seed finder
  /// @param state State of the seed finder
  /// @param cache Cache object to store intermediate results
  /// @param bottomFinder Finder for bottom doublets
  /// @param topFinder Finder for top doublets
  /// @param tripletCuts Derived cuts for the triplet space points
  /// @param filter Triplet seed filter that defines the filtering criteria
  /// @param spacePoints Space point container
  /// @param bottomSpGroups Groups of space points to be used as innermost SP in a seed
  /// @param middleSpGroup Group of space points to be used as middle SP in a seed
  /// @param topSpGroups Groups of space points to be used as outermost SP in a seed
  /// @param radiusRangeForMiddle Range of radii for the middle space points
  /// @param outputSeeds Output container for the seeds
  void createSeedsFromGroups(
      const Options& options, State& state, Cache& cache,
      const DoubletSeedFinder& bottomFinder, const DoubletSeedFinder& topFinder,
      const DerivedTripletCuts& tripletCuts,
      const BroadTripletSeedFilter& filter,
      const SpacePointContainer2& spacePoints,
      const std::span<SpacePointContainer2::ConstRange>& bottomSpGroups,
      const SpacePointContainer2::ConstRange& middleSpGroup,
      const std::span<SpacePointContainer2::ConstRange>& topSpGroups,
      const std::pair<float, float>& radiusRangeForMiddle,
      SeedContainer2& outputSeeds) const;

 private:
  std::unique_ptr<const Logger> m_logger;

  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts::Experimental
