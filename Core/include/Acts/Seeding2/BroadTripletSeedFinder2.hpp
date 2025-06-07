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
#include "Acts/Seeding2/BroadTripletSeedFilter2.hpp"
#include "Acts/Seeding2/DoubletSeedFinder2.hpp"
#include "Acts/Seeding2/SpacePointContainerPointers2.hpp"
#include "Acts/Seeding2/detail/CandidatesForMiddleSp2.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstdint>
#include <memory>
#include <vector>

namespace Acts {

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
class BroadTripletSeedFinder2 {
 public:
  struct Config {};

  struct Options {
    float bFieldInZ = 2 * UnitConstants::T;

    /// Delegates for accessors to detailed information on double measurement
    /// that produced the space point. This is mainly referring to space points
    /// produced when combining measurement from strips on back-to-back modules.
    /// Enables setting of the following delegates.
    bool useDetailedDoubleMeasurementInfo = false;

    BroadTripletSeedFilter2::Options filter;
  };

  struct DerivedTripletCuts;

  struct TripletCuts {
    /// Minimum transverse momentum (pT) used to check the r-z slope
    /// compatibility of triplets with maximum multiple scattering effect
    /// (produced by the minimum allowed pT particle) + a certain uncertainty
    /// term. Check the documentation for more information
    /// https://acts.readthedocs.io/en/latest/core/reconstruction/pattern_recognition/seeding.html
    float minPt = 400. * UnitConstants::MeV;
    /// Number of sigmas of scattering angle to be considered in the minimum pT
    /// scattering term
    float sigmaScattering = 5;
    /// Term that accounts for the thickness of scattering medium in radiation
    /// lengths in the Lynch & Dahl correction to the Highland equation default
    /// is 5%
    /// TODO: necessary to make amount of material dependent on detector region?
    float radLengthPerSeed = 0.05;
    /// Maximum transverse momentum for scattering calculation
    float maxPtScattering = 10 * UnitConstants::GeV;
    /// Maximum value of impact parameter estimation of the seed candidates
    float impactMax = 20. * UnitConstants::mm;
    /// Parameter which can loosen the tolerance of the track seed to form a
    /// helix. This is useful for e.g. misaligned seeding.
    float helixCutTolerance = 1.;

    /// Tolerance parameter used to check the compatibility of space-point
    /// coordinates in xyz. This is only used in a detector specific check for
    /// strip modules
    float toleranceParam = 1.1 * UnitConstants::mm;

    DerivedTripletCuts derive(float bFieldInZ) const;
  };

  struct DerivedTripletCuts : public TripletCuts {
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
    BroadTripletSeedFilter2::Cache filter;

    DoubletSeedFinder2::DoubletsForMiddleSp bottomDoublets;
    DoubletSeedFinder2::DoubletsForMiddleSp topDoublets;

    TripletCache tripletCache;
    TripletTopCandidates tripletTopCandidates;

    CandidatesForMiddleSp2 candidatesCollector;
    std::vector<TripletCandidate2> sortedCandidates;
  };

  struct State {
    BroadTripletSeedFilter2::State filter;
  };

  explicit BroadTripletSeedFinder2(
      const Config& config,
      std::unique_ptr<const Logger> logger =
          getDefaultLogger("BroadTripletSeedFinder2", Logging::Level::INFO));

  const Config& config() const { return m_cfg; }

  /// Create all possible seeds from bottom, middle, and top space points.
  ///
  /// @param state State of the seed finder
  /// @param cache Cache object to store intermediate results
  /// @param containerPointers Space point container and its extra columns
  /// @param bottomSps Group of space points to be used as innermost SP in a seed
  /// @param middleSp Space point candidate to be used as middle SP in a seed
  /// @param topSps Group of space points to be used as outermost SP in a seed
  /// @param outputSeeds Output container for the seeds
  void createSeedsFromGroup(
      const Options& options, State& state, Cache& cache,
      const DoubletSeedFinder2::DerivedCuts& bottomCuts,
      const DoubletSeedFinder2::DerivedCuts& topCuts,
      const DerivedTripletCuts& tripletCuts,
      const BroadTripletSeedFilter2& filter,
      const SpacePointContainerPointers2& containerPointers,
      std::span<const SpacePointIndex2> bottomSps, SpacePointIndex2 middleSp,
      std::span<const SpacePointIndex2> topSps,
      SeedContainer2& outputSeeds) const;

 private:
  /// Create triplets from the bottom, middle, and top space points.
  ///
  /// @param cache Cache object to store intermediate results
  /// @param cuts Triplet cuts that define the compatibility of space points
  /// @param filter Triplet seed filter that defines the filtering criteria
  /// @param filterOptions Options for the triplet seed filter
  /// @param filterState State object that holds the state of the filter
  /// @param filterCache Cache object that holds memory used in SeedFilter
  /// @param containerPointers Space point container and its extra columns
  /// @param spM Space point candidate to be used as middle SP in a seed
  /// @param bottomDoublets Bottom doublets to be used for triplet creation
  /// @param topDoublets Top doublets to be used for triplet creation
  /// @param tripletTopCandidates Cache for triplet top candidates
  /// @param candidatesCollector Collector for candidates for middle space points
  static void createTriplets(
      TripletCache& cache, const DerivedTripletCuts& cuts,
      const BroadTripletSeedFilter2& filter,
      const BroadTripletSeedFilter2::Options& filterOptions,
      BroadTripletSeedFilter2::State& filterState,
      BroadTripletSeedFilter2::Cache& filterCache,
      const SpacePointContainerPointers2& containerPointers,
      const ConstSpacePointProxy2& spM,
      const DoubletSeedFinder2::DoubletsForMiddleSp& bottomDoublets,
      const DoubletSeedFinder2::DoubletsForMiddleSp& topDoublets,
      TripletTopCandidates& tripletTopCandidates,
      CandidatesForMiddleSp2& candidatesCollector);

  /// Create triplets from the bottom, middle, and top space points.
  ///
  /// @param cuts Triplet cuts that define the compatibility of space points
  /// @param filter Triplet seed filter that defines the filtering criteria
  /// @param filterOptions Options for the triplet seed filter
  /// @param filterState State object that holds the state of the filter
  /// @param filterCache Cache object that holds memory used in SeedFilter
  /// @param containerPointers Space point container and its extra columns
  /// @param spM Space point candidate to be used as middle SP in a seed
  /// @param bottomDoublets Bottom doublets to be used for triplet creation
  /// @param topDoublets Top doublets to be used for triplet creation
  /// @param tripletTopCandidates Cache for triplet top candidates
  /// @param candidatesCollector Collector for candidates for middle space points
  static void createTripletsDetailed(
      const DerivedTripletCuts& cuts, const BroadTripletSeedFilter2& filter,
      const BroadTripletSeedFilter2::Options& filterOptions,
      BroadTripletSeedFilter2::State& filterState,
      BroadTripletSeedFilter2::Cache& filterCache,
      const SpacePointContainerPointers2& containerPointers,
      const ConstSpacePointProxy2& spM,
      const DoubletSeedFinder2::DoubletsForMiddleSp& bottomDoublets,
      const DoubletSeedFinder2::DoubletsForMiddleSp& topDoublets,
      TripletTopCandidates& tripletTopCandidates,
      CandidatesForMiddleSp2& candidatesCollector);

  /// Check the compatibility of strip space point coordinates in xyz assuming
  /// the Bottom-Middle direction with the strip measurement details
  static bool stripCoordinateCheck(
      float tolerance, const ConstSpacePointProxy2& sp,
      const SpacePointContainerPointers2& containerPointers,
      const Vector3& spacePointPosition, Vector3& outputCoordinates);

 private:
  const Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts
