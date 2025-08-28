// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding2/DoubletSeedFinder.hpp"
#include "Acts/Utilities/ContainerIterator.hpp"

#include <vector>

namespace Acts::Experimental {

/// Container for triplet candidates found by the triplet seed finder.
///
/// This implementation uses partial AoS/SoA depending on the access pattern in
/// the triplet finding process.
class TripletTopCandidates {
 public:
  using Index = std::uint32_t;

  Index size() const { return static_cast<Index>(m_topSpacePoints.size()); }

  void reserve(Index size) {
    m_topSpacePoints.reserve(size);
    m_curvatures.reserve(size);
    m_impactParameters.reserve(size);
  }

  void clear() {
    m_topSpacePoints.clear();
    m_curvatures.clear();
    m_impactParameters.clear();
  }

  void emplace_back(SpacePointIndex2 spT, float curvature,
                    float impactParameter) {
    m_topSpacePoints.emplace_back(spT);
    m_curvatures.emplace_back(curvature);
    m_impactParameters.emplace_back(impactParameter);
  }

  const std::vector<SpacePointIndex2>& topSpacePoints() const {
    return m_topSpacePoints;
  }
  const std::vector<float>& curvatures() const { return m_curvatures; }
  const std::vector<float>& impactParameters() const {
    return m_impactParameters;
  }

  class Proxy {
   public:
    Proxy(const TripletTopCandidates* container, Index index)
        : m_container(container), m_index(index) {}

    SpacePointIndex2 spacePoint() const {
      return m_container->m_topSpacePoints[m_index];
    }

    float curvature() const { return m_container->m_curvatures[m_index]; }

    float impactParameter() const {
      return m_container->m_impactParameters[m_index];
    }

   private:
    const TripletTopCandidates* m_container{};
    Index m_index{};
  };

  Proxy operator[](Index index) const { return Proxy(this, index); }

  using const_iterator =
      ContainerIterator<TripletTopCandidates, Proxy, Index, true>;

  const_iterator begin() const { return const_iterator(*this, 0); }
  const_iterator end() const { return const_iterator(*this, size()); }

 private:
  std::vector<SpacePointIndex2> m_topSpacePoints;
  std::vector<float> m_curvatures;
  std::vector<float> m_impactParameters;
};

/// Interface and a collection of standard implementations for a triplet seed
/// finder.
///
/// @note The standard implementations rely on virtual function dispatch which
/// did not turn out to affect the performance after measurement.
class TripletSeedFinder {
 public:
  /// Collection of configuration parameters for the triplet seed finder. This
  /// includes triplet cuts, steering switches, and assumptions about the space
  /// points.
  struct Config {
    /// Delegates for accessors to detailed information on double strip
    /// measurement that produced the space point. This is mainly referring to
    /// space points produced when combining measurement from strips on
    /// back-to-back modules. Enables setting of the following delegates.
    bool useStripInfo = false;

    /// Whether the input doublets are sorted by cotTheta
    bool sortedByCotTheta = true;

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

  /// Derived configuration for the triplet seed finder using a magnetic field.
  struct DerivedConfig : public Config {
    DerivedConfig(const Config& config, float bFieldInZ);

    float bFieldInZ = std::numeric_limits<float>::quiet_NaN();
    float highland = std::numeric_limits<float>::quiet_NaN();
    float pTPerHelixRadius = std::numeric_limits<float>::quiet_NaN();
    float minHelixDiameter2 = std::numeric_limits<float>::quiet_NaN();
    float sigmapT2perRadius = std::numeric_limits<float>::quiet_NaN();
    float multipleScattering2 = std::numeric_limits<float>::quiet_NaN();
  };

  /// Creates a new triplet seed finder instance given the configuration.
  static std::unique_ptr<TripletSeedFinder> create(const DerivedConfig& config);

  virtual ~TripletSeedFinder() = default;

  /// Returns the configuration of the triplet seed finder.
  virtual const DerivedConfig& config() const = 0;

  /// Create triplets from the bottom, middle, and top space points.
  ///
  /// @param spacePoints Space point container
  /// @param spM Space point candidate to be used as middle SP in a seed
  /// @param bottomDoublet Bottom doublet to be used for triplet creation
  /// @param topDoublets Top doublets to be used for triplet creation
  /// @param tripletTopCandidates Cache for triplet top candidates
  virtual void createTripletTopCandidates(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp::Proxy& bottomDoublet,
      DoubletsForMiddleSp::Range& topDoublets,
      TripletTopCandidates& tripletTopCandidates) const = 0;

  /// Create triplets from the bottom, middle, and top space points.
  ///
  /// @param spacePoints Space point container
  /// @param spM Space point candidate to be used as middle SP in a seed
  /// @param bottomDoublet Bottom doublet to be used for triplet creation
  /// @param topDoublets Top doublets to be used for triplet creation
  /// @param tripletTopCandidates Cache for triplet top candidates
  virtual void createTripletTopCandidates(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp::Proxy& bottomDoublet,
      DoubletsForMiddleSp::Subset& topDoublets,
      TripletTopCandidates& tripletTopCandidates) const = 0;

  /// Create triplets from the bottom, middle, and top space points.
  ///
  /// @param spacePoints Space point container
  /// @param spM Space point candidate to be used as middle SP in a seed
  /// @param bottomDoublet Bottom doublet to be used for triplet creation
  /// @param topDoublets Top doublets to be used for triplet creation
  /// @param tripletTopCandidates Cache for triplet top candidates
  virtual void createTripletTopCandidates(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp::Proxy& bottomDoublet,
      DoubletsForMiddleSp::Subset2& topDoublets,
      TripletTopCandidates& tripletTopCandidates) const = 0;
};

}  // namespace Acts::Experimental
