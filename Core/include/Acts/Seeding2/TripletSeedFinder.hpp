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

class TripletTopCandidates {
 public:
  using Index = std::uint32_t;

  Index size() const { return m_topSpacePoints.size(); }

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

  Proxy operator[](std::size_t index) const { return Proxy(this, index); }

  using const_iterator =
      ContainerIterator<TripletTopCandidates, Proxy, Index, true>;

  const_iterator begin() const { return const_iterator(*this, 0); }
  const_iterator end() const { return const_iterator(*this, size()); }

 private:
  std::vector<SpacePointIndex2> m_topSpacePoints;
  std::vector<float> m_curvatures;
  std::vector<float> m_impactParameters;
};

class TripletSeedFinder {
 public:
  struct Config {
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

    /// Delegates for accessors to detailed information on double measurement
    /// that produced the space point. This is mainly referring to space points
    /// produced when combining measurement from strips on back-to-back modules.
    /// Enables setting of the following delegates.
    bool useStripInfo = false;

    /// Whether the input doublets are sorted by cotTheta
    bool sortedByCotTheta = false;
  };

  struct DerivedConfig : public Config {
    DerivedConfig(const Config& config, float bFieldInZ);

    float bFieldInZ = std::numeric_limits<float>::quiet_NaN();
    float highland = std::numeric_limits<float>::quiet_NaN();
    float pTPerHelixRadius = std::numeric_limits<float>::quiet_NaN();
    float minHelixDiameter2 = std::numeric_limits<float>::quiet_NaN();
    float sigmapT2perRadius = std::numeric_limits<float>::quiet_NaN();
    float multipleScattering2 = std::numeric_limits<float>::quiet_NaN();
  };

  explicit TripletSeedFinder(const DerivedConfig& config);

  const DerivedConfig& config() const { return m_impl->config(); }

  /// Create triplets from the bottom, middle, and top space points.
  ///
  /// @param spacePoints Space point container
  /// @param spM Space point candidate to be used as middle SP in a seed
  /// @param bottomDoublet Bottom doublet to be used for triplet creation
  /// @param topDoublets Top doublets to be used for triplet creation
  /// @param tripletTopCandidates Cache for triplet top candidates
  void createTripletTopCandidates(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp::Proxy& bottomDoublet,
      DoubletsForMiddleSp::Range& topDoublets,
      TripletTopCandidates& tripletTopCandidates) const {
    m_impl->createTripletTopCandidates(spacePoints, spM, bottomDoublet,
                                       topDoublets, tripletTopCandidates);
  }

  /// Create triplets from the bottom, middle, and top space points.
  ///
  /// @param spacePoints Space point container
  /// @param spM Space point candidate to be used as middle SP in a seed
  /// @param bottomDoublet Bottom doublet to be used for triplet creation
  /// @param topDoublets Top doublets to be used for triplet creation
  /// @param tripletTopCandidates Cache for triplet top candidates
  void createTripletTopCandidates(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp::Proxy& bottomDoublet,
      DoubletsForMiddleSp::Subset& topDoublets,
      TripletTopCandidates& tripletTopCandidates) const {
    m_impl->createTripletTopCandidates(spacePoints, spM, bottomDoublet,
                                       topDoublets, tripletTopCandidates);
  }

  /// Create triplets from the bottom, middle, and top space points.
  ///
  /// @param spacePoints Space point container
  /// @param spM Space point candidate to be used as middle SP in a seed
  /// @param bottomDoublet Bottom doublet to be used for triplet creation
  /// @param topDoublets Top doublets to be used for triplet creation
  /// @param tripletTopCandidates Cache for triplet top candidates
  void createTripletTopCandidates(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp::Proxy& bottomDoublet,
      DoubletsForMiddleSp::Subset2& topDoublets,
      TripletTopCandidates& tripletTopCandidates) const {
    m_impl->createTripletTopCandidates(spacePoints, spM, bottomDoublet,
                                       topDoublets, tripletTopCandidates);
  }

 private:
  class ImplBase {
   public:
    explicit ImplBase(const DerivedConfig& config) : m_cfg(config) {}
    virtual ~ImplBase() = default;

    const DerivedConfig& config() const { return m_cfg; }

    virtual void createTripletTopCandidates(
        const SpacePointContainer2& spacePoints,
        const ConstSpacePointProxy2& spM,
        const DoubletsForMiddleSp::Proxy& bottomDoublet,
        DoubletsForMiddleSp::Range& topDoublets,
        TripletTopCandidates& tripletTopCandidates) const = 0;

    virtual void createTripletTopCandidates(
        const SpacePointContainer2& spacePoints,
        const ConstSpacePointProxy2& spM,
        const DoubletsForMiddleSp::Proxy& bottomDoublet,
        DoubletsForMiddleSp::Subset& topDoublets,
        TripletTopCandidates& tripletTopCandidates) const = 0;

    virtual void createTripletTopCandidates(
        const SpacePointContainer2& spacePoints,
        const ConstSpacePointProxy2& spM,
        const DoubletsForMiddleSp::Proxy& bottomDoublet,
        DoubletsForMiddleSp::Subset2& topDoublets,
        TripletTopCandidates& tripletTopCandidates) const = 0;

   protected:
    DerivedConfig m_cfg;
  };
  template <bool useStripInfo, bool sortedInCotTheta>
  class Impl;

  static std::shared_ptr<ImplBase> makeImpl(const DerivedConfig& config);

  std::shared_ptr<ImplBase> m_impl;

  DerivedConfig m_cfg;
};

}  // namespace Acts::Experimental
