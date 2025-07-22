// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding/SeedFinderUtils.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <memory>
#include <vector>

namespace Acts::Experimental {

class DoubletSeedFinder {
 public:
  struct Config {
    /// Direction of the doublet candidate space points. Either forward, also
    /// called top doublet, or backward, also called bottom doublet.
    Direction candidateDirection = Direction::Forward();

    /// Minimum radial distance between two doublet components
    float deltaRMin = 5 * Acts::UnitConstants::mm;
    /// Maximum radial distance between two doublet components
    float deltaRMax = 270 * Acts::UnitConstants::mm;

    /// Minimal z distance between two doublet components
    float deltaZMin = -std::numeric_limits<float>::infinity();
    /// Maximum z distance between two doublet components
    float deltaZMax = std::numeric_limits<float>::infinity();

    /// Maximum value of impact parameter estimation of the seed candidates
    float impactMax = 20 * UnitConstants::mm;

    /// Enable cut on the compatibility between interaction point and doublet,
    /// this is an useful approximation to speed up the seeding
    bool interactionPointCut = false;

    /// Limiting location of collision region in z-axis used to check if doublet
    /// origin is within reasonable bounds
    float collisionRegionMin = -150 * UnitConstants::mm;
    float collisionRegionMax = +150 * UnitConstants::mm;

    /// Maximum allowed cotTheta between two space-points in doublet, used to
    /// check if forward angle is within bounds
    float cotThetaMax = 10.01788;  // equivalent to eta = 3 (pseudorapidity)

    /// Minimum transverse momentum (pT) used to check the r-z slope
    /// compatibility of triplets with maximum multiple scattering effect
    /// (produced by the minimum allowed pT particle) + a certain uncertainty
    /// term. Check the documentation for more information
    /// https://acts.readthedocs.io/en/latest/core/reconstruction/pattern_recognition/seeding.html
    float minPt = 400 * UnitConstants::MeV;
    /// Parameter which can loosen the tolerance of the track seed to form a
    /// helix. This is useful for e.g. misaligned seeding.
    float helixCutTolerance = 1;

    /// Delegate to apply experiment specific cuts during doublet finding
    Delegate<bool(float /*bottomRadius*/, float /*cotTheta*/)> experimentCuts;

    /// Whether the input space points are sorted by radius
    bool spacePointsSortedByRadius = false;
  };

  struct DerivedConfig : public Config {
    DerivedConfig(const Config& config, float bFieldInZ);

    float bFieldInZ = 0;
    float minHelixDiameter2 = std::numeric_limits<float>::quiet_NaN();
  };

  struct DoubletsForMiddleSp {
    std::vector<SpacePointIndex2> spacePoints;
    /// contains parameters required to calculate a circle with linear equation
    std::vector<LinCircle> linCircles;
    std::vector<float> cotTheta;

    [[nodiscard]] bool empty() const { return spacePoints.empty(); }

    [[nodiscard]] std::size_t size() const { return spacePoints.size(); }

    void clear() {
      spacePoints.clear();
      linCircles.clear();
      cotTheta.clear();
    }

    void emplace_back(SpacePointIndex2 sp, const LinCircle& linCircle) {
      spacePoints.emplace_back(sp);
      linCircles.emplace_back(linCircle);
      cotTheta.emplace_back(linCircle.cotTheta);
    }
  };

  struct MiddleSpInfo {
    /// minus one over radius of middle SP
    float uIP{};
    /// square of uIP
    float uIP2{};
    /// ratio between middle SP x position and radius
    float cosPhiM{};
    /// ratio between middle SP y position and radius
    float sinPhiM{};
  };

  static MiddleSpInfo computeMiddleSpInfo(const ConstSpacePointProxy2& spM);

  explicit DoubletSeedFinder(const DerivedConfig& config);

  const DerivedConfig& config() const { return m_impl->config(); }

  /// Creates compatible dublets by applying a series of cuts that can be
  /// tested with only two SPs.
  ///
  /// @param middleSp Space point candidate to be used as middle SP in a seed
  /// @param middleSpInfo Information about the middle space point
  /// @param candidateSps Subset of space points to be used as candidates for
  ///   middle SP in a seed
  /// @param compatibleDoublets Output container for compatible doublets
  void createDoublets(const ConstSpacePointProxy2& middleSp,
                      const MiddleSpInfo& middleSpInfo,
                      SpacePointContainer2::ConstSubset& candidateSps,
                      DoubletsForMiddleSp& compatibleDoublets) const {
    m_impl->createDoublets(middleSp, middleSpInfo, candidateSps,
                           compatibleDoublets);
  }

  /// Creates compatible dublets by applying a series of cuts that can be
  /// tested with only two SPs.
  ///
  /// @param middleSp Space point candidate to be used as middle SP in a seed
  /// @param middleSpInfo Information about the middle space point
  /// @param candidateSps Range of space points to be used as candidates for
  ///   middle SP in a seed
  /// @param compatibleDoublets Output container for compatible doublets
  void createDoublets(const ConstSpacePointProxy2& middleSp,
                      const MiddleSpInfo& middleSpInfo,
                      SpacePointContainer2::ConstRange& candidateSps,
                      DoubletsForMiddleSp& compatibleDoublets) const {
    m_impl->createDoublets(middleSp, middleSpInfo, candidateSps,
                           compatibleDoublets);
  }

 private:
  class ImplBase {
   public:
    explicit ImplBase(const DerivedConfig& config) : m_cfg(config) {}
    virtual ~ImplBase() = default;

    const DerivedConfig& config() const { return m_cfg; }

    virtual void createDoublets(
        const ConstSpacePointProxy2& middleSp, const MiddleSpInfo& middleSpInfo,
        SpacePointContainer2::ConstSubset& candidateSps,
        DoubletsForMiddleSp& compatibleDoublets) const = 0;

    virtual void createDoublets(
        const ConstSpacePointProxy2& middleSp, const MiddleSpInfo& middleSpInfo,
        SpacePointContainer2::ConstRange& candidateSpRange,
        DoubletsForMiddleSp& compatibleDoublets) const = 0;

   protected:
    DerivedConfig m_cfg;
  };
  template <bool isBottomCandidate, bool interactionPointCut, bool sortedInR>
  class Impl;

  static std::shared_ptr<ImplBase> makeImpl(const DerivedConfig& config);

  std::shared_ptr<ImplBase> m_impl;
};

}  // namespace Acts::Experimental
