// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"

#include <cstdint>
#include <vector>

namespace Acts::Experimental {

class DoubletsForMiddleSp {
 public:
  using Index = std::uint32_t;
  using IndexRange = std::pair<Index, Index>;
  using IndexSubset = std::span<const Index>;

  [[nodiscard]] bool empty() const { return m_spacePoints.empty(); }
  [[nodiscard]] Index size() const {
    return static_cast<Index>(m_spacePoints.size());
  }

  void clear() {
    m_spacePoints.clear();
    m_cotTheta.clear();
    m_er_iDeltaR.clear();
    m_uv.clear();
    m_xy.clear();
  }

  void emplace_back(SpacePointIndex2 sp, float cotTheta, float iDeltaR,
                    float er, float u, float v, float x, float y) {
    m_spacePoints.push_back(sp);
    m_cotTheta.push_back(cotTheta);
    m_er_iDeltaR.push_back({er, iDeltaR});
    m_uv.push_back({u, v});
    m_xy.push_back({x, y});
  }

  const std::vector<SpacePointIndex2>& spacePoints() const {
    return m_spacePoints;
  }
  const std::vector<float>& cotTheta() const { return m_cotTheta; }

  struct IndexAndCotTheta {
    Index index{};
    float cotTheta{};
  };

  using IndexAndCotThetaSubset = std::span<const IndexAndCotTheta>;

  void sortByCotTheta(const IndexRange& range,
                      std::vector<IndexAndCotTheta>& indexAndCotTheta) const {
    indexAndCotTheta.clear();
    indexAndCotTheta.reserve(range.second - range.first);
    for (Index i = range.first; i < range.second; ++i) {
      indexAndCotTheta.emplace_back(i, m_cotTheta[i]);
    }
    std::ranges::sort(indexAndCotTheta, {}, [](const IndexAndCotTheta& item) {
      return item.cotTheta;
    });
  }

  class Proxy {
   public:
    Proxy(const DoubletsForMiddleSp* container, Index index)
        : m_container(container), m_index(index) {}

    const DoubletsForMiddleSp& container() const { return *m_container; }
    Index index() const { return m_index; }

    SpacePointIndex2 spacePointIndex() const {
      return m_container->m_spacePoints[m_index];
    }

    float cotTheta() const { return m_container->m_cotTheta[m_index]; }
    float er() const { return m_container->m_er_iDeltaR[m_index][0]; }
    float iDeltaR() const { return m_container->m_er_iDeltaR[m_index][1]; }
    float u() const { return m_container->m_uv[m_index][0]; }
    float v() const { return m_container->m_uv[m_index][1]; }
    float x() const { return m_container->m_xy[m_index][0]; }
    float y() const { return m_container->m_xy[m_index][1]; }

   private:
    const DoubletsForMiddleSp* m_container{};
    Index m_index{};
  };
  class Proxy2 : public Proxy {
   public:
    Proxy2(const DoubletsForMiddleSp* container,
           IndexAndCotTheta indexAndCotTheta)
        : Proxy(container, indexAndCotTheta.index),
          m_cotTheta(indexAndCotTheta.cotTheta) {}

    float cotTheta() const { return m_cotTheta; }

   private:
    float m_cotTheta{};
  };

  Proxy operator[](Index index) const { return Proxy(this, index); }
  Proxy2 operator[](IndexAndCotTheta indexAndCotTheta) const {
    return Proxy2(this, indexAndCotTheta);
  }

  using const_iterator =
      ContainerIterator<DoubletsForMiddleSp, Proxy, Index, true>;

  const_iterator begin() const { return const_iterator(*this, 0); }
  const_iterator end() const { return const_iterator(*this, size()); }

  class Range
      : public ContainerRange<Range, Range, DoubletsForMiddleSp, Index, true> {
   public:
    using Base = ContainerRange<Range, Range, DoubletsForMiddleSp, Index, true>;

    using Base::Base;
  };

  Range range() const noexcept { return Range(*this, {0, size()}); }
  Range range(const IndexRange& range) const noexcept {
    return Range(*this, range);
  }

  class Subset : public ContainerSubset<Subset, DoubletsForMiddleSp, Proxy,
                                        Index, true> {
   public:
    using Base =
        ContainerSubset<Subset, DoubletsForMiddleSp, Proxy, Index, true>;

    using Base::Base;
  };
  class Subset2 : public ContainerSubset<Subset2, DoubletsForMiddleSp, Proxy2,
                                         IndexAndCotTheta, true> {
   public:
    using Base = ContainerSubset<Subset2, DoubletsForMiddleSp, Proxy2,
                                 IndexAndCotTheta, true>;

    using Base::Base;
  };

  Subset subset(const IndexSubset& subset) const noexcept {
    return Subset(*this, subset);
  }
  Subset2 subset(const IndexAndCotThetaSubset& subset) const noexcept {
    return Subset2(*this, subset);
  }

 private:
  std::vector<SpacePointIndex2> m_spacePoints;

  // parameters required to calculate a circle with linear equation
  std::vector<float> m_cotTheta;
  std::vector<std::array<float, 2>> m_er_iDeltaR;
  std::vector<std::array<float, 2>> m_uv;
  std::vector<std::array<float, 2>> m_xy;
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

class DoubletSeedFinder {
 public:
  struct Config {
    /// Direction of the doublet candidate space points. Either forward, also
    /// called top doublet, or backward, also called bottom doublet.
    Direction candidateDirection = Direction::Forward();

    /// Minimum radial distance between two doublet components
    float deltaRMin = 5 * UnitConstants::mm;
    /// Maximum radial distance between two doublet components
    float deltaRMax = 270 * UnitConstants::mm;

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
    Delegate<bool(const ConstSpacePointProxy2& /*middle*/,
                  const ConstSpacePointProxy2& /*other*/, float /*cotTheta*/,
                  bool /*isBottomCandidate*/)>
        experimentCuts;

    /// Whether the input space points are sorted by radius
    bool spacePointsSortedByRadius = false;
  };

  struct DerivedConfig : public Config {
    DerivedConfig(const Config& config, float bFieldInZ);

    float bFieldInZ = std::numeric_limits<float>::quiet_NaN();
    float minHelixDiameter2 = std::numeric_limits<float>::quiet_NaN();
  };

  static MiddleSpInfo computeMiddleSpInfo(const ConstSpacePointProxy2& spM);

  static std::unique_ptr<DoubletSeedFinder> create(const DerivedConfig& config);

  virtual ~DoubletSeedFinder() = default;

  virtual const DerivedConfig& config() const = 0;

  /// Creates compatible dublets by applying a series of cuts that can be
  /// tested with only two SPs.
  ///
  /// @param middleSp Space point candidate to be used as middle SP in a seed
  /// @param middleSpInfo Information about the middle space point
  /// @param candidateSps Subset of space points to be used as candidates for
  ///   middle SP in a seed
  /// @param compatibleDoublets Output container for compatible doublets
  virtual void createDoublets(
      const ConstSpacePointProxy2& middleSp, const MiddleSpInfo& middleSpInfo,
      SpacePointContainer2::ConstSubset& candidateSps,
      DoubletsForMiddleSp& compatibleDoublets) const = 0;

  /// Creates compatible dublets by applying a series of cuts that can be
  /// tested with only two SPs.
  ///
  /// @param middleSp Space point candidate to be used as middle SP in a seed
  /// @param middleSpInfo Information about the middle space point
  /// @param candidateSps Range of space points to be used as candidates for
  ///   middle SP in a seed
  /// @param compatibleDoublets Output container for compatible doublets
  virtual void createDoublets(
      const ConstSpacePointProxy2& middleSp, const MiddleSpInfo& middleSpInfo,
      SpacePointContainer2::ConstRange& candidateSps,
      DoubletsForMiddleSp& compatibleDoublets) const = 0;
};

}  // namespace Acts::Experimental
