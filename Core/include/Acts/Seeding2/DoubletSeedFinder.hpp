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
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/detail/ContainerIterator.hpp"

#include <cstdint>
#include <vector>

namespace Acts {

/// Container for doublets found by the doublet seed finder.
///
/// This implementation uses partial AoS/SoA depending on the access pattern in
/// the doublet finding process.
class DoubletsForMiddleSp {
 public:
  /// Type alias for index type used in doublets container
  using Index = std::uint32_t;
  /// Type alias for range of indices in doublets container
  using IndexRange = std::pair<Index, Index>;
  /// Type alias for subset of indices in doublets container
  using IndexSubset = std::span<const Index>;

  /// Check if the doublets container is empty
  /// @return True if container has no doublets
  [[nodiscard]] bool empty() const { return m_spacePoints.empty(); }
  /// Get the number of doublets in container
  /// @return Number of doublets stored
  [[nodiscard]] Index size() const {
    return static_cast<Index>(m_spacePoints.size());
  }

  /// Clear all stored doublets and associated data
  void clear() {
    m_spacePoints.clear();
    m_cotTheta.clear();
    m_er_iDeltaR.clear();
    m_uv.clear();
    m_xy.clear();
  }

  /// Add a new doublet with associated parameters
  /// @param sp spacepoint index for the doublet
  /// @param cotTheta Cotangent of polar angle
  /// @param iDeltaR Inverse delta R parameter
  /// @param er Error in R coordinate
  /// @param u U coordinate parameter
  /// @param v V coordinate parameter
  /// @param x X coordinate
  /// @param y Y coordinate
  void emplace_back(SpacePointIndex2 sp, float cotTheta, float iDeltaR,
                    float er, float u, float v, float x, float y) {
    m_spacePoints.push_back(sp);
    m_cotTheta.push_back(cotTheta);
    m_er_iDeltaR.push_back({er, iDeltaR});
    m_uv.push_back({u, v});
    m_xy.push_back({x, y});
  }

  /// Get reference to spacepoint indices container
  /// @return Const reference to spacepoint indices vector
  const std::vector<SpacePointIndex2>& spacePoints() const {
    return m_spacePoints;
  }
  /// Get reference to cotTheta values container
  /// @return Const reference to cotTheta values vector
  const std::vector<float>& cotTheta() const { return m_cotTheta; }

  /// Pair of doublet index and cotTheta value.
  struct IndexAndCotTheta {
    /// Doublet index
    Index index{};
    /// Cotangent of theta
    float cotTheta{};
  };

  /// Type alias for subset of index and cotTheta pairs
  using IndexAndCotThetaSubset = std::span<const IndexAndCotTheta>;

  /// Sort doublets by cotTheta within given range
  /// @param range Index range to sort within
  /// @param indexAndCotTheta Output vector containing sorted index and cotTheta pairs
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

  /// Proxy accessor for a single doublet entry.
  class Proxy {
   public:
    /// Constructor
    /// @param container Pointer to the doublet container
    /// @param index Index of the doublet
    Proxy(const DoubletsForMiddleSp* container, Index index)
        : m_container(container), m_index(index) {}

    /// Get the doublet container
    /// @return Reference to the container
    const DoubletsForMiddleSp& container() const { return *m_container; }
    /// Get the doublet index
    /// @return The index
    Index index() const { return m_index; }

    /// Get spacepoint index pair
    /// @return The spacepoint index
    SpacePointIndex2 spacePointIndex() const {
      return m_container->m_spacePoints[m_index];
    }

    /// Get cotangent of theta
    /// @return The cotTheta value
    float cotTheta() const { return m_container->m_cotTheta[m_index]; }
    /// Get er value
    /// @return The er value
    float er() const { return m_container->m_er_iDeltaR[m_index][0]; }
    /// Get inverse delta r
    /// @return The inverse delta r value
    float iDeltaR() const { return m_container->m_er_iDeltaR[m_index][1]; }
    /// Get u coordinate
    /// @return The u value
    float u() const { return m_container->m_uv[m_index][0]; }
    /// Get v coordinate
    /// @return The v value
    float v() const { return m_container->m_uv[m_index][1]; }
    /// Get x coordinate
    /// @return The x value
    float x() const { return m_container->m_xy[m_index][0]; }
    /// Get y coordinate
    /// @return The y value
    float y() const { return m_container->m_xy[m_index][1]; }

   private:
    const DoubletsForMiddleSp* m_container{};
    Index m_index{};
  };
  /// Same as `Proxy` but also contains `cotTheta`. This is useful after sorting
  /// doublets by `cotTheta` to avoid indirect access.
  class Proxy2 : public Proxy {
   public:
    /// Constructor for Proxy2 with precomputed cotTheta
    /// @param container Pointer to the doublets container
    /// @param indexAndCotTheta Index and cotTheta pair
    Proxy2(const DoubletsForMiddleSp* container,
           IndexAndCotTheta indexAndCotTheta)
        : Proxy(container, indexAndCotTheta.index),
          m_cotTheta(indexAndCotTheta.cotTheta) {}

    /// Get precomputed cotTheta value (avoids indirect access)
    /// @return Cotangent of polar angle
    float cotTheta() const { return m_cotTheta; }

   private:
    float m_cotTheta{};
  };

  /// Access doublet by index
  /// @param index Index of the doublet to access
  /// @return Proxy object for the doublet
  Proxy operator[](Index index) const { return Proxy(this, index); }
  /// Access doublet by index and cotTheta pair
  /// @param indexAndCotTheta Index and cotTheta pair for the doublet
  /// @return Proxy2 object for the doublet with precomputed cotTheta
  Proxy2 operator[](IndexAndCotTheta indexAndCotTheta) const {
    return Proxy2(this, indexAndCotTheta);
  }

  /// Type alias for const iterator over doublets in container
  using const_iterator =
      detail::ContainerIterator<DoubletsForMiddleSp, Proxy, Index, true>;

  /// Get iterator to beginning of doublets container
  /// @return Const iterator to first doublet
  const_iterator begin() const { return const_iterator(*this, 0); }
  /// Get iterator to end of doublets container
  /// @return Const iterator past the last doublet
  const_iterator end() const { return const_iterator(*this, size()); }

  /// Range view over doublets in the container.
  class Range : public detail::ContainerRange<Range, Range, DoubletsForMiddleSp,
                                              Index, true> {
   public:
    /// Base class type alias
    using Base =
        detail::ContainerRange<Range, Range, DoubletsForMiddleSp, Index, true>;

    using Base::Base;
  };

  /// Get range view of all doublets
  /// @return Range object covering all doublets
  Range range() const noexcept { return Range(*this, {0, size()}); }
  /// Get range view of doublets within specified index range
  /// @param range Index range to create view for
  /// @return Range object covering specified doublets
  Range range(const IndexRange& range) const noexcept {
    return Range(*this, range);
  }

  /// Subset view of doublets addressed by indices.
  class Subset
      : public detail::ContainerSubset<Subset, Subset, DoubletsForMiddleSp,
                                       Proxy, Index, true> {
   public:
    /// Base class type alias
    using Base = detail::ContainerSubset<Subset, Subset, DoubletsForMiddleSp,
                                         Proxy, Index, true>;

    using Base::Base;
  };
  /// Subset view of doublets addressed by index and cotTheta pairs.
  class Subset2
      : public detail::ContainerSubset<Subset2, Subset2, DoubletsForMiddleSp,
                                       Proxy2, IndexAndCotTheta, true> {
   public:
    /// Base class type alias
    using Base = detail::ContainerSubset<Subset2, Subset2, DoubletsForMiddleSp,
                                         Proxy2, IndexAndCotTheta, true>;

    using Base::Base;
  };

  /// Create subset view from index subset
  /// @param subset Span of indices to include in subset
  /// @return Subset object for the specified indices
  Subset subset(const IndexSubset& subset) const noexcept {
    return Subset(*this, subset);
  }
  /// Create subset view from index and cotTheta subset
  /// @param subset Span of index and cotTheta pairs to include
  /// @return Subset2 object with precomputed cotTheta values
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

/// Derived quantities for the middle spacepoint in a doublet.
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

/// Interface and a collection of standard implementations for a doublet seed
/// finder. Given a starting spacepoint and a collection of candidates, it
/// finds all doublets that satisfy the selection criteria. For the standard
/// implementations the criteria are given by interaction point cuts.
///
/// @note The standard implementations rely on virtual function dispatch which
/// did not turn out to affect the performance after measurement.
class DoubletSeedFinder {
 public:
  /// Collection of configuration parameters for the doublet seed finder. This
  /// includes doublet cuts, steering switches, and assumptions about the space
  /// points.
  struct Config {
    /// Whether the input spacepoints are sorted by radius
    bool spacePointsSortedByRadius = false;

    /// Direction of the doublet candidate spacepoints. Either forward, also
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
    /// Maximum collision region boundary in z-axis for doublet origin checks
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

    /// Type alias for delegate to apply experiment specific cuts during doublet
    /// finding
    using ExperimentCuts =
        Delegate<bool(const ConstSpacePointProxy2& /*middle*/,
                      const ConstSpacePointProxy2& /*other*/,
                      float /*cotTheta*/, bool /*isBottomCandidate*/)>;

    /// Delegate to apply experiment specific cuts during doublet finding
    ExperimentCuts experimentCuts;
  };

  /// Derived configuration for the doublet seed finder using a magnetic field.
  struct DerivedConfig : public Config {
    /// Constructor for derived configuration with magnetic field
    /// @param config Base configuration to derive from
    /// @param bFieldInZ Magnetic field strength in z-direction
    DerivedConfig(const Config& config, float bFieldInZ);

    /// Magnetic field strength in z-direction for helix calculation
    float bFieldInZ = std::numeric_limits<float>::quiet_NaN();
    /// Squared minimum helix diameter derived from magnetic field and minimum
    /// pT
    float minHelixDiameter2 = std::numeric_limits<float>::quiet_NaN();
  };

  /// Computes additional quantities from the middle spacepoint which can be
  /// reused during doublet finding.
  /// @param spM Middle spacepoint for doublet computation
  /// @return MiddleSpInfo structure with computed quantities
  static MiddleSpInfo computeMiddleSpInfo(const ConstSpacePointProxy2& spM);

  /// Creates a new doublet seed finder instance given the configuration.
  /// @param config Configuration for the doublet seed finder
  /// @return Unique pointer to new DoubletSeedFinder instance
  static std::unique_ptr<DoubletSeedFinder> create(const DerivedConfig& config);

  virtual ~DoubletSeedFinder() = default;

  /// Returns the configuration of the doublet seed finder.
  /// @return Reference to the configuration object
  virtual const DerivedConfig& config() const = 0;

  /// Creates compatible dublets by applying a series of cuts that can be
  /// tested with only two SPs.
  ///
  /// @param middleSp spacepoint candidate to be used as middle SP in a seed
  /// @param middleSpInfo Information about the middle spacepoint
  /// @param candidateSps Subset of spacepoints to be used as candidates for
  ///   middle SP in a seed
  /// @param compatibleDoublets Output container for compatible doublets
  virtual void createDoublets(
      const ConstSpacePointProxy2& middleSp, const MiddleSpInfo& middleSpInfo,
      SpacePointContainer2::ConstSubset& candidateSps,
      DoubletsForMiddleSp& compatibleDoublets) const = 0;

  /// Creates compatible dublets by applying a series of cuts that can be
  /// tested with only two SPs.
  ///
  /// @param middleSp spacepoint candidate to be used as middle SP in a seed
  /// @param middleSpInfo Information about the middle spacepoint
  /// @param candidateSps Range of spacepoints to be used as candidates for
  ///   middle SP in a seed
  /// @param compatibleDoublets Output container for compatible doublets
  virtual void createDoublets(
      const ConstSpacePointProxy2& middleSp, const MiddleSpInfo& middleSpInfo,
      SpacePointContainer2::ConstRange& candidateSps,
      DoubletsForMiddleSp& compatibleDoublets) const = 0;
};

}  // namespace Acts
