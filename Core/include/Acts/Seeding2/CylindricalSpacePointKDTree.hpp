// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Utilities/KDTree.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <vector>

namespace Acts::Experimental {

/// A cylindrical space point KD-tree used for seeding in a cylindrical detector
/// geometry.
/// The tree is defined in cylindrical coordinates (phi, r, z) and allows for
/// efficient access to space points based on their azimuthal angle,
/// radial distance, and z-coordinate.
class CylindricalSpacePointKDTree {
 public:
  /// Space point index type used in the grid.
  using SpacePointIndex = std::uint32_t;

  /// @brief Set the number of dimensions in which to embed points. This is just
  /// 3 for now (phi, r, and z), but we might want to increase or decrease this
  /// number in the future.
  static constexpr std::size_t NDims = 3;

  /// @brief Enumeration of the different dimensions in which we can apply cuts.
  enum Dim { DimPhi = 0, DimR = 1, DimZ = 2 };

  /// @brief The k-d tree type used by this seeder internally, which is
  /// three-dimensional, contains internal spacepoint pointers, uses the Acts
  /// scalar type for coordinates, stores its coordinates in std::arrays, and
  /// has leaf size 4.
  using Tree = KDTree<NDims, SpacePointIndex, float, std::array, 4>;

  struct Options {
    /// maximum extension of sensitive detector layer relevant for seeding as
    /// distance from x=y=0 (i.e. in r)
    float rMax = 600 * UnitConstants::mm;
    /// minimum extension of sensitive detector layer relevant for seeding in
    /// negative direction in z
    float zMin = -2800 * UnitConstants::mm;
    /// maximum extension of sensitive detector layer relevant for seeding in
    /// positive direction in z
    float zMax = 2800 * UnitConstants::mm;
    /// minimum phi value for phiAxis construction
    float phiMin = -std::numbers::pi_v<float>;
    /// maximum phi value for phiAxis construction
    float phiMax = std::numbers::pi_v<float>;

    /// Minimum radial distance between two doublet components
    float deltaRMin = 5 * UnitConstants::mm;
    /// Maximum radial distance between two doublet components
    float deltaRMax = 270 * UnitConstants::mm;

    /// Minimal z distance between two doublet components
    float deltaZMin = -std::numeric_limits<float>::infinity();
    /// Maximum z distance between two doublet components
    float deltaZMax = std::numeric_limits<float>::infinity();

    /// Limiting location of collision region in z-axis used to check if doublet
    /// origin is within reasonable bounds
    float collisionRegionMin = -150 * UnitConstants::mm;
    float collisionRegionMax = +150 * UnitConstants::mm;

    /// Maximum allowed cotTheta between two space-points in doublet, used to
    /// check if forward angle is within bounds
    float cotThetaMax = 10.01788;  // equivalent to eta = 3 (pseudorapidity)

    /// Shrink the phi range of middle space-point (analogous to phi bin size in
    /// grid from default seeding + number of phi bins used in search)
    float deltaPhiMax = 0.085;
  };

  struct Candidates {
    /// denotes the candidates bottom seed points, assuming that the track has
    /// monotonically _increasing_ z position
    std::vector<SpacePointIndex> bottom_lh_v;
    /// denotes the candidate bottom points assuming that the track has
    /// monotonically _decreasing_ z position
    std::vector<SpacePointIndex> bottom_hl_v;
    /// are the candidate top points for an increasing z track
    std::vector<SpacePointIndex> top_lh_v;
    /// are the candidate top points for a decreasing z track
    std::vector<SpacePointIndex> top_hl_v;

    void reserve(std::size_t n) {
      bottom_lh_v.reserve(n);
      bottom_hl_v.reserve(n);
      top_lh_v.reserve(n);
      top_hl_v.reserve(n);
    }

    void clear() {
      bottom_lh_v.clear();
      bottom_hl_v.clear();
      top_lh_v.clear();
      top_hl_v.clear();
    }
  };

  /// Construct a cylindrical space point grid with the given configuration and
  /// an optional logger.
  explicit CylindricalSpacePointKDTree(
      Tree tree, std::unique_ptr<const Logger> logger = getDefaultLogger(
                     "CylindricalSpacePointKDTree", Logging::Level::INFO));

  std::size_t size() const { return m_tree.size(); }

  auto begin() const { return m_tree.begin(); }
  auto end() const { return m_tree.end(); }

  Tree::range_t validTupleOrthoRangeLH(const Options& options,
                                       const ConstSpacePointProxy2& low) const;
  Tree::range_t validTupleOrthoRangeHL(const Options& options,
                                       const ConstSpacePointProxy2& high) const;

  void validTuples(const Options& lhOptions, const Options& hlOptions,
                   const ConstSpacePointProxy2& spM, std::size_t nTopSeedConf,
                   Candidates& candidates) const;

 private:
  Tree m_tree;

  std::unique_ptr<const Logger> m_logger;

  const Logger& logger() const { return *m_logger; }
};

class CylindricalSpacePointKDTreeBuilder {
 public:
  /// Space point index type used in the grid.
  using SpacePointIndex = CylindricalSpacePointKDTree::SpacePointIndex;

  /// @brief Set the number of dimensions in which to embed points. This is just
  /// 3 for now (phi, r, and z), but we might want to increase or decrease
  /// this number in the future.
  static constexpr std::size_t NDims = CylindricalSpacePointKDTree::NDims;

  /// @brief Enumeration of the different dimensions in which we can apply cuts.
  using Dim = CylindricalSpacePointKDTree::Dim;

  /// @brief The k-d tree type used by this seeder internally, which is
  /// three-dimensional, contains internal spacepoint pointers, uses the Acts
  /// scalar type for coordinates, stores its coordinates in std::arrays, and
  /// has leaf size 4.
  using Tree = CylindricalSpacePointKDTree::Tree;

  /// Construct a cylindrical space point grid with the given configuration
  /// and an optional logger.
  explicit CylindricalSpacePointKDTreeBuilder(
      std::unique_ptr<const Logger> logger = getDefaultLogger(
          "CylindricalSpacePointKDTree", Logging::Level::INFO));

  /// Get the number of space points in the grid.
  /// @return The number of space points in the grid
  std::size_t size() const { return m_points.size(); }

  void reserve(std::size_t n) { m_points.reserve(n); }

  /// Clear the grid and drop all state. The object will behave like a newly
  /// constructed one.
  void clear() { m_points.clear(); }

  /// Insert a space point into the grid.
  /// @param index The index of the space point to insert
  /// @param phi The azimuthal angle of the space point in radians
  /// @param r The radial distance of the space point from the origin
  /// @param z The z-coordinate of the space point
  void insert(SpacePointIndex index, float phi, float r, float z);

  /// Insert a space point into the grid.
  /// @param sp The space point to insert
  void insert(const ConstSpacePointProxy2& sp) {
    return insert(sp.index(), sp.phi(), sp.r(), sp.z());
  }

  /// Fill the grid with space points from the container.
  /// @param spacePoints The space point container to fill the grid with
  void extend(const SpacePointContainer2::ConstRange& spacePoints);

  CylindricalSpacePointKDTree build();

 private:
  std::unique_ptr<const Logger> m_logger;

  std::vector<Tree::pair_t> m_points;

  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts::Experimental
