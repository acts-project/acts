// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/IAxis.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

// Forward declared friend class for private access in tests.
namespace ActsTests {
struct SurfaceArrayCreatorFixture;
}

namespace Acts {

/// @brief Provides Surface binning in 2 dimensions
///
/// Uses @c Grid under the hood to implement the storage and lookup
/// Contains a lookup struct which talks to the @c Grid
/// and performs utility actions. This struct needs to be initialised
/// externally and passed to @c SurfaceArray on construction.
class SurfaceArray {
 public:
  /// Bounds on the neighbor window the lookup derives from the crossing angle,
  /// in bins per grid axis. The floor is for lookups that cannot see the
  /// surfaces the fill saw, e.g. under an aligned geometry context.
  struct NeighborWindow {
    /// Fewest bins served in each direction, regardless of the crossing angle
    std::array<std::uint8_t, 2> min = {0, 0};
    /// Most bins served in each direction; also sizes the neighbor cache
    std::array<std::uint8_t, 2> max = {2, 2};

    /// Default window: no floor, at most two bins per axis
    constexpr NeighborWindow() = default;

    /// Window with an explicit floor and bound per axis
    /// @param windowMin Fewest bins served in each direction
    /// @param windowMax Most bins served in each direction
    constexpr NeighborWindow(std::array<std::uint8_t, 2> windowMin,
                             std::array<std::uint8_t, 2> windowMax)
        : min(windowMin), max(windowMax) {}

    /// Isotropic window bounded by a single neighbor distance, so that a call
    /// written against the scalar interface keeps the window it asked for: at
    /// least one bin unless the bound itself forbids it, at most
    /// @p maxNeighborDistance.
    /// @deprecated Give the bounds per axis, the two axes of a layer grid are
    ///             not equivalent. See @ref NeighborWindow.
    /// @param maxNeighborDistance Most bins served in each direction
    [[deprecated(
        "Pass an Acts::SurfaceArray::NeighborWindow with per-axis bounds "
        "instead of a single neighbor distance")]]
    constexpr NeighborWindow(  // NOLINT(google-explicit-constructor)
        std::uint8_t maxNeighborDistance)
        : min({std::min<std::uint8_t>(1, maxNeighborDistance),
               std::min<std::uint8_t>(1, maxNeighborDistance)}),
          max({maxNeighborDistance, maxNeighborDistance}) {}
  };

  /// Constructor with a single surface
  /// @param srf The one and only surface
  explicit SurfaceArray(std::shared_ptr<const Surface> srf);

  /// Constructor to create a surface grid lookup for a given representative
  /// surface, tolerance, and axes.
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param surfaces The input vector of surfaces that will be accessible
  ///                 through this @ref SurfaceArray.
  /// @param representative The surface which is used as representative
  /// @param tolerance The tolerance used for intersection checks
  /// @param axes The axes used for the grid
  /// @param neighborWindow Bounds on the neighbor window the lookup derives
  ///        from the crossing angle, in bins per grid axis
  SurfaceArray(const GeometryContext& gctx,
               std::vector<std::shared_ptr<const Surface>> surfaces,
               std::shared_ptr<RegularSurface> representative, double tolerance,
               std::tuple<const IAxis&, const IAxis&> axes,
               NeighborWindow neighborWindow = NeighborWindow{{0, 0}, {2, 2}});

  // non-copyable but movable due to unique_ptr member. deferred implementation
  // to source since the pimpl is not fully defined in the header.

  /// @param other the other SurfaceArray to copy from
  SurfaceArray(const SurfaceArray& other) = delete;

  /// @param other the other SurfaceArray to move from
  SurfaceArray(SurfaceArray&& other) noexcept;

  /// @param other the other SurfaceArray to copy-assign from
  /// @return reference to this SurfaceArray
  SurfaceArray& operator=(const SurfaceArray& other) = delete;

  /// @param other the other SurfaceArray to move-assign from
  /// @return reference to this SurfaceArray
  SurfaceArray& operator=(SurfaceArray&& other) noexcept;

  ~SurfaceArray();

  /// Get all surfaces in bin given by the global bin index
  /// @param bin the global bin index
  /// @return span of surface pointers of the bin at that position
  std::span<const Surface* const> at(std::size_t bin) const;

  /// Get all surfaces in bin given by position @p pos.
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position the lookup position
  /// @param direction the lookup direction
  /// @return span of surface pointers of the bin at that position
  std::span<const Surface* const> at(const GeometryContext& gctx,
                                     const Vector3& position,
                                     const Vector3& direction) const;

  /// Get all surfaces in bin given by local grid indices and neighbor
  /// distance.
  /// @param gridIndices the local grid indices
  /// @param neighborDistance the neighbor distance to include in the lookup,
  ///        the same along both axes
  /// @return span of surface pointers of the bin at that position and its neighbors
  std::span<const Surface* const> neighbors(
      std::array<std::size_t, 2> gridIndices,
      std::uint8_t neighborDistance) const;

  /// Get all surfaces in bin given by local grid indices and a neighbor
  /// distance per axis.
  /// @param gridIndices the local grid indices
  /// @param neighborDistance the neighbor distance to include in the lookup,
  ///        per axis
  /// @return span of surface pointers of the bin at that position and its neighbors
  std::span<const Surface* const> neighbors(
      std::array<std::size_t, 2> gridIndices,
      std::array<std::uint8_t, 2> neighborDistance) const;

  /// Get all surfaces in bin at @p pos and its neighbors
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The position to lookup
  /// @param direction The direction to lookup
  /// @return span of surface pointers of neighbors and nominal
  std::span<const Surface* const> neighbors(const GeometryContext& gctx,
                                            const Vector3& position,
                                            const Vector3& direction) const;

  /// Get the size of the underlying grid structure including under/overflow
  /// bins
  /// @return the size
  std::size_t size() const;

  /// Get the center of the bin identified by global bin index @p bin
  /// @param bin the global bin index
  /// @return Center position of the bin in global coordinates
  Vector3 getBinCenter(std::size_t bin) const;

  /// Get all surfaces attached to this @c SurfaceArray
  /// @return Reference to vector of all surfaces
  /// @note This does not reflect the actual state of the grid. It only
  ///       returns what was given in the constructor, without any checks
  ///       if that is actually what's in the grid.
  const std::vector<const Surface*>& surfaces() const {
    return m_surfacesRawPointers;
  }

  /// Get vector of axes spanning the grid as @c AnyAxis
  /// @return vector of @c AnyAxis
  /// @note The axes in the vector are copies. Only use for introspection and
  ///       querying.
  std::vector<const IAxis*> getAxes() const;

  /// Checks if global bin is valid
  /// @param bin the global bin index
  /// @return bool if the bin is valid
  /// @note Valid means that the index points to a bin which is not a under
  ///       or overflow bin or out of range in any axis.
  bool isValidBin(std::size_t bin) const;

  /// The binning values described by this surface grid lookup. They are in
  /// order of the axes
  /// @return Vector of axis directions for binning
  std::vector<AxisDirection> binningValues() const;

  /// Get string representation of this @c SurfaceArray
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param sl Output stream to write to
  /// @return the output stream given as @p sl
  std::ostream& toStream(const GeometryContext& gctx, std::ostream& sl) const;

  /// Get the representative surface used for this surface array
  /// @return Surface pointer
  const Surface* surfaceRepresentation() const;

  /// Get the number of local bins in each dimension. This is used to
  /// determine the size of the grid for neighbor lookups.
  /// @return Array of number of local bins in each dimension
  std::array<std::size_t, 2> numLocalBins() const;

  /// Get the bounds on the neighbor window this lookup serves. The window
  /// itself is derived per axis from the crossing angle and clamped to them.
  /// @return Neighbor window bounds per axis
  NeighborWindow neighborWindow() const;

  /// Get the maximum neighbor distance that is supported by this lookup.
  /// @deprecated The window is bounded per axis, use @ref neighborWindow.
  ///             This reports the larger of the two bounds.
  /// @return Maximum neighbor distance over both axes
  [[deprecated(
      "Use Acts::SurfaceArray::neighborWindow(), the window is bounded per "
      "axis")]]
  std::uint8_t maxNeighborDistance() const;

  /// Forward declaration of the internal lookup struct. The actual definition
  /// is in the source file.
  struct ISurfaceGridLookup;

 private:
  /// The actual grid lookup implementation
  std::unique_ptr<ISurfaceGridLookup> m_gridLookup;
  /// this vector makes sure we have shared ownership over the surfaces
  std::vector<std::shared_ptr<const Surface>> m_surfaces;
  /// this vector is returned, so that (expensive) copying of the shared_ptr
  /// vector does not happen by default
  std::vector<const Surface*> m_surfacesRawPointers;
};

}  // namespace Acts
