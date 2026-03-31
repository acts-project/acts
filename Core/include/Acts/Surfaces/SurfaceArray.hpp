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
#include "Acts/Utilities/Diagnostics.hpp"
#include "Acts/Utilities/IAxis.hpp"

#include <iostream>
#include <vector>

namespace Acts {

using SurfaceVector = std::vector<const Surface*>;

/// @brief Provides Surface binning in 2 dimensions
///
/// Uses @c Grid under the hood to implement the storage and lookup
/// Contains a lookup struct which talks to the @c Grid
/// and performs utility actions. This struct needs to be initialised
/// externally and passed to @c SurfaceArray on construction.
class SurfaceArray {
 public:
  /// Base interface for all surface lookups.
  struct ISurfaceGridLookup {
    virtual ~ISurfaceGridLookup() = default;

    /// Fill provided surfaces into the contained @c Grid.
    /// @param gctx The current geometry context object, e.g. alignment
    /// @param surfaces Input surface pointers
    virtual void fill(const GeometryContext& gctx,
                      std::span<const Surface* const> surfaces) = 0;

    /// Performs lookup at @c pos and returns bin content as const reference
    /// @param position Lookup position
    /// @param direction Lookup direction
    /// @return A vector of surfaces at given bin
    [[deprecated("Use lookup with GeometryContext instead")]]
    virtual const std::vector<const Surface*>& lookup(
        const Vector3& position, const Vector3& direction) const = 0;

    /// Performs lookup at @c pos and returns bin content as const reference
    /// @param gctx The current geometry context object, e.g. alignment
    /// @param position Lookup position
    /// @param direction Lookup direction
    /// @return A span of surface pointers
    virtual std::span<const Surface* const> lookup(
        const GeometryContext& gctx, const Vector3& position,
        const Vector3& direction) const = 0;

    /// Performs lookup at global bin and returns bin content as reference
    /// @param bin Global lookup bin
    /// @return A vector of surfaces at given bin
    [[deprecated(
        "Mutable access will be removed without replacement in the future")]]
    virtual std::vector<const Surface*>& lookup(std::size_t bin) = 0;

    /// Performs lookup at global bin and returns bin content as const reference
    /// @param bin Global lookup bin
    /// @return A vector of surfaces at given bin
    virtual const std::vector<const Surface*>& lookup(
        std::size_t bin) const = 0;

    /// Performs a lookup at @c pos, but returns neighbors as well
    /// @param position Lookup position
    /// @param direction Lookup direction
    /// @return A vector of surfaces at given bin. Copy of all bins selected
    [[deprecated("Use neighbors with GeometryContext instead")]]
    virtual const std::vector<const Surface*>& neighbors(
        const Vector3& position, const Vector3& direction) const = 0;

    /// Performs a lookup at @c pos, but returns neighbors as well
    /// @param gctx The current geometry context object, e.g. alignment
    /// @param position Lookup position
    /// @param direction Lookup direction
    /// @return A span of surface pointers
    virtual std::span<const Surface* const> neighbors(
        const GeometryContext& gctx, const Vector3& position,
        const Vector3& direction) const = 0;

    /// Returns the total size of the grid (including under/overflow bins)
    /// @return Size of the grid data structure
    virtual std::size_t size() const = 0;

    /// Gets the center position of bin @c bin in global coordinates
    /// @param bin the global bin index
    /// @return The bin center
    virtual Vector3 getBinCenter(std::size_t bin) const = 0;

    /// Returns copies of the axes used in the grid as @c AnyAxis
    /// @return The axes
    /// @note This returns copies. Use for introspection and querying.
    virtual std::vector<const IAxis*> getAxes() const = 0;

    /// Get the representative surface used for this lookup
    /// @return Surface pointer
    virtual const Surface* surfaceRepresentation() const = 0;

    /// Checks if global bin is valid
    /// @param bin the global bin index
    /// @return bool if the bin is valid
    /// @note Valid means that the index points to a bin which is not a under
    ///       or overflow bin or out of range in any axis.
    virtual bool isValidBin(std::size_t bin) const = 0;

    /// The binning values described by this surface grid lookup. They are in
    /// order of the axes (optional) and empty for eingle lookups
    /// @return Vector of axis directions for binning
    virtual std::vector<AxisDirection> binningValues() const { return {}; }

    /// Get the number of local bins in each dimension. This is used to
    /// determine the size of the grid for neighbor lookups.
    /// @return Array of number of local bins in each dimension
    virtual std::array<std::size_t, 2> numLocalBins() const = 0;

    /// Get the maximum neighbor distance that is supported by this lookup. This
    /// is used to determine how many neighbors to include in neighbor lookups.
    /// @return Maximum neighbor distance
    virtual std::uint8_t maxNeighborDistance() const = 0;

    /// Get all surfaces in bin given by local grid indices and neighbor
    /// distance.
    /// @param gridIndices the local grid indices
    /// @param neighborDistance the neighbor distance to include in the lookup
    /// @return span of surface pointers of the bin at that position and its neighbors
    virtual std::span<const Surface* const> at(
        std::array<std::size_t, 2> gridIndices,
        std::uint8_t neighborDistance) const = 0;
  };

 private:
  /// Factory method to create a surface grid lookup for a given representative
  /// surface, tolerance, and axes. This will internally create the appropriate
  /// lookup class based on the axes and concrete @ref Grid.
  /// @param representative The surface which is used as representative
  /// @param tolerance The tolerance used for intersection checks
  /// @param axes The axes used for the grid
  /// @param maxNeighborDistance Maximum next neighbor distance to be included in neighbor lookups
  /// @return A unique pointer to the surface grid lookup
  static std::unique_ptr<ISurfaceGridLookup> makeSurfaceGridLookup(
      std::shared_ptr<RegularSurface> representative, double tolerance,
      std::tuple<const IAxis&, const IAxis&> axes,
      std::uint8_t maxNeighborDistance = 1);

  // This is temporary until Gen1 is removed
  friend class SurfaceArrayCreator;

 public:
  /// Lookup helper which encapsulates a @c Grid
  /// @tparam Axis1 The first axis
  /// @tparam Axis2 The second axis
  /// @deprecated This is deprecated in favor of direct @ref SurfaceArray constructors.
  template <class Axis1, class Axis2>
  struct [[deprecated("Use makeSurfaceGridLookup instead")]] SurfaceGridLookup
      : ISurfaceGridLookup {
    /// Construct a surface grid lookup
    /// @param representative The surface which is used as representative
    /// @param tolerance The tolerance used for intersection checks
    /// @param axes The axes used for the grid
    /// @param bValues Optional vector of axis directions for binning
    /// @param maxNeighborDistance Maximum next neighbor distance to be included in neighbor lookups
    SurfaceGridLookup(std::shared_ptr<RegularSurface> representative,
                      double tolerance, std::tuple<Axis1, Axis2> axes,
                      const std::vector<AxisDirection>& bValues = {},
                      std::uint8_t maxNeighborDistance = 1)
        : m_impl(makeSurfaceGridLookup(std::move(representative), tolerance,
                                       std::move(axes), maxNeighborDistance)) {
      static_cast<void>(bValues);
    }

    /// Fill provided surfaces into the contained @c Grid.
    /// @param gctx The current geometry context object, e.g. alignment
    /// @param surfaces Input surface pointers
    void fill(const GeometryContext& gctx,
              std::span<const Surface* const> surfaces) override {
      m_impl->fill(gctx, surfaces);
    }

    /// Performs lookup at @c pos and returns bin content as const reference
    /// @param position Lookup position
    /// @param direction Lookup direction
    /// @return A vector of surfaces at given bin
    const std::vector<const Surface*>& lookup(
        const Vector3& position, const Vector3& direction) const override {
      ACTS_PUSH_IGNORE_DEPRECATED()
      return m_impl->lookup(position, direction);
      ACTS_POP_IGNORE_DEPRECATED()
    }

    /// Performs lookup at @c pos and returns bin content as const reference
    /// @param gctx The current geometry context object, e.g. alignment
    /// @param position Lookup position
    /// @param direction Lookup direction
    /// @return A span of surface pointers
    std::span<const Surface* const> lookup(
        const GeometryContext& gctx, const Vector3& position,
        const Vector3& direction) const override {
      return m_impl->lookup(gctx, position, direction);
    }

    /// Performs lookup at global bin and returns bin content as reference
    /// @param bin Global lookup bin
    /// @return A vector of surfaces at given bin
    std::vector<const Surface*>& lookup(std::size_t bin) override {
      return m_impl->lookup(bin);
    }

    /// Performs lookup at global bin and returns bin content as const reference
    /// @param bin Global lookup bin
    /// @return A vector of surfaces at given bin
    const std::vector<const Surface*>& lookup(std::size_t bin) const override {
      return m_impl->lookup(bin);
    }

    /// Performs a lookup at @c pos, but returns neighbors as well
    /// @param position Lookup position
    /// @param direction Lookup direction
    /// @return A vector of surfaces at given bin. Copy of all bins selected
    const std::vector<const Surface*>& neighbors(
        const Vector3& position, const Vector3& direction) const override {
      ACTS_PUSH_IGNORE_DEPRECATED()
      return m_impl->neighbors(position, direction);
      ACTS_POP_IGNORE_DEPRECATED()
    }

    /// Performs a lookup at @c pos, but returns neighbors as well
    /// @param gctx The current geometry context object, e.g. alignment
    /// @param position Lookup position
    /// @param direction Lookup direction
    /// @return A span of surface pointers
    std::span<const Surface* const> neighbors(
        const GeometryContext& gctx, const Vector3& position,
        const Vector3& direction) const override {
      return m_impl->neighbors(gctx, position, direction);
    }

    /// Returns the total size of the grid (including under/overflow bins)
    /// @return Size of the grid data structure
    std::size_t size() const override { return m_impl->size(); }

    /// Gets the center position of bin @c bin in global coordinates
    /// @param bin the global bin index
    /// @return The bin center
    Vector3 getBinCenter(std::size_t bin) const override {
      return m_impl->getBinCenter(bin);
    }

    /// Returns copies of the axes used in the grid as @c AnyAxis
    /// @return The axes
    /// @note This returns copies. Use for introspection and querying.
    std::vector<const IAxis*> getAxes() const override {
      return m_impl->getAxes();
    }

    /// Get the representative surface used for this lookup
    /// @return Surface pointer
    const Surface* surfaceRepresentation() const override {
      return m_impl->surfaceRepresentation();
    }

    /// Checks if global bin is valid
    /// @param bin the global bin index
    /// @return bool if the bin is valid
    /// @note Valid means that the index points to a bin which is not a under
    ///       or overflow bin or out of range in any axis.
    bool isValidBin(std::size_t bin) const override {
      return m_impl->isValidBin(bin);
    }

    /// The binning values described by this surface grid lookup. They are in
    /// order of the axes (optional) and empty for eingle lookups
    /// @return Vector of axis directions for binning
    std::vector<AxisDirection> binningValues() const override {
      return m_impl->binningValues();
    }

    std::array<std::size_t, 2> numLocalBins() const override {
      return m_impl->numLocalBins();
    }

    std::uint8_t maxNeighborDistance() const override {
      return m_impl->maxNeighborDistance();
    }

    std::span<const Surface* const> at(
        std::array<std::size_t, 2> gridIndices,
        std::uint8_t neighborDistance) const override {
      return m_impl->at(gridIndices, neighborDistance);
    }

   private:
    std::unique_ptr<ISurfaceGridLookup> m_impl;
  };

  /// Lookup implementation which wraps one element and always returns this
  /// element when lookup is called
  /// @deprecated Construct the @ref SurfaceArray directly with a single surface
  struct [[deprecated(
      "Construct the SurfaceArray directly with a single "
      "surface")]] SingleElementLookup : ISurfaceGridLookup {
    /// Default constructor.
    /// @param element the one and only element.
    explicit SingleElementLookup(const Surface* element)
        : m_element({element}) {}

    /// Default constructor.
    /// @param elements the surfaces that are provided through a single lookup
    explicit SingleElementLookup(const std::vector<const Surface*>& elements)
        : m_element(elements) {}

    /// Lookup, always returns @c element
    /// @return reference to vector containing only @c element
    const std::vector<const Surface*>& lookup(
        const Vector3& /*position*/,
        const Vector3& /*direction*/) const override {
      return m_element;
    }

    /// Lookup, always returns @c element
    /// @return span of surface pointers containing only @c element
    std::span<const Surface* const> lookup(
        const GeometryContext& /*gctx*/, const Vector3& /*position*/,
        const Vector3& /*direction*/) const override {
      return m_element;
    }

    /// Lookup, always returns @c element
    /// @return reference to vector containing only @c element
    std::vector<const Surface*>& lookup(std::size_t /*bin*/) override {
      return m_element;
    }

    /// Lookup, always returns @c element
    /// @return reference to vector containing only @c element
    const std::vector<const Surface*>& lookup(
        std::size_t /*bin*/) const override {
      return m_element;
    }

    /// Lookup, always returns @c element
    /// @return reference to vector containing only @c element
    const std::vector<const Surface*>& neighbors(
        const Vector3& /*position*/,
        const Vector3& /*direction*/) const override {
      return m_element;
    }

    /// Lookup, always returns @c element
    /// @return span of surface pointers containing only @c element
    std::span<const Surface* const> neighbors(
        const GeometryContext& /*gctx*/, const Vector3& /*position*/,
        const Vector3& /*direction*/) const override {
      return m_element;
    }

    /// Returns the total size of the grid (including under/overflow bins)
    /// @return 1
    std::size_t size() const override { return 1; }

    /// Gets the bin center, but always returns (0, 0, 0)
    /// @return (0, 0, 0)
    Vector3 getBinCenter(std::size_t /*bin*/) const override {
      return Vector3(0, 0, 0);
    }

    /// Returns an empty vector of @c AnyAxis
    /// @return empty vector
    std::vector<const IAxis*> getAxes() const override { return {}; }

    const Surface* surfaceRepresentation() const override { return nullptr; }

    /// Comply with concept and provide fill method
    /// @note Does nothing
    void fill(const GeometryContext& /*gctx*/,
              std::span<const Surface* const> /*surfaces*/) override {}

    /// Returns if the bin is valid (it is)
    /// @return always true
    bool isValidBin(std::size_t /*bin*/) const override { return true; }

    std::array<std::size_t, 2> numLocalBins() const override { return {1, 1}; }

    std::uint8_t maxNeighborDistance() const override { return 0; }

    std::span<const Surface* const> at(
        std::array<std::size_t, 2> gridIndices,
        std::uint8_t neighborDistance) const override {
      if (gridIndices != std::array<std::size_t, 2>{0, 0} ||
          neighborDistance != 0) {
        throw std::out_of_range(
            "SingleElementLookupImpl only contains one bin with zero neighbor "
            "distance");
      }
      return m_element;
    }

   private:
    std::vector<const Surface*> m_element;
  };

  /// Default constructor which takes a @c SurfaceLookup and a vector of surfaces
  /// @param gridLookup The grid storage. @c SurfaceArray does not fill it on its own
  /// @param surfaces The input vector of surfaces. This is only for bookkeeping, so we can ask
  /// @param transform Optional additional transform for this SurfaceArray
  [[deprecated("Use the constructor with axes instead")]]
  SurfaceArray(std::unique_ptr<ISurfaceGridLookup> gridLookup,
               std::vector<std::shared_ptr<const Surface>> surfaces,
               const Transform3& transform = Transform3::Identity());

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
  SurfaceArray(const GeometryContext& gctx,
               std::vector<std::shared_ptr<const Surface>> surfaces,
               std::shared_ptr<RegularSurface> representative, double tolerance,
               std::tuple<const IAxis&, const IAxis&> axes);

  /// Get all surfaces in bin given by position @p pos.
  /// @param position the lookup position
  /// @param direction the lookup direction
  /// @return const reference to surface vector contained in bin at that position
  const std::vector<const Surface*>& at(const Vector3& position,
                                        const Vector3& direction) const {
    ACTS_PUSH_IGNORE_DEPRECATED()
    return m_gridLookup->lookup(position, direction);
    ACTS_POP_IGNORE_DEPRECATED()
  }

  /// Get all surfaces in bin given by position @p pos.
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position the lookup position
  /// @param direction the lookup direction
  /// @return span of surface pointers of the bin at that position
  std::span<const Surface* const> at(const GeometryContext& gctx,
                                     const Vector3& position,
                                     const Vector3& direction) const {
    return m_gridLookup->lookup(gctx, position, direction);
  }

  /// Get all surfaces in bin given by global bin index @p bin.
  /// @param bin the global bin index
  /// @return mutable reference to surface vector contained in bin
  [[deprecated(
      "Mutable access will be removed without replacement in the future")]]
  std::vector<const Surface*>& at(std::size_t bin) {
    ACTS_PUSH_IGNORE_DEPRECATED()
    return m_gridLookup->lookup(bin);
    ACTS_POP_IGNORE_DEPRECATED()
  }

  /// Get all surfaces in bin given by global bin index.
  /// @param bin the global bin index
  /// @return const reference to surface vector contained in bin
  const std::vector<const Surface*>& at(std::size_t bin) const {
    return const_cast<const ISurfaceGridLookup*>(m_gridLookup.get())
        ->lookup(bin);
  }

  /// Get all surfaces in bin at @p pos and its neighbors
  /// @param position The position to lookup
  /// @param direction The direction to lookup
  /// @return Merged surface vector of neighbors and nominal
  [[deprecated("Use neighbors with GeometryContext instead")]]
  const std::vector<const Surface*>& neighbors(const Vector3& position,
                                               const Vector3& direction) const {
    ACTS_PUSH_IGNORE_DEPRECATED()
    return m_gridLookup->neighbors(position, direction);
    ACTS_POP_IGNORE_DEPRECATED()
  }

  /// Get all surfaces in bin at @p pos and its neighbors
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The position to lookup
  /// @param direction The direction to lookup
  /// @return span of surface pointers of neighbors and nominal
  std::span<const Surface* const> neighbors(const GeometryContext& gctx,
                                            const Vector3& position,
                                            const Vector3& direction) const {
    return m_gridLookup->neighbors(gctx, position, direction);
  }

  /// Get the size of the underlying grid structure including under/overflow
  /// bins
  /// @return the size
  std::size_t size() const { return m_gridLookup->size(); }

  /// Get the center of the bin identified by global bin index @p bin
  /// @param bin the global bin index
  /// @return Center position of the bin in global coordinates
  Vector3 getBinCenter(std::size_t bin) const {
    return m_gridLookup->getBinCenter(bin);
  }

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
  std::vector<const IAxis*> getAxes() const { return m_gridLookup->getAxes(); }

  /// Checks if global bin is valid
  /// @param bin the global bin index
  /// @return bool if the bin is valid
  /// @note Valid means that the index points to a bin which is not a under
  ///       or overflow bin or out of range in any axis.
  bool isValidBin(std::size_t bin) const {
    return m_gridLookup->isValidBin(bin);
  }

  /// Get the transform of this surface array.
  /// @return Reference to the transformation matrix
  /// @deprecated This is an implementation detail and will be removed soon
  [[deprecated("This is an implementation detail and will be removed soon")]]
  const Transform3& transform() const {
    return m_transform;
  }

  /// The binning values described by this surface grid lookup. They are in
  /// order of the axes
  /// @return Vector of axis directions for binning
  std::vector<AxisDirection> binningValues() const {
    return m_gridLookup->binningValues();
  }

  /// Get string representation of this @c SurfaceArray
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param sl Output stream to write to
  /// @return the output stream given as @p sl
  std::ostream& toStream(const GeometryContext& gctx, std::ostream& sl) const;

  /// Return the lookup object
  /// @return Reference to the surface grid lookup interface
  [[deprecated(
      "Grid lookup is an implementation detail and will be removed soon")]]
  const ISurfaceGridLookup& gridLookup() const {
    return *m_gridLookup;
  }

  /// Get the representative surface used for this surface array
  /// @return Surface pointer
  const Surface* surfaceRepresentation() const {
    return m_gridLookup->surfaceRepresentation();
  }

  /// Get the number of local bins in each dimension. This is used to
  /// determine the size of the grid for neighbor lookups.
  /// @return Array of number of local bins in each dimension
  std::array<std::size_t, 2> numLocalBins() const {
    return m_gridLookup->numLocalBins();
  }

  /// Get the maximum neighbor distance that is supported by this lookup. This
  /// is used to determine how many neighbors to include in neighbor lookups.
  /// @return Maximum neighbor distance
  std::uint8_t maxNeighborDistance() const {
    return m_gridLookup->maxNeighborDistance();
  }

  /// Get all surfaces in bin given by local grid indices and neighbor
  /// distance.
  /// @param gridIndices the local grid indices
  /// @param neighborDistance the neighbor distance to include in the lookup
  /// @return span of surface pointers of the bin at that position and its neighbors
  std::span<const Surface* const> at(std::array<std::size_t, 2> gridIndices,
                                     std::uint8_t neighborDistance) const {
    return m_gridLookup->at(gridIndices, neighborDistance);
  }

 private:
  /// The actual grid lookup implementation
  std::unique_ptr<ISurfaceGridLookup> m_gridLookup;
  /// this vector makes sure we have shared ownership over the surfaces
  std::vector<std::shared_ptr<const Surface>> m_surfaces;
  /// this vector is returned, so that (expensive) copying of the shared_ptr
  /// vector does not happen by default
  std::vector<const Surface*> m_surfacesRawPointers;
  /// this is only used to keep info on transform applied by l2g and g2l
  Transform3 m_transform;
};

}  // namespace Acts
