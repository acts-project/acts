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
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AnyGridView.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <iostream>
#include <limits>
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
  /// @brief Base interface for all surface lookups.
  struct ISurfaceGridLookup {
    /// @brief Fill provided surfaces into the contained @c Grid.
    /// @param gctx The current geometry context object, e.g. alignment
    /// @param surfaces Input surface pointers
    virtual void fill(const GeometryContext& gctx,
                      const SurfaceVector& surfaces) = 0;

    /// @brief Performs lookup at @c pos and returns bin content as const
    /// reference
    /// @param position Lookup position
    /// @param direction Lookup direction
    /// @return @c SurfaceVector at given bin
    virtual const SurfaceVector& lookup(const Vector3& position,
                                        const Vector3& direction) const = 0;

    /// @brief Performs lookup at global bin and returns bin content as
    /// reference
    /// @param bin Global lookup bin
    /// @return @c SurfaceVector at given bin
    virtual SurfaceVector& lookup(std::size_t bin) = 0;

    /// @brief Performs lookup at global bin and returns bin content as const
    /// reference
    /// @param bin Global lookup bin
    /// @return @c SurfaceVector at given bin
    virtual const SurfaceVector& lookup(std::size_t bin) const = 0;

    /// @brief Performs a lookup at @c pos, but returns neighbors as well
    ///
    /// @param position Lookup position
    /// @param direction Lookup direction
    /// @return @c SurfaceVector at given bin. Copy of all bins selected
    virtual const SurfaceVector& neighbors(const Vector3& position,
                                           const Vector3& direction) const = 0;

    /// @brief Returns the total size of the grid (including under/overflow
    /// bins)
    /// @return Size of the grid data structure
    virtual std::size_t size() const = 0;

    /// @brief Gets the center position of bin @c bin in global coordinates
    /// @param bin the global bin index
    /// @return The bin center
    virtual Vector3 getBinCenter(std::size_t bin) const = 0;

    /// @brief Returns copies of the axes used in the grid as @c AnyAxis
    /// @return The axes
    /// @note This returns copies. Use for introspection and querying.
    virtual std::vector<const IAxis*> getAxes() const = 0;

    /// @brief Get a view of the grid for inspection
    /// @return Optional grid view containing surface vectors
    virtual std::optional<AnyGridConstView<SurfaceVector>> getGridView()
        const = 0;

    /// @brief Get the representative surface used for this lookup
    /// @return Surface pointer
    virtual const Surface* surfaceRepresentation() const = 0;

    /// @brief Checks if global bin is valid
    /// @param bin the global bin index
    /// @return bool if the bin is valid
    /// @note Valid means that the index points to a bin which is not a under
    ///       or overflow bin or out of range in any axis.
    virtual bool isValidBin(std::size_t bin) const = 0;

    /// @brief The binning values described by this surface grid lookup
    /// They are in order of the axes (optional) and empty for eingle lookups
    /// @return Vector of axis directions for binning
    virtual std::vector<AxisDirection> binningValues() const { return {}; };

    /// Pure virtual destructor
    virtual ~ISurfaceGridLookup() = 0;
  };

 private:
  /// Factory method to create a surface grid lookup for a given representative
  /// surface, tolerance, and axes. This will internally create the appropriate
  /// lookup class based on the axes and concrete @ref Grid.
  /// @param representative The surface which is used as representative
  /// @param tolerance The tolerance used for intersection checks
  /// @param axes The axes used for the grid
  /// @param bValues Optional vector of axis directions for binning
  /// @return A unique pointer to the surface grid lookup
  static std::unique_ptr<ISurfaceGridLookup> makeSurfaceGridLookup(
      std::shared_ptr<RegularSurface> representative, double tolerance,
      std::tuple<const IAxis&, const IAxis&> axes);

  // This is temporary until Gen1 is removed
  friend class SurfaceArrayCreator;

 public:
  /// @deprecated Use @ref makeSurfaceGridLookup instead. The construction of
  ///             this class is now handled internally and becomes an
  ///             implementation detail.
  template <class Axis1, class Axis2>
  struct [[deprecated("Use makeSurfaceGridLookup instead")]] SurfaceGridLookup
      : ISurfaceGridLookup {
    /// Construct a surface grid lookup
    /// @param representative The surface which is used as representative
    /// @param tolerance The tolerance used for intersection checks
    /// @param axes The axes used for the grid
    /// @param bValues Optional vector of axis directions for binning
    SurfaceGridLookup(std::shared_ptr<RegularSurface> representative,
                      double tolerance, std::tuple<Axis1, Axis2> axes,
                      const std::vector<AxisDirection>& bValues = {})

        : m_impl(makeSurfaceGridLookup(std::move(representative), tolerance,
                                       axes)) {
      static_cast<void>(bValues);
    }

    /// @brief Fill provided surfaces into the contained @c Grid.
    /// @param gctx The current geometry context object, e.g. alignment
    /// @param surfaces Input surface pointers
    void fill(const GeometryContext& gctx,
              const SurfaceVector& surfaces) override {
      m_impl->fill(gctx, surfaces);
    }

    /// @brief Performs lookup at @c pos and returns bin content as const
    /// reference
    /// @param position Lookup position
    /// @param direction Lookup direction
    /// @return @c SurfaceVector at given bin
    const SurfaceVector& lookup(const Vector3& position,
                                const Vector3& direction) const override {
      return m_impl->lookup(position, direction);
    }

    /// @brief Performs lookup at global bin and returns bin content as
    /// reference
    /// @param bin Global lookup bin
    /// @return @c SurfaceVector at given bin
    SurfaceVector& lookup(std::size_t bin) override {
      return m_impl->lookup(bin);
    }

    /// @brief Performs lookup at global bin and returns bin content as const
    /// reference
    /// @param bin Global lookup bin
    /// @return @c SurfaceVector at given bin
    const SurfaceVector& lookup(std::size_t bin) const override {
      return m_impl->lookup(bin);
    }

    /// @brief Performs a lookup at @c pos, but returns neighbors as well
    ///
    /// @param position Lookup position
    /// @param direction Lookup direction
    /// @return @c SurfaceVector at given bin. Copy of all bins selected
    const SurfaceVector& neighbors(const Vector3& position,
                                   const Vector3& direction) const override {
      return m_impl->neighbors(position, direction);
    }

    /// @brief Returns the total size of the grid (including under/overflow
    /// bins)
    /// @return Size of the grid data structure
    std::size_t size() const override { return m_impl->size(); }

    /// @brief Gets the center position of bin @c bin in global coordinates
    /// @param bin the global bin index
    /// @return The bin center
    Vector3 getBinCenter(std::size_t bin) const override {
      return m_impl->getBinCenter(bin);
    }

    /// @brief Returns copies of the axes used in the grid as @c AnyAxis
    /// @return The axes
    /// @note This returns copies. Use for introspection and querying.
    std::vector<const IAxis*> getAxes() const override {
      return m_impl->getAxes();
    }

    /// @brief Get a view of the grid for inspection
    /// @return Optional grid view containing surface vectors
    std::optional<AnyGridConstView<SurfaceVector>> getGridView()
        const override {
      return m_impl->getGridView();
    }

    /// @brief Get the representative surface used for this lookup
    /// @return Surface pointer
    const Surface* surfaceRepresentation() const override {
      return m_impl->surfaceRepresentation();
    }

    /// @brief Checks if global bin is valid
    /// @param bin the global bin index
    /// @return bool if the bin is valid
    /// @note Valid means that the index points to a bin which is not a under
    ///       or overflow bin or out of range in any axis.
    bool isValidBin(std::size_t bin) const override {
      return m_impl->isValidBin(bin);
    }

    /// @brief The binning values described by this surface grid lookup
    /// They are in order of the axes (optional) and empty for eingle lookups
    /// @return Vector of axis directions for binning
    std::vector<AxisDirection> binningValues() const override {
      return m_impl->binningValues();
    }

    /// Pure virtual destructor
    ~SurfaceGridLookup() override = default;

   private:
    std::unique_ptr<ISurfaceGridLookup> m_impl;
  };

  /// @brief Lookup implementation which wraps one element and always returns
  ///        this element when lookup is called
  /// @deprecated Construct the @ref SurfaceArray directly with a single surface
  struct [[deprecated(
      "Construct the SurfaceArray directly with a single "
      "surface")]] SingleElementLookup : ISurfaceGridLookup {
    /// @brief Default constructor.
    /// @param element the one and only element.
    explicit SingleElementLookup(SurfaceVector::value_type element)
        : m_element({element}) {}

    /// @brief Default constructor.
    /// @param elements the surfaces that are provided through a single lookup
    explicit SingleElementLookup(const SurfaceVector& elements)
        : m_element(elements) {}

    /// @brief Lookup, always returns @c element
    /// @return reference to vector containing only @c element
    const SurfaceVector& lookup(const Vector3& /*position*/,
                                const Vector3& /*direction*/) const override {
      return m_element;
    }

    /// @brief Lookup, always returns @c element
    /// @return reference to vector containing only @c element
    SurfaceVector& lookup(std::size_t /*bin*/) override { return m_element; }

    /// @brief Lookup, always returns @c element
    /// @return reference to vector containing only @c element
    const SurfaceVector& lookup(std::size_t /*bin*/) const override {
      return m_element;
    }

    /// @brief Lookup, always returns @c element
    /// @return reference to vector containing only @c element
    const SurfaceVector& neighbors(
        const Vector3& /*position*/,
        const Vector3& /*direction*/) const override {
      return m_element;
    }

    /// @brief returns 1
    /// @return 1
    std::size_t size() const override { return 1; }

    /// @brief Gets the bin center, but always returns (0, 0, 0)
    /// @return (0, 0, 0)
    Vector3 getBinCenter(std::size_t /*bin*/) const override {
      return Vector3(0, 0, 0);
    }

    /// @brief Returns an empty vector of @c AnyAxis
    /// @return empty vector
    std::vector<const IAxis*> getAxes() const override { return {}; }

    std::optional<AnyGridConstView<SurfaceVector>> getGridView()
        const override {
      return std::nullopt;
    }

    const Surface* surfaceRepresentation() const override { return nullptr; }

    /// @brief Comply with concept and provide fill method
    /// @note Does nothing
    void fill(const GeometryContext& /*gctx*/,
              const SurfaceVector& /*surfaces*/) override {}

    /// @brief Returns if the bin is valid (it is)
    /// @return always true
    bool isValidBin(std::size_t /*bin*/) const override { return true; }

   private:
    SurfaceVector m_element;
  };

  /// @brief Default constructor which takes a @c SurfaceLookup and a vector of
  /// surfaces
  /// @param gridLookup The grid storage. @c SurfaceArray does not fill it on
  /// its own
  /// @param surfaces The input vector of surfaces. This is only for
  /// bookkeeping, so we can ask
  /// @param transform Optional additional transform for this SurfaceArray
  [[deprecated("Use the constructor with axes instead")]]
  SurfaceArray(std::unique_ptr<ISurfaceGridLookup> gridLookup,
               std::vector<std::shared_ptr<const Surface>> surfaces,
               const Transform3& transform = Transform3::Identity());

  /// @brief Constructor with a single surface
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

  /// @brief Get all surfaces in bin given by position @p pos.
  /// @param position the lookup position
  /// @param direction the lookup direction
  /// @return const reference to @c SurfaceVector contained in bin at that
  /// position
  const SurfaceVector& at(const Vector3& position,
                          const Vector3& direction) const {
    return p_gridLookup->lookup(position, direction);
  }

  /// @brief Get all surfaces in bin given by global bin index @p bin.
  /// @param bin the global bin index
  /// @return reference to @c SurfaceVector contained in bin
  SurfaceVector& at(std::size_t bin) { return p_gridLookup->lookup(bin); }

  /// @brief Get all surfaces in bin given by global bin index.
  /// @param bin the global bin index
  /// @return const reference to @c SurfaceVector contained in bin
  const SurfaceVector& at(std::size_t bin) const {
    return p_gridLookup->lookup(bin);
  }

  /// @brief Get all surfaces in bin at @p pos and its neighbors
  /// @param position The position to lookup
  /// @param direction The direction to lookup
  /// @return Merged @c SurfaceVector of neighbors and nominal
  /// @note The @c SurfaceVector will be combined. For technical reasons, the
  ///       different bin content vectors have to be copied, so the resulting
  ///       vector contains copies.
  const SurfaceVector& neighbors(const Vector3& position,
                                 const Vector3& direction) const {
    return p_gridLookup->neighbors(position, direction);
  }

  /// @brief Get the size of the underlying grid structure including
  /// under/overflow bins
  /// @return the size
  std::size_t size() const { return p_gridLookup->size(); }

  /// @brief Get the center of the bin identified by global bin index @p bin
  /// @param bin the global bin index
  /// @return Center position of the bin in global coordinates
  Vector3 getBinCenter(std::size_t bin) const {
    return p_gridLookup->getBinCenter(bin);
  }

  /// @brief Get all surfaces attached to this @c SurfaceArray
  /// @return Reference to @c SurfaceVector containing all surfaces
  /// @note This does not reflect the actual state of the grid. It only
  ///       returns what was given in the constructor, without any checks
  ///       if that is actually what's in the grid.
  const SurfaceVector& surfaces() const { return m_surfacesRawPointers; }

  /// @brief Get vector of axes spanning the grid as @c AnyAxis
  /// @return vector of @c AnyAxis
  /// @note The axes in the vector are copies. Only use for introspection and
  ///       querying.
  std::vector<const IAxis*> getAxes() const { return p_gridLookup->getAxes(); }

  /// @brief Checks if global bin is valid
  /// @param bin the global bin index
  /// @return bool if the bin is valid
  /// @note Valid means that the index points to a bin which is not a under
  ///       or overflow bin or out of range in any axis.
  bool isValidBin(std::size_t bin) const {
    return p_gridLookup->isValidBin(bin);
  }

  /// Get the transform of this surface array.
  /// @return Reference to the transformation matrix
  /// @deprecated This is an implementation detail and will be removed soon
  [[deprecated("This is an implementation detail and will be removed soon")]]
  const Transform3& transform() const {
    return m_transform;
  }

  /// @brief The binning values described by this surface grid lookup
  /// They are in order of the axes
  /// @return Vector of axis directions for binning
  std::vector<AxisDirection> binningValues() const {
    return p_gridLookup->binningValues();
  };

  /// @brief String representation of this @c SurfaceArray
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param sl Output stream to write to
  /// @return the output stream given as @p sl
  std::ostream& toStream(const GeometryContext& gctx, std::ostream& sl) const;

  /// Return the lookup object
  /// @return Reference to the surface grid lookup interface
  [[deprecated(
      "Grid lookup is an implementation detail and will be removed soon")]]
  const ISurfaceGridLookup& gridLookup() const {
    return *p_gridLookup;
  }

  /// Get the representative surface used for this surface array
  /// @return Surface pointer
  const Surface* surfaceRepresentation() const {
    return p_gridLookup->surfaceRepresentation();
  }

  /// Get a view of the grid for inspection
  /// @return Optional grid view containing surface vectors
  std::optional<AnyGridConstView<SurfaceVector>> getGridView() const {
    return p_gridLookup->getGridView();
  }

 private:
  /// Check consistency between provided surfaces and grid contents.
  ///
  /// Iterates over all local grid bins, collects every surface pointer seen in
  /// the bins, and compares that set against the surfaces provided to this
  /// array. Throws if the sets differ (e.g. a provided surface is not present
  /// in the grid).
  ///
  /// @param grid The grid to check
  void checkGrid(AnyGridConstView<SurfaceVector> grid);

  std::unique_ptr<ISurfaceGridLookup> p_gridLookup;
  // this vector makes sure we have shared ownership over the surfaces
  std::vector<std::shared_ptr<const Surface>> m_surfaces;
  // this vector is returned, so that (expensive) copying of the shared_ptr
  // vector does not happen by default
  SurfaceVector m_surfacesRawPointers;
  // this is only used to keep info on transform applied
  // by l2g and g2l
  Transform3 m_transform;
};

}  // namespace Acts
