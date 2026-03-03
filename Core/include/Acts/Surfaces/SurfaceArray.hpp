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

  /// @brief Lookup helper which encapsulates a @c Grid
  /// @tparam Axes The axes used for the grid
  template <class Axis1, class Axis2>
  struct SurfaceGridLookup : ISurfaceGridLookup {
    /// Grid type storing surface vectors with two axes
    using Grid_t = Grid<SurfaceVector, Axis1, Axis2>;

    /// Construct a surface grid lookup
    /// @param representative The surface which is used as representative
    /// @param tolerance The tolerance used for intersection checks
    /// @param axes The axes used for the grid
    /// @param bValues Optional vector of axis directions for binning
    SurfaceGridLookup(std::shared_ptr<RegularSurface> representative,
                      double tolerance, std::tuple<Axis1, Axis2> axes,
                      std::vector<AxisDirection> bValues = {})
        : m_representative(std::move(representative)),
          m_tolerance(tolerance),
          m_grid(std::move(axes)),
          m_binValues(std::move(bValues)) {
      m_neighborMap.resize(m_grid.size());
    }

    /// @brief Fill provided surfaces into the contained @c Grid.
    ///
    /// This is done by iterating, accessing the referencePosition, lookup
    /// and append.
    /// Also populates the neighbor map by combining the filled bins of
    /// all bins around a given one.
    ///
    /// @param gctx The current geometry context object, e.g. alignment
    /// @param surfaces Input surface pointers
    void fill(const GeometryContext& gctx,
              const SurfaceVector& surfaces) override {
      for (const Surface* surface : surfaces) {
        const std::size_t globalBin = fillSurfaceToBinMapping(gctx, *surface);
        if (globalBin == 0) {
          continue;
        }

        fillBinToSurfaceMapping(gctx, *surface, globalBin);
      }

      populateNeighborCache();
    }

    const SurfaceVector& lookup(const Vector3& position,
                                const Vector3& direction) const override {
      return m_grid.at(findGlobalBin(position, direction,
                                     std::numeric_limits<double>::infinity()));
    }

    /// @brief Performs lookup at global bin and returns bin content as
    /// reference
    /// @param bin Global lookup bin
    /// @return @c SurfaceVector at given bin
    SurfaceVector& lookup(std::size_t bin) override { return m_grid.at(bin); }

    /// @brief Performs lookup at global bin and returns bin content as const
    /// reference
    /// @param bin Global lookup bin
    /// @return @c SurfaceVector at given bin
    const SurfaceVector& lookup(std::size_t bin) const override {
      return m_grid.at(bin);
    }

    /// @brief Performs a lookup at @c pos, but returns neighbors as well
    ///
    /// @param position Lookup position
    /// @param direction Lookup direction
    /// @return @c SurfaceVector at given bin. Copy of all bins selected
    const SurfaceVector& neighbors(const Vector3& position,
                                   const Vector3& direction) const override {
      return m_neighborMap.at(findGlobalBin(
          position, direction, std::numeric_limits<double>::infinity()));
    }

    /// @brief Returns the total size of the grid (including under/overflow
    /// bins)
    /// @return Size of the grid data structure
    std::size_t size() const override { return m_grid.size(); }

    /// @brief The binning values described by this surface grid lookup
    /// They are in order of the axes
    /// @return Vector of axis directions for binning
    std::vector<AxisDirection> binningValues() const override {
      return m_binValues;
    }

    /// @brief Gets the center position of bin @c bin in global coordinates
    /// @param bin the global bin index
    /// @return The bin center
    Vector3 getBinCenter(std::size_t bin) const override {
      auto gctx = GeometryContext::dangerouslyDefaultConstruct();
      return getBinCenterImpl(gctx, bin);
    }

    /// @brief Returns copies of the axes used in the grid as @c AnyAxis
    /// @return The axes
    /// @note This returns copies. Use for introspection and querying.
    std::vector<const IAxis*> getAxes() const override {
      auto arr = m_grid.axes();
      return std::vector<const IAxis*>(arr.begin(), arr.end());
    }

    std::optional<AnyGridConstView<SurfaceVector>> getGridView()
        const override {
      return AnyGridConstView<SurfaceVector>{m_grid};
    }

    const Surface* surfaceRepresentation() const override {
      return m_representative.get();
    }

    /// @brief Checks if global bin is valid
    /// @param bin the global bin index
    /// @return bool if the bin is valid
    /// @note Valid means that the index points to a bin which is not a under
    ///       or overflow bin or out of range in any axis.
    bool isValidBin(std::size_t bin) const override {
      std::array<std::size_t, 2> indices = m_grid.localBinsFromGlobalBin(bin);
      std::array<std::size_t, 2> nBins = m_grid.numLocalBins();
      for (std::size_t i = 0; i < indices.size(); ++i) {
        std::size_t idx = indices.at(i);
        if (idx <= 0 || idx >= nBins.at(i) + 1) {
          return false;
        }
      }
      return true;
    }

   private:
    /// map surface center to grid
    std::size_t fillSurfaceToBinMapping(const GeometryContext& gctx,
                                        const Surface& surface) {
      const Vector3 pos = surface.referencePosition(gctx, AxisDirection::AxisR);
      const Vector3 normal = m_representative->normal(gctx, pos);
      const std::size_t globalBin = findGlobalBin(pos, normal, m_tolerance);
      if (globalBin != 0) {
        m_grid.at(globalBin).push_back(&surface);
      }
      return globalBin;
    };

    /// flood fill neighboring bins given a starting bin
    void fillBinToSurfaceMapping(const GeometryContext& gctx,
                                 const Surface& surface, std::size_t startBin) {
      const std::array<std::size_t, 2> startIndices =
          m_grid.localBinsFromGlobalBin(startBin);
      const auto startNeighborIndices =
          m_grid.neighborHoodIndices(startIndices, 1u);

      std::set<std::size_t> visited({startBin});
      std::vector<std::size_t> queue(startNeighborIndices.begin(),
                                     startNeighborIndices.end());

      while (!queue.empty()) {
        const std::size_t current = queue.back();
        queue.pop_back();
        if (visited.contains(current)) {
          continue;
        }

        const std::array<std::size_t, 2> currentIndices =
            m_grid.localBinsFromGlobalBin(current);
        visited.insert(current);

        const std::array<double, 2> gridLocal =
            m_grid.binCenter(currentIndices);
        const Vector2 surfaceLocal = gridToSurfaceLocal(gridLocal);
        const Vector3 normal = m_representative->normal(gctx, surfaceLocal);
        const Vector3 global =
            m_representative->localToGlobal(gctx, surfaceLocal, normal);

        const Intersection3D intersection =
            surface.intersect(gctx, global, normal, BoundaryTolerance::None())
                .closest();
        if (!intersection.isValid() ||
            std::abs(intersection.pathLength()) > m_tolerance) {
          continue;
        }
        m_grid.at(current).push_back(&surface);

        const auto neighborIndices =
            m_grid.neighborHoodIndices(currentIndices, 1u);
        queue.insert(queue.end(), neighborIndices.begin(),
                     neighborIndices.end());
      }
    };

    void populateNeighborCache() {
      // calculate neighbors for every bin and store in map
      for (std::size_t i = 0; i < m_grid.size(); i++) {
        if (!isValidBin(i)) {
          continue;
        }
        const std::array<std::size_t, 2> indices =
            m_grid.localBinsFromGlobalBin(i);
        std::vector<const Surface*>& neighbors = m_neighborMap.at(i);
        neighbors.clear();

        for (std::size_t idx : m_grid.neighborHoodIndices(indices, 1u)) {
          const std::vector<const Surface*>& binContent = m_grid.at(idx);
          std::copy(binContent.begin(), binContent.end(),
                    std::back_inserter(neighbors));
        }

        std::ranges::sort(neighbors);
        auto last = std::ranges::unique(neighbors);
        neighbors.erase(last.begin(), last.end());
        neighbors.shrink_to_fit();
      }
    }

    Vector3 getBinCenterImpl(const GeometryContext& gctx,
                             std::size_t bin) const {
      const std::array<double, 2> gridLocal =
          m_grid.binCenter(m_grid.localBinsFromGlobalBin(bin));
      const Vector2 surfaceLocal = gridToSurfaceLocal(gridLocal);
      return m_representative->localToGlobal(gctx, surfaceLocal);
    }

    const CylinderBounds* getCylinderBounds() const {
      return dynamic_cast<const CylinderBounds*>(&m_representative->bounds());
    }

    Vector2 gridToSurfaceLocal(std::array<double, 2> gridLocal) const {
      Vector2 surfaceLocal = Eigen::Map<Vector2>(gridLocal.data());
      if (const CylinderBounds* bounds = getCylinderBounds();
          bounds != nullptr) {
        surfaceLocal[0] *= bounds->get(CylinderBounds::eR);
      }
      return surfaceLocal;
    }
    std::array<double, 2> surfaceToGridLocal(Vector2 local) const {
      std::array<double, 2> gridLocal = {local[0], local[1]};
      if (const CylinderBounds* bounds = getCylinderBounds();
          bounds != nullptr) {
        gridLocal[0] /= bounds->get(CylinderBounds::eR);
      }
      return gridLocal;
    }

    std::size_t findGlobalBin(const Vector3& position, const Vector3& direction,
                              double tolerance) const {
      auto gctx = GeometryContext::dangerouslyDefaultConstruct();

      const Intersection3D intersection =
          m_representative
              ->intersect(gctx, position, direction,
                          BoundaryTolerance::Infinite())
              .closest();
      if (!intersection.isValid() ||
          std::abs(intersection.pathLength()) > tolerance) {
        return 0;  // overflow bin
      }
      const Vector2 surfaceLocal =
          m_representative
              ->globalToLocal(gctx, intersection.position(), direction)
              .value();
      const std::array<double, 2> gridLocal = surfaceToGridLocal(surfaceLocal);
      return m_grid.globalBinFromPosition(gridLocal);
    }

    std::shared_ptr<RegularSurface> m_representative;
    double m_tolerance{};
    Grid_t m_grid;
    std::vector<AxisDirection> m_binValues;
    std::vector<SurfaceVector> m_neighborMap;
  };

  /// @brief Lookup implementation which wraps one element and always returns
  ///        this element when lookup is called
  struct SingleElementLookup : ISurfaceGridLookup {
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
  explicit SurfaceArray(std::unique_ptr<ISurfaceGridLookup> gridLookup,
                        std::vector<std::shared_ptr<const Surface>> surfaces,
                        const Transform3& transform = Transform3::Identity());

  /// @brief Constructor with a single surface
  /// @param srf The one and only surface
  explicit SurfaceArray(std::shared_ptr<const Surface> srf);

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
  const Transform3& transform() const { return m_transform; }

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
  const ISurfaceGridLookup& gridLookup() const { return *p_gridLookup; }

 private:
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
