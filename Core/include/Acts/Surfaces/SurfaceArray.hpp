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
#include "Acts/Utilities/AnyGridView.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/IAxis.hpp"

#include <iostream>
#include <type_traits>
#include <vector>

namespace Acts {

using SurfaceVector = std::vector<const Surface*>;

/// @brief Provides Surface binning in N dimensions
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

    /// @brief Attempts to fix sub-optimal binning by filling closest
    ///        Surfaces into empty bin
    ///
    /// @param gctx The current geometry context object, e.g. alignment

    /// @param surfaces The surface pointers to fill
    /// @return number of bins that were filled
    virtual std::size_t completeBinning(const GeometryContext& gctx,
                                        const SurfaceVector& surfaces) = 0;

    /// @brief Performs lookup at @c pos and returns bin content as const
    /// reference
    /// @param position Lookup position
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

    virtual std::optional<AnyGridConstView<SurfaceVector>> getGridView()
        const = 0;

    virtual Surface::SurfaceType surfaceType() const = 0;

    /// @brief Get the number of dimensions of the grid.
    /// @return number of dimensions
    virtual std::size_t dimensions() const = 0;

    /// @brief Checks if global bin is valid
    /// @param bin the global bin index
    /// @return bool if the bin is valid
    /// @note Valid means that the index points to a bin which is not a under
    ///       or overflow bin or out of range in any axis.
    virtual bool isValidBin(std::size_t bin) const = 0;

    /// @brief The binning values described by this surface grid lookup
    /// They are in order of the axes (optional) and empty for eingle lookups
    virtual std::vector<AxisDirection> binningValues() const { return {}; };

    /// Pure virtual destructor
    virtual ~ISurfaceGridLookup() = 0;
  };

  /// @brief Lookup helper which encapsulates a @c Grid
  /// @tparam Axes The axes used for the grid
  template <class... Axes>
  struct SurfaceGridLookup : ISurfaceGridLookup {
    static constexpr std::size_t DIM = sizeof...(Axes);

   public:
    /// @brief Specifies the local coordinate type.
    /// This resolves to @c ActsVector<DIM> for DIM > 1, else @c
    /// std::array<double, 1>
    using point_t =
        std::conditional_t<DIM == 1, std::array<double, 1>, ActsVector<DIM>>;
    using Grid_t = Grid<SurfaceVector, Axes...>;

    SurfaceGridLookup(std::shared_ptr<RegularSurface> representative,
                      double tolerance, std::tuple<Axes...> axes,
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
      // surface to bin matching
      for (const Surface* srf : surfaces) {
        Vector3 pos = srf->referencePosition(gctx, AxisDirection::AxisR);
        Vector3 normal = m_representative->normal(gctx, pos);
        auto intersection =
            m_representative->intersect(gctx, pos, normal).closest();
        if (!intersection.isValid() ||
            std::abs(intersection.pathLength()) > m_tolerance) {
          continue;
        }
        Vector2 lposition =
            m_representative
                ->globalToLocal(gctx, intersection.position(), normal)
                .value();
        m_grid.atPosition(lposition).push_back(srf);
      }

      // bin to surface matching
      for (std::size_t i = 0; i < m_grid.size(); i++) {
        auto j = m_grid.localBinsFromGlobalBin(i);
        Vector2 local = Eigen::Map<Vector2>(m_grid.binCenter(j).data());
        Vector3 normal = m_representative->normal(gctx, local);
        Vector3 global = m_representative->localToGlobal(gctx, local, normal);

        for (const Surface* srf : surfaces) {
          auto intersection = srf->intersect(gctx, global, normal).closest();
          if (!intersection.isValid() ||
              std::abs(intersection.pathLength()) > m_tolerance) {
            continue;
          }
          m_grid.at(i).push_back(srf);
        }
      }

      populateNeighborCache();

      // deduplicate
      for (std::size_t i = 0; i < m_grid.size(); i++) {
        auto& binContent = m_grid.at(i);
        std::ranges::sort(binContent);
        auto last = std::ranges::unique(binContent);
        binContent.erase(last.begin(), last.end());
        binContent.shrink_to_fit();
      }
    }

    /// @brief Attempts to fix sub-optimal binning by filling closest
    ///        Surfaces into empty bins
    /// @note This does not always do what you want.
    ///
    /// @param gctx The current geometry context object, e.g. alignment
    /// @param surfaces The surface pointers to fill
    /// @return number of bins that were filled
    std::size_t completeBinning(const GeometryContext& gctx,
                                const SurfaceVector& surfaces) override {
      std::size_t binCompleted = 0;
      std::size_t nBins = size();
      double minPath = 0;
      double curPath = 0;
      const Surface* minSrf = nullptr;

      for (std::size_t b = 0; b < nBins; ++b) {
        if (!isValidBin(b)) {
          continue;
        }
        std::vector<const Surface*>& binContent = lookup(b);
        // only complete if we have an empty bin
        if (!binContent.empty()) {
          continue;
        }

        Vector3 binCtr = getBinCenter(b);
        minPath = std::numeric_limits<double>::max();
        for (const auto& srf : surfaces) {
          curPath =
              (binCtr - srf->referencePosition(gctx, AxisDirection::AxisR))
                  .norm();

          if (curPath < minPath) {
            minPath = curPath;
            minSrf = srf;
          }
        }

        binContent.push_back(minSrf);
        ++binCompleted;
      }

      // recreate neighborcache
      populateNeighborCache();
      return binCompleted;
    }

    const SurfaceVector& lookup(const Vector3& position,
                                const Vector3& direction) const override {
      GeometryContext gctx;

      auto intersection =
          m_representative->intersect(gctx, position, direction).closest();
      Vector2 lposition =
          m_representative
              ->globalToLocal(gctx, intersection.position(), direction)
              .value();
      return m_grid.atPosition(lposition);
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
    /// @return @c SurfaceVector at given bin. Copy of all bins selected
    const SurfaceVector& neighbors(const Vector3& position,
                                   const Vector3& direction) const override {
      GeometryContext gctx;

      auto intersection =
          m_representative->intersect(gctx, position, direction).closest();
      Vector2 lposition =
          m_representative
              ->globalToLocal(gctx, intersection.position(), direction)
              .value();
      return m_neighborMap.at(m_grid.globalBinFromPosition(lposition));
    }

    /// @brief Returns the total size of the grid (including under/overflow
    /// bins)
    /// @return Size of the grid data structure
    std::size_t size() const override { return m_grid.size(); }

    /// @brief The binning values described by this surface grid lookup
    /// They are in order of the axes
    std::vector<AxisDirection> binningValues() const override {
      return m_binValues;
    }

    /// @brief Gets the center position of bin @c bin in global coordinates
    /// @param bin the global bin index
    /// @return The bin center
    Vector3 getBinCenter(std::size_t bin) const override {
      GeometryContext gctx;
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

    Surface::SurfaceType surfaceType() const override {
      return m_representative->type();
    }

    /// @brief Get the number of dimensions of the grid.
    /// @return number of dimensions
    std::size_t dimensions() const override { return DIM; }

    /// @brief Checks if global bin is valid
    /// @param bin the global bin index
    /// @return bool if the bin is valid
    /// @note Valid means that the index points to a bin which is not a under
    ///       or overflow bin or out of range in any axis.
    bool isValidBin(std::size_t bin) const override {
      std::array<std::size_t, DIM> indices = m_grid.localBinsFromGlobalBin(bin);
      std::array<std::size_t, DIM> nBins = m_grid.numLocalBins();
      for (std::size_t i = 0; i < indices.size(); ++i) {
        std::size_t idx = indices.at(i);
        if (idx <= 0 || idx >= nBins.at(i) + 1) {
          return false;
        }
      }

      return true;
    }

   private:
    void populateNeighborCache() {
      // calculate neighbors for every bin and store in map
      for (std::size_t i = 0; i < m_grid.size(); i++) {
        if (!isValidBin(i)) {
          continue;
        }
        typename Grid_t::index_t loc = m_grid.localBinsFromGlobalBin(i);
        auto neighborIdxs = m_grid.neighborHoodIndices(loc, 1u);
        std::vector<const Surface*>& neighbors = m_neighborMap.at(i);
        neighbors.clear();

        for (const auto idx : neighborIdxs) {
          const std::vector<const Surface*>& binContent = m_grid.at(idx);
          std::copy(binContent.begin(), binContent.end(),
                    std::back_inserter(neighbors));
        }
      }
    }

    /// Internal method.
    /// This is here, because apparently Eigen doesn't like Vector1.
    /// So SurfaceGridLookup internally uses std::array<double, 1> instead
    /// of Vector1 (see the point_t typedef). This needs to be switched here,
    /// so as not to
    /// attempt an initialization of Vector1 that Eigen will complain about.
    /// The SFINAE is hidden in this private method so the public
    /// interface stays the same, since we don't care what happens
    /// here on the callers end
    /// This is the version for DIM>1
    Vector3 getBinCenterImpl(const GeometryContext& gctx, std::size_t bin) const
      requires(DIM != 1)
    {
      return m_representative->localToGlobal(
          gctx,
          ActsVector<DIM>(
              m_grid.binCenter(m_grid.localBinsFromGlobalBin(bin)).data()));
    }

    /// Internal method, see above.
    /// This is the version for DIM==1
    Vector3 getBinCenterImpl(const GeometryContext& gctx, std::size_t bin) const
      requires(DIM == 1)
    {
      return m_representative->localToGlobal(
          gctx, m_grid.binCenter(m_grid.localBinsFromGlobalBin(bin)));
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

    Surface::SurfaceType surfaceType() const override {
      return Surface::SurfaceType::Other;
    }

    /// @brief Get the number of dimensions
    /// @return always 0
    std::size_t dimensions() const override { return 0; }

    /// @brief Comply with concept and provide fill method
    /// @note Does nothing
    void fill(const GeometryContext& /*gctx*/,
              const SurfaceVector& /*surfaces*/) override {}

    /// @brief Comply with concept and provide completeBinning method
    /// @note Does nothing
    std::size_t completeBinning(const GeometryContext& /*gctx*/,
                                const SurfaceVector& /*surfaces*/) override {
      return 0;
    }

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
  /// @param position The position to lookup as nominal
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

  const Transform3& transform() const { return m_transform; }

  /// @brief The binning values described by this surface grid lookup
  /// They are in order of the axes
  std::vector<AxisDirection> binningValues() const {
    return p_gridLookup->binningValues();
  };

  /// @brief String representation of this @c SurfaceArray
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param sl Output stream to write to
  /// @return the output stream given as @p sl
  std::ostream& toStream(const GeometryContext& gctx, std::ostream& sl) const;

  /// Return the lookup object
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
