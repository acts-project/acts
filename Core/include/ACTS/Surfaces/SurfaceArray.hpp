// This file is part of the ACTS project.
//
// Copyright (C) 2017-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_SURFACES_SURFACEARRAY_H
#define ACTS_SURFACES_SURFACEARRAY_H

#include <iostream>
#include <type_traits>
#include <vector>
#include "ACTS/Surfaces/concept/AnySurfaceGridLookup.hpp"
#include "ACTS/Utilities/concept/AnyGrid.hpp"
#include "ACTS/Utilities/detail/Axis.hpp"
#include "ACTS/Utilities/detail/Grid.hpp"
//#include <boost/any.hpp>
#include <boost/type_erasure/any.hpp>
#include <boost/type_erasure/any_cast.hpp>
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/BinningType.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

namespace bte = boost::type_erasure;

using SurfaceVector = std::vector<const Surface*>;
template <class... Axes>
using SurfaceGrid = detail::Grid<SurfaceVector, Axes...>;

/// @brief Provides Surface binning in N dimensions
///
/// Uses @c Grid under the hood to implement the storage and lookup
/// Contains a type-erased lookup struct which talks to the @c Grid
/// and performs utility actions. This struct needs to be initialised
/// externally and passed to @c SurfaceArray on construction.
class SurfaceArray
{
  // typedef to the Grid type used
  template <class Point, size_t DIM>
  using AnyGrid_t = concept::AnyNDimGrid<SurfaceVector, Point, DIM>;

  friend std::ostream&
  operator<<(std::ostream& sl, const SurfaceArray& sa)
  {
    return sa.dump(sl);
  }

public:
  using AnySurfaceGridLookup_t = concept::AnySurfaceGridLookup<SurfaceVector>;

  /// @brief Lookup helper which encapsulates a @c Grid
  /// @tparam Axes The axes used for the grid
  template <class... Axes>
  struct SurfaceGridLookup
  {
    static constexpr size_t DIM = sizeof...(Axes);

  public:
    /// @brief Specifies the local coordinate type.
    /// This resolves to @c ActsVector<DIM> for DIM > 1, else @c
    /// std::array<double, 1>
    using point_t
        = std::conditional_t<DIM == 1, std::array<double, 1>, ActsVectorD<DIM>>;
    using Grid_t = detail::Grid<SurfaceVector, Axes...>;

    /// @brief Default constructor
    ///
    /// @param globalToLocal Callable that converts from global to local
    /// @param localToGlobal Callable that converts from local to global
    /// @param grid The grid data structur. Will be type-erased at this point
    /// @note Signature of localToGlobal and globalToLocal depends on @c DIM.
    ///       If DIM > 1, local coords are @c ActsVectorD<DIM> else
    ///       @c std::array<double, 1>.
    SurfaceGridLookup(std::function<point_t(const Vector3D&)> globalToLocal,
                      std::function<Vector3D(const point_t&)> localToGlobal,
                      std::tuple<Axes...>                     axes)
      : m_globalToLocal(std::move(globalToLocal))
      , m_localToGlobal(std::move(localToGlobal))
      , m_grid(std::move(axes))
    {
      m_neighborMap.resize(m_grid.size());
    }

    /// @brief Fill provided surfaces into the contained @c Grid.
    ///
    /// This is done by iterating, accessing the binningPosition, lookup
    /// and append.
    /// Also populates the neighbor map by combining the filled bins of
    /// all bins around a given one
    ///
    /// @param surfaces Input surface pointers
    void
    fill(const SurfaceVector& surfaces)
    {
      for (const auto& srf : surfaces) {
        Vector3D pos = srf->binningPosition(binR);
        lookup(pos).push_back(srf);
      }

      populateNeighborCache();
    }

    /// @brief Attempts to fix sub-optimal binning by filling closest
    ///        Surfaces into empty bins
    ///
    /// @param surfaces The surface pointers to fill
    /// @return number of bins that were filled
    size_t
    completeBinning(const SurfaceVector& surfaces)
    {
      size_t         binCompleted = 0;
      size_t         nBins        = size();
      double         minPath, curPath;
      const Surface* minSrf;

      for (size_t b = 0; b < nBins; ++b) {
        if (!isValidBin(b)) continue;
        std::vector<const Surface*>& binContent = lookup(b);
        // only complete if we have an empty bin
        if (binContent.size() > 0) continue;

        Vector3D binCtr = getBinCenter(b);
        minPath         = std::numeric_limits<double>::max();
        for (const auto& srf : surfaces) {
          curPath = (binCtr - srf->binningPosition(binR)).mag();

          if (curPath < minPath) {
            minPath = curPath;
            minSrf  = srf;
          }
        }

        binContent.push_back(minSrf);
        ++binCompleted;
      }

      // recreate neighborcache
      populateNeighborCache();
      return binCompleted;
    }

    /// @brief Performs lookup at @c pos and returns bin content as reference
    /// @param pos Lookup position
    /// @return @c SurfaceVector at given bin
    SurfaceVector&
    lookup(const Vector3D& pos)
    {
      return m_grid.at(m_globalToLocal(pos));
    }

    /// @brief Performs lookup at @c pos and returns bin content as const
    /// reference
    /// @param pos Lookup position
    /// @return @c SurfaceVector at given bin
    const SurfaceVector&
    lookup(const Vector3D& pos) const
    {
      return m_grid.at(m_globalToLocal(pos));
    }

    /// @brief Performs lookup at global bin and returns bin content as
    /// reference
    /// @param bin Global lookup bin
    /// @return @c SurfaceVector at given bin
    SurfaceVector&
    lookup(size_t bin)
    {
      return m_grid.at(bin);
    }

    /// @brief Performs lookup at global bin and returns bin content as const
    /// reference
    /// @param bin Global lookup bin
    /// @return @c SurfaceVector at given bin
    const SurfaceVector&
    lookup(size_t bin) const
    {
      return m_grid.at(bin);
    }

    /// @brief Performs a lookup at @c pos, but returns neighbors as well
    ///
    /// @param pos Lookup position
    /// @return @c SurfaceVector at given bin. Copy of all bins selected
    const SurfaceVector&
    neighbors(const Vector3D& pos) const
    {
      auto loc = m_globalToLocal(pos);
      return m_neighborMap.at(m_grid.getGlobalBinIndex(loc));
    }

    /// @brief Returns the total size of the grid (including under/overflow
    /// bins)
    /// @return Size of the grid data structure
    size_t
    size() const
    {
      return m_grid.size();
    }

#ifdef DOXYGEN
    /// @brief Gets the center position of bin @c bin in global coordinates
    /// @param bin the global bin index
    /// @return The bin center
    Vector3D
    getBinCenter(size_t bin) const;
#endif

    /// @cond
    template <size_t D = DIM, std::enable_if_t<D != 1, int> = 0>
    Vector3D
    getBinCenter(size_t bin) const
    {
      return m_localToGlobal(ActsVectorD<DIM>(
          m_grid.getBinCenter(m_grid.getLocalBinIndices(bin)).data()));
    }

    template <size_t D = DIM, std::enable_if_t<D == 1, int> = 0>
    Vector3D
    getBinCenter(size_t bin) const
    {
      point_t pos = m_grid.getBinCenter(m_grid.getLocalBinIndices(bin));
      return m_localToGlobal(pos);
    }
    /// @endcond

    /// @brief Returns copies of the axes used in the grid as @c AnyAxis
    /// @return The axes
    /// @note This returns copies. Use for introspection and querying.
    std::vector<concept::AnyAxis<>>
    getAxes() const
    {
      auto arr = m_grid.getAxes();
      return std::vector<concept::AnyAxis<>>(arr.begin(), arr.end());
    }

    /// @brief Get the number of dimensions of the grid.
    /// @return number of dimensions
    static constexpr size_t
    dimensions()
    {
      return DIM;
    }

    /// @brief Checks if global bin is valid
    /// @param bin the global bin index
    /// @return bool if the bin is valid
    /// @note Valid means that the index points to a bin which is not a under
    ///       or overflow bin or out of range in any axis.
    bool
    isValidBin(size_t bin) const
    {
      std::array<size_t, DIM> indices = m_grid.getLocalBinIndices(bin);
      std::array<size_t, DIM> nBins   = m_grid.getNBins();
      for (size_t i = 0; i < indices.size(); ++i) {
        size_t idx = indices.at(i);
        if (idx <= 0 || idx >= nBins.at(i) + 1) return false;
      }

      return true;
    }

  private:
    void
    populateNeighborCache()
    {
      // calculate neighbors for every bin and store in map
      for (size_t i = 0; i < m_grid.size(); i++) {
        if (!isValidBin(i)) continue;
        typename Grid_t::index_t loc  = m_grid.getLocalBinIndices(i);
        std::set<size_t> neighborIdxs = m_grid.neighborHoodIndices(loc, 1u);
        std::vector<const Surface*>& neighbors = m_neighborMap.at(i);
        neighbors.clear();

        for (const auto& idx : neighborIdxs) {
          const std::vector<const Surface*>& binContent = m_grid.at(idx);
          std::copy(binContent.begin(),
                    binContent.end(),
                    std::back_inserter(neighbors));
        }
      }
    }

    std::function<point_t(const Vector3D&)> m_globalToLocal;
    std::function<Vector3D(const point_t&)> m_localToGlobal;
    Grid_t                                  m_grid;
    std::vector<SurfaceVector>              m_neighborMap;
  };

  /// @brief Lookup implementation which wraps one element and always returns
  /// this
  ///        element when lookup is called
  struct SingleElementLookup
  {

    /// @brief Default constructor.
    /// @param element the one and only element.
    SingleElementLookup(SurfaceVector::value_type element)
      : m_element({element})
    {
    }

    /// @brief Lookup, always returns @c element
    /// @param pos is ignored
    /// @return reference to vector containing only @c element
    SurfaceVector&
    lookup(const Vector3D&)
    {
      return m_element;
    }

    /// @brief Lookup, always returns @c element
    /// @param pos is ignored
    /// @return reference to vector containing only @c element
    const SurfaceVector&
    lookup(const Vector3D&) const
    {
      return m_element;
    }

    /// @brief Lookup, always returns @c element
    /// @param bin is ignored
    /// @return reference to vector containing only @c element
    SurfaceVector& lookup(size_t) { return m_element; }

    /// @brief Lookup, always returns @c element
    /// @param bin is ignored
    /// @return reference to vector containing only @c element
    const SurfaceVector& lookup(size_t) const { return m_element; }

    /// @brief Lookup, always returns @c element
    /// @param pos is ignored
    /// @return reference to vector containing only @c element
    const SurfaceVector&
    neighbors(const Vector3D&) const
    {
      return m_element;
    }

    /// @brief returns 1
    /// @return 1
    size_t
    size() const
    {
      return 1;
    }

    /// @brief Gets the bin center, but always returns (0, 0, 0)
    /// @param bin is ignored
    /// @return (0, 0, 0)
    Vector3D getBinCenter(size_t) const { return Vector3D(0, 0, 0); }

    /// @brief Returns an empty vector of @c AnyAxis
    /// @return empty vector
    std::vector<concept::AnyAxis<>>
    getAxes() const
    {
      return {};
    }

    /// @brief Get the number of dimensions
    /// @return always 0
    static constexpr size_t
    dimensions()
    {
      return 0;
    }

    /// @brief Comply with concept and provide fill method
    /// @note Does nothing
    void
    fill(const SurfaceVector&)
    {
    }

    /// @brief Comply with concept and provide completeBinning method
    /// @note Does nothing
    size_t
    completeBinning(const SurfaceVector&)
    {
      return 0;
    }

    /// @brief Returns if the bin is valid (it is)
    /// @param bin is ignored
    /// @return always true
    static constexpr bool isValidBin(size_t) { return true; }

  private:
    SurfaceVector m_element;
  };

  /// @brief Default constructor which takes a @c SurfaceLookup and a vector of
  /// surfaces
  /// @param gridLookup The grid storage. @c SurfaceArray does not fill it on
  /// its own
  /// @param surfaces The input vector of surfaces. This is only for
  /// bookkeeping, so we can ask
  ///                 it for 'all contained surfaces'
  SurfaceArray(AnySurfaceGridLookup_t gridLookup, SurfaceVector surfaces)
    : m_gridLookup(std::move(gridLookup)), m_surfaces(surfaces)
  {
  }

  /// @brief Convenience constructor for single element mode. Uses the @c
  /// SingleElementLookup
  /// @param srf The one and only surface
  SurfaceArray(const Surface* srf) : m_gridLookup(SingleElementLookup(srf)) {}

  /// @brief Get all surfaces in bin given by position.
  /// @param pos the lookup position
  /// @return reference to @c SurfaceVector contained in bin at that position
  SurfaceVector&
  at(const Vector3D& pos)
  {
    return m_gridLookup.lookup(pos);
  }

  /// @brief Get all surfaces in bin given by position @p pos.
  /// @param pos the lookup position
  /// @return const reference to @c SurfaceVector contained in bin at that
  /// position
  const SurfaceVector&
  at(const Vector3D& pos) const
  {
    return m_gridLookup.lookup(pos);
  }

  /// @brief Get all surfaces in bin given by global bin index @p bin.
  /// @param bin the global bin index
  /// @return reference to @c SurfaceVector contained in bin
  SurfaceVector&
  at(size_t bin)
  {
    return m_gridLookup.lookup(bin);
  }

  /// @brief Get all surfaces in bin given by global bin index.
  /// @param bin the global bin index
  /// @return const reference to @c SurfaceVector contained in bin
  const SurfaceVector&
  at(size_t bin) const
  {
    return m_gridLookup.lookup(bin);
  }

  /// @brief Get all surfaces in bin at @p pos and its neighbors
  /// @param pos The position to lookup as nominal
  /// @param size How many neighbors we want in each direction. (default: 1)
  /// @return Merged @c SurfaceVector of neighbors and nominal
  /// @note The @c SurfaceVector will be combined. For technical reasons, the
  ///       different bin content vectors have to be copied, so the resulting
  ///       vector contains copies.
  SurfaceVector
  neighbors(const Vector3D& pos) const
  {
    return m_gridLookup.neighbors(pos);
  }

  /// @brief Get the size of the underlying grid structure including
  /// under/overflow bins
  /// @return the size
  size_t
  size() const
  {
    return m_gridLookup.size();
  }

  /// @brief Get the center of the bin identified by global bin index @p bin
  /// @param bin the global bin index
  /// @return Center position of the bin in global coordinates
  Vector3D
  getBinCenter(size_t bin)
  {
    return m_gridLookup.getBinCenter(bin);
  }

  /// @brief Get all surfaces attached to this @c SurfaceArray
  /// @return Reference to @c SurfaceVector containing all surfaces
  /// @note This does not reflect the actual state of the grid. It only
  ///       returns what was given in the constructor, without any checks
  ///       if that is actually whats in the grid.
  const SurfaceVector&
  surfaces() const
  {
    return m_surfaces;
  }

  /// @brief Get vector of axes spanning the grid as @c AnyAxis
  /// @return vector of @c AnyAxis
  /// @note The axes in the vector are copies. Only use for introspection and
  ///       querying.
  std::vector<concept::AnyAxis<>>
  getAxes() const
  {
    return m_gridLookup.getAxes();
  }

  /// @brief Checks if global bin is valid
  /// @param bin the global bin index
  /// @return bool if the bin is valid
  /// @note Valid means that the index points to a bin which is not a under
  ///       or overflow bin or out of range in any axis.
  bool
  isValidBin(size_t bin) const
  {
    return m_gridLookup.isValidBin(bin);
  }

  /// @brief String representation of this @c SurfaceArray
  /// @param sl Output stream to write to
  /// @return the output stream given as @p sl
  std::ostream&
  dump(std::ostream& sl) const
  {
    sl << "SurfaceArray:" << std::endl;
    sl << " - no surfaces: " << m_surfaces.size() << std::endl;
    sl << " - grid dim:    " << m_gridLookup.dimensions() << std::endl;

    auto axes = m_gridLookup.getAxes();

    for (size_t j = 0; j < axes.size(); ++j) {
      detail::AxisBoundaryType bdt = axes.at(j).getBoundaryType();
      sl << " - axis " << (j + 1) << std::endl;
      sl << "   - boundary type: ";
      if (bdt == detail::AxisBoundaryType::Open) sl << "open";
      if (bdt == detail::AxisBoundaryType::Bound) sl << "bound";
      if (bdt == detail::AxisBoundaryType::Closed) sl << "closed";
      sl << std::endl;
      sl << "   - type: "
         << (axes.at(j).isEquidistant() ? "equidistant" : "variable")
         << std::endl;
      sl << "   - n bins: " << axes.at(j).getNBins() << std::endl;
      sl << "   - bin edges: [ ";
      auto binEdges = axes.at(j).getBinEdges();
      for (size_t i = 0; i < binEdges.size(); ++i) {
        if (i > 0) sl << ", ";
        sl << binEdges.at(i);
      }
      sl << " ]" << std::endl;
    }
    return sl;
  }

private:
  AnySurfaceGridLookup_t m_gridLookup;
  SurfaceVector          m_surfaces;
};

}  // namespace Acts

#endif  // ACTS_SURFACES_SURFACEARRAY_H
