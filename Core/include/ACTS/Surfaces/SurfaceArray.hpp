// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
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

template <typename T, typename U>
constexpr bool
eqTypes()
{
  return std::is_same<T, U>::value;
}

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

public:
  using AnySurfaceGridLookup_t = concept::AnySurfaceGridLookup<SurfaceVector>;

  /// @brief Lookup helper which encapsulates a @c Grid
  template <size_t DIM>
  struct SurfaceGridLookup
  {
  public:

    /// @brief Default constructor
    ///
    /// @param globalToLocal Callable that converts from global to local
    /// @param localToGlobal Callable that converts from local to global
    /// @param grid The grid data structur. Will be type-erased at this point
    SurfaceGridLookup(
        std::function<ActsVectorD<DIM>(const Vector3D&)> globalToLocal,
        std::function<Vector3D(const ActsVectorD<DIM>&)> localToGlobal,
        AnyGrid_t<ActsVectorD<DIM>, DIM>                 grid)
      : m_globalToLocal(std::move(globalToLocal))
      , m_localToGlobal(std::move(localToGlobal))
      , m_grid(std::move(grid))
    {
    }

    /// @brief Fill provided surfaces into the contained @c Grid.
    ///
    /// This is done by iterating, accessing the binningPosition, lookup
    /// and append.
    /// 
    /// @param surfaces Input surface pointers
    void
    fill(const SurfaceVector surfaces)
    {
      for (const auto& srf : surfaces) {
        Vector3D pos = srf->binningPosition(binR);
        lookup(pos).push_back(srf);
      }
    }

    /// @brief Performs lookup at @c pos and returns bin content as reference
    /// @param pos Lookup position
    /// @return @c SurfaceVector at given bin
    SurfaceVector&
    lookup(const Vector3D& pos)
    {
      return m_grid.at(m_globalToLocal(pos));
    }

    /// @brief Performs lookup at @c pos and returns bin content as const reference
    /// @param pos Lookup position
    /// @return @c SurfaceVector at given bin
    const SurfaceVector&
    lookup(const Vector3D& pos) const
    {
      return m_grid.at(m_globalToLocal(pos));
    }

    /// @brief Performs lookup at global bin and returns bin content as reference
    /// @param bin Global lookup bin
    /// @return @c SurfaceVector at given bin
    SurfaceVector&
    lookup(size_t bin)
    {
      return m_grid.at(bin);
    }

    /// @brief Performs lookup at global bin and returns bin content as const reference
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
    /// @param size Optional number of neighbors, default is 1 (= next neighbor)
    /// @return @c SurfaceVector at given bin. Copy of all bins selected
    ///
    /// @note The resulting SurfaceVectors for each bin will be merged. Sine @c std::vector
    ///       cannot store references, the source bin content vector entries have to be copied.
    ///       This should be fine for the pointer value type used here.
    std::vector<const Surface*>
    neighbors(const Vector3D& pos, size_t size = 1) const
    {
      auto loc = m_globalToLocal(pos);
      std::set<size_t> neighborIdxs = m_grid.neighborHoodIndices(loc, size);
      std::vector<std::vector<const Surface*>> binContents
          = m_grid.atBins(neighborIdxs);
      std::vector<const Surface*> out;

      for (const auto& binContent : binContents) {
        // move b/c bincontents is already a copy and not needed after this
        std::move(
            binContent.begin(), binContent.end(), std::back_inserter(out));
      }
      return out;
    }

    /// @brief Returns the total size of the grid (including under/overflow bins)
    /// @return Size of the grid data structure
    size_t
    size() const
    {
      return m_grid.size();
    }

    /// @brief Gets the center position of bin @c bin in global coordinates
    /// @param bin the global bin index
    /// @return The bin center
    Vector3D
    getBinCenter(size_t bin) const
    {
      return m_localToGlobal(ActsVectorD<DIM>(
          m_grid.getBinCenter(m_grid.getLocalBinIndices(bin)).data()));
    }

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
    std::function<ActsVectorD<DIM>(const Vector3D&)> m_globalToLocal;
    std::function<Vector3D(const ActsVectorD<DIM>&)> m_localToGlobal;
    AnyGrid_t<ActsVectorD<DIM>, DIM>                 m_grid;
  };

  /// @brief Lookup implementation which wraps one element and always returns this
  ///        element when lookup is called
  struct SingleElementLookup
  {

    /// @brief Default constructor.
    /// @param element the one and only element.
    SingleElementLookup(SurfaceVector::value_type element) : m_element({element})
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
    SurfaceVector
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

    /// @brief Returns if the bin is valid (it is)
    /// @param bin is ignored
    /// @return always true
    static constexpr bool
    isValidBin(size_t)
    {
      return true;
    }

  private:
    SurfaceVector m_element;
  };

  // useful typedefs
  using SurfaceGridLookup1D = SurfaceGridLookup<1>;
  using SurfaceGridLookup2D = SurfaceGridLookup<2>;
  using SurfaceGridLookup3D = SurfaceGridLookup<3>;

  /// @brief Default constructor which takes a @c SurfaceLookup and a vector of surfaces
  /// @param gridLookup The grid storage. @c SurfaceArray does not fill it on its own
  /// @param surfaces The input vector of surfaces. This is only for bookkeeping, so we can ask
  ///                 it for 'all contained surfaces'
  SurfaceArray(AnySurfaceGridLookup_t gridLookup, SurfaceVector surfaces)
    : m_gridLookup(std::move(gridLookup)), m_surfaces(surfaces)
  {
  }

  /// @brief Convenience constructor for single element mode. Uses the @c SingleElementLookup
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
  /// @return const reference to @c SurfaceVector contained in bin at that position
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
  neighbors(const Vector3D& pos, size_t size = 1) const
  {
    return m_gridLookup.neighbors(pos, size);
  }

  /// @brief Get the size of the underlying grid structure including under/overflow bins
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
      // auto boundaries = axes.boundaries();
      // std::cout << "axis is " << axis.isEquidistant() << std::endl;
      detail::AxisWrapping wrap = axes.at(j).getWrapping();
      std::cout << " - axis " << (j + 1) << std::endl;
      std::cout << " - wrapping: ";
      if (wrap == detail::AxisWrapping::UnderOverflow)
        std::cout << "under/overflow";
      if (wrap == detail::AxisWrapping::Open) std::cout << "open";
      if (wrap == detail::AxisWrapping::Closed) std::cout << "closed";
      std::cout << std::endl;
      std::cout << "   - bin edges: [ ";
      auto binEdges = axes.at(j).getBinEdges();
      for (size_t i = 0; i < binEdges.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << binEdges.at(i);
      }
      std::cout << " ]" << std::endl;
    }
    return sl;
  }

private:
  AnySurfaceGridLookup_t m_gridLookup;
  SurfaceVector             m_surfaces;
};

}  // namespace Acts

#endif  // ACTS_SURFACES_SURFACEARRAY_H
