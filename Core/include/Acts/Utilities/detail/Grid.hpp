// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <numeric>
#include <set>
#include <tuple>
#include <type_traits>
#include <vector>
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/Interpolation.hpp"
#include "Acts/Utilities/detail/grid_helper.hpp"

namespace Acts {

namespace detail {

  /// @brief class for describing a regular multi-dimensional grid
  ///
  /// @tparam T    type of values stored inside the bins of the grid
  /// @tparam Axes parameter pack of axis types defining the grid
  ///
  /// Class describing a multi-dimensional, regular grid which can store objects
  /// in its multi-dimensional bins. Bins are hyper-boxes and can be accessed
  /// either by global bin index, local bin indices or position.
  ///
  /// @note @c T must be default-constructible.
  template <typename T, class... Axes>
  class Grid final
  {
    /// number of dimensions of the grid
    static constexpr size_t DIM = sizeof...(Axes);

  public:
    /// type of values stored
    using value_type = T;
    /// reference type to values stored
    using reference = value_type&;
    /// constant reference type to values stored
    using const_reference = const value_type&;
    /// type for points in d-dimensional grid space
    using point_t = std::array<double, DIM>;
    /// index type using local bin indices along each axis
    using index_t = std::array<size_t, DIM>;

    /// @brief default constructor
    ///
    /// @param [in] axes actual axis objects spanning the grid
    Grid(std::tuple<Axes...> axes) : m_axes(std::move(axes))
    {
      m_values.resize(size());
    }

    /// @brief access value stored in bin for a given point
    ///
    /// @tparam Point any type with point semantics supporting component access
    ///               through @c operator[]
    /// @param [in] point point used to look up the corresponding bin in the
    ///                   grid
    /// @return reference to value stored in bin containing the given point
    ///
    /// @pre The given @c Point type must represent a point in d (or higher)
    ///      dimensions where d is dimensionality of the grid.
    ///
    /// @note The look-up considers under-/overflow bins along each axis.
    ///       Therefore, the look-up will never fail.
    //
    template <class Point>
    reference
    at(const Point& point)
    {
      return m_values.at(getGlobalBinIndex(point));
    }

    /// @brief access value stored in bin for a given point
    ///
    /// @tparam Point any type with point semantics supporting component access
    ///               through @c operator[]
    /// @param [in] point point used to look up the corresponding bin in the
    ///                   grid
    /// @return const-reference to value stored in bin containing the given
    ///         point
    ///
    /// @pre The given @c Point type must represent a point in d (or higher)
    ///      dimensions where d is dimensionality of the grid.
    ///
    /// @note The look-up considers under-/overflow bins along each axis.
    ///       Therefore, the look-up will never fail.
    template <class Point>
    const_reference
    at(const Point& point) const
    {
      return m_values.at(getGlobalBinIndex(point));
    }

    /// @brief access value stored in bin with given global bin number
    ///
    /// @param  [in] bin global bin number
    /// @return reference to value stored in bin containing the given
    ///         point
    reference
    at(size_t bin)
    {
      return m_values.at(bin);
    }

    /// @brief access value stored in bin with given global bin number
    ///
    /// @param  [in] bin global bin number
    /// @return const-reference to value stored in bin containing the given
    ///         point
    const_reference
    at(size_t bin) const
    {
      return m_values.at(bin);
    }

    /// @brief access value stored in bin with given local bin numbers
    ///
    /// @param  [in] localBins local bin indices along each axis
    /// @return reference to value stored in bin containing the given
    ///         point
    ///
    /// @pre All local bin indices must be a valid index for the corresponding
    ///      axis (including the under-/overflow bin for this axis).
    reference
    at(const index_t& localBins)
    {
      return m_values.at(getGlobalBinIndex(localBins));
    }

    /// @brief access value stored in bin with given local bin numbers
    ///
    /// @param  [in] localBins local bin indices along each axis
    /// @return const-reference to value stored in bin containing the given
    ///         point
    ///
    /// @pre All local bin indices must be a valid index for the corresponding
    ///      axis (including the under-/overflow bin for this axis).
    const_reference
    at(const index_t& localBins) const
    {
      return m_values.at(getGlobalBinIndex(localBins));
    }

    /// @brief get global bin indices for closest points on grid
    ///
    /// @tparam Point any type with point semantics supporting component access
    ///               through @c operator[]
    /// @param [in] position point of interest
    /// @return set of global bin indices for bins whose lower-left corners are
    ///         the closest points on the grid to the given point
    ///
    /// @pre The given @c Point type must represent a point in d (or higher)
    ///      dimensions where d is dimensionality of the grid. It must lie
    ///      within the grid range (i.e. not within a under-/overflow bin).
    template <class Point>
    std::set<size_t>
    closestPointsIndices(const Point& position) const
    {
      return grid_helper::closestPointsIndices(getGlobalBinIndex(position),
                                               m_axes);
    }

    /// @brief dimensionality of grid
    ///
    /// @return number of axes spanning the grid
    static constexpr size_t
    dimensions()
    {
      return DIM;
    }

    /// @brief get center position of bin with given local bin numbers
    ///
    /// @param  [in] localBins local bin indices along each axis
    /// @return center position of bin
    ///
    /// @pre All local bin indices must be a valid index for the corresponding
    ///      axis (excluding the under-/overflow bins for each axis).
    std::array<double, DIM>
    getBinCenter(const index_t& localBins) const
    {
      return grid_helper::getBinCenter(localBins, m_axes);
    }

    /// @brief determine global index for bin containing the given point
    ///
    /// @tparam Point any type with point semantics supporting component access
    ///               through @c operator[]
    ///
    /// @param  [in] point point to look up in the grid
    /// @return global index for bin containing the given point
    ///
    /// @pre The given @c Point type must represent a point in d (or higher)
    ///      dimensions where d is dimensionality of the grid.
    /// @note This could be a under-/overflow bin along one or more axes.
    template <class Point>
    size_t
    getGlobalBinIndex(const Point& point) const
    {
      return grid_helper::getGlobalBin(point, m_axes);
    }

    /// @brief determine global bin index from local bin indices along each axis
    ///
    /// @param  [in] localBins local bin indices along each axis
    /// @return global index for bin defined by the local bin indices
    ///
    /// @pre All local bin indices must be a valid index for the corresponding
    ///      axis (including the under-/overflow bin for this axis).
    size_t
    getGlobalBinIndex(const index_t& localBins) const
    {
      return grid_helper::getGlobalBin(localBins, m_axes);
    }

    /// @brief determine local bin index for each axis from global bin index
    ///
    /// @param  [in] bin global bin index
    /// @return array with local bin indices along each axis (in same order as
    ///         given @c axes object)
    ///
    /// @note Local bin indices can contain under-/overflow bins along the
    ///       corresponding axis.
    index_t
    getLocalBinIndices(size_t bin) const
    {
      return grid_helper::getLocalBinIndices(bin, m_axes);
    }

    /// @brief retrieve lower-left bin edge from set of local bin indices
    ///
    /// @param  [in] localBins local bin indices along each axis
    /// @return generalized lower-left bin edge position
    ///
    /// @pre @c localBins must only contain valid bin indices (excluding
    ///      underflow bins).
    point_t
    getLowerLeftBinEdge(const index_t& localBins) const
    {
      return grid_helper::getLowerLeftBinEdge(localBins, m_axes);
    }

    /// @brief get number of bins along each specific axis
    ///
    /// @return array giving the number of bins along all axes
    ///
    /// @note Not including under- and overflow bins
    index_t
    getNBins() const
    {
      return grid_helper::getNBins(m_axes);
    }

    /// @brief get the minimum value of all axes of one grid
    ///
    /// @return array returning the minima of all given axes
    point_t
    getMin() const
    {
      return grid_helper::getMin(m_axes);
    }

    /// @brief get the maximum value of all axes of one grid
    ///
    /// @return array returning the maxima of all given axes
    point_t
    getMax() const
    {
      return grid_helper::getMax(m_axes);
    }

    /// @brief retrieve upper-right bin edge from set of local bin indices
    ///
    /// @param  [in] localBins local bin indices along each axis
    /// @return generalized upper-right bin edge position
    ///
    /// @pre @c localBins must only contain valid bin indices (excluding
    ///      overflow bins).
    point_t
    getUpperRightBinEdge(const index_t& localBins) const
    {
      return grid_helper::getUpperRightBinEdge(localBins, m_axes);
    }

    /// @brief set all overflow and underflow bins to a certain value
    ///
    /// @param [in] value value to be inserted in every overflow and underflow
    ///                   bin of the grid.
    ///
    void
    setExteriorBins(const value_type& value)
    {
      for (size_t index : grid_helper::exteriorBinIndices(m_axes)) {
        at(index) = value;
      }
    }

    /// @brief interpolate grid values to given position
    ///
    /// @tparam Point type specifying geometric positions
    /// @tparam U     dummy template parameter identical to @c T
    ///
    /// @param [in] point location to which to interpolate grid values. The
    ///                   position must be within the grid dimensions and not
    ///                   lie in an under-/overflow bin along any axis.
    ///
    /// @return interpolated value at given position
    ///
    /// @pre The given @c Point type must represent a point in d (or higher)
    ///      dimensions where d is dimensionality of the grid.
    ///
    /// @note This function is available only if the following conditions are
    /// fulfilled:
    /// - Given @c U and @c V of value type @c T as well as two @c double @c a
    /// and @c b, then the following must be a valid expression <tt>a * U + b *
    /// V</tt> yielding an object which is (implicitly) convertible to @c T.
    /// - @c Point must represent a d-dimensional position and support
    /// coordinate access using @c operator[] which should return a @c double
    /// (or a value which is implicitly convertible). Coordinate indices must
    /// start at 0.
    /// @note Bin values are interpreted as being the field values at the
    /// lower-left corner of the corresponding hyper-box.
    template <class Point,
              typename U = T,
              typename
              = std::enable_if_t<can_interpolate<Point,
                                                 std::array<double, DIM>,
                                                 std::array<double, DIM>,
                                                 U>::value>>
    T
    interpolate(const Point& point) const
    {
      // there are 2^DIM corner points used during the interpolation
      constexpr size_t nCorners = 1 << DIM;

      // construct vector of pairs of adjacent bin centers and values
      std::array<value_type, nCorners> neighbors;

      // get local indices for current bin
      // value of bin is interpreted as being the field value at its lower left
      // corner
      const auto& llIndices = getLocalBinIndices(getGlobalBinIndex(point));

      // get local indices for "upper right" bin (i.e. all local indices
      // incremented by one but avoid overflows)
      auto urIndices = grid_helper::getUpperRightBinIndices(llIndices, m_axes);

      // get global indices for all surrounding corner points
      const auto& closestIndices = closestPointsIndices(point);

      // get values on grid points
      size_t i = 0;
      for (size_t index : closestIndices) {
        // In case it is the last bin return the last value as neighbour
        neighbors.at(i++) = (index < size()) ? at(index) : at(size() - 1);
      }

      return Acts::interpolate(point,
                               getLowerLeftBinEdge(llIndices),
                               getLowerLeftBinEdge(urIndices),
                               neighbors);
    }

    /// @brief check whether given point is inside grid limits
    ///
    /// @return @c true if \f$\text{xmin_i} \le x_i < \text{xmax}_i \forall i=0,
    ///         \dots, d-1\f$, otherwise @c false
    ///
    /// @pre The given @c Point type must represent a point in d (or higher)
    ///      dimensions where d is dimensionality of the grid.
    ///
    /// @post If @c true is returned, the global bin containing the given point
    ///       is a valid bin, i.e. it is neither a underflow nor an overflow bin
    ///       along any axis.
    template <class Point>
    bool
    isInside(const Point& position) const
    {
      return grid_helper::isInside(position, m_axes);
    }

    /// @brief get global bin indices for neighborhood
    ///
    /// @param [in] localBins center bin defined by local bin indices along each
    ///                       axis
    /// @param [in] size      size of neighborhood determining how many adjacent
    ///                       bins along each axis are considered
    /// @return set of global bin indices for all bins in neighborhood
    ///
    /// @note Over-/underflow bins are included in the neighborhood.
    /// @note The @c size parameter sets the range by how many units each local
    ///       bin index is allowed to be varied. All local bin indices are
    ///       varied independently, that is diagonal neighbors are included.
    ///       Ignoring the truncation of the neighborhood size reaching beyond
    ///       over-/underflow bins, the neighborhood is of size \f$2 \times
    ///       \text{size}+1\f$ along each dimension.
    std::set<size_t>
    neighborHoodIndices(const index_t& localBins, size_t size = 1u) const
    {
      return grid_helper::neighborHoodIndices(localBins, size, m_axes);
    }

    /// @brief get global bin indices for neighborhood of bin identified by @p
    /// pos
    /// @param pos position around which to look
    /// @param size how many neighbors (defaul: 1)
    /// @return set of global bin indices pointing to the neighbors
    template <class Point>
    std::set<size_t>
    neighborHoodIndices(const Point& pos, size_t size = 1u) const
    {
      const size_t  bin       = getGlobalBinIndex(pos);
      const index_t localBins = getLocalBinIndices(bin);
      return grid_helper::neighborHoodIndices(localBins, size, m_axes);
    }

    /// @brief total number of bins
    ///
    /// @return total number of bins in the grid
    ///
    /// @note This number contains under-and overflow bins along all axes.
    size_t
    size() const
    {
      index_t nBinsArray = getNBins();
      // add under-and overflow bins for each axis and multiply all bins
      return std::accumulate(
          nBinsArray.begin(),
          nBinsArray.end(),
          1,
          [](const size_t& a, const size_t& b) { return a * (b + 2); });
    }

    std::array<const IAxis*, DIM>
    getAxes() const
    {
      return grid_helper::getAxes(m_axes);
    }

  private:
    /// set of axis defining the multi-dimensional grid
    std::tuple<Axes...> m_axes;
    /// linear value store for each bin
    std::vector<T> m_values;
  };
}  // namespace detail

}  // namespace Acts
