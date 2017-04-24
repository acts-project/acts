// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <ACTS/Utilities/detail/grid_bins_helper.hpp>
#include <array>
#include <tuple>
#include "ACTS/Utilities/detail/global_bin_helper.hpp"

namespace Acts {

namespace detail {

  /// @brief class for describing a regular multi-dimensional grid
  ///
  /// @tparam T    type of values stored inside the bins of the grid
  /// @tparam Axes parameter pack of axis types defining the grid
  ///
  /// Class describing a multi-dimensional, regular grid which can store objects
  /// in its multi-dimensional bins. Bins are hyper-cubes and can be accessed
  /// either by global bin index, local bin indices or position.
  ///
  /// @note @c T must be default-constructible.
  template <typename T, class... Axes>
  class Grid final
  {
    /// number of dimensions of the grid
    static constexpr size_t DIM = sizeof...(Axes);

  public:
    typedef T                 value_type;
    typedef value_type&       reference;
    typedef const value_type& const_reference;

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
    reference
    at(const std::array<size_t, DIM>& localBins)
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
    at(const std::array<size_t, DIM>& localBins) const
    {
      return m_values.at(getGlobalBinIndex(localBins));
    }

    /// @brief dimensionality of grid
    ///
    /// @return number of axes spanning the grid
    static constexpr size_t
    dimension()
    {
      return DIM;
    }

    /// @brief get center position of bin with given global bin number
    ///
    /// @param  [in] bin global bin number
    /// @return center position of bin
    ///
    /// @pre The specified global bin must not be an under- or overflow bin
    ///      along any axis.
    std::array<double, DIM>
    getBinCenter(size_t bin) const
    {
      const auto& localIndices = getLocalBinIndices(bin);
      return getBinCenter(localIndices);
    }

    /// @brief get center position of bin with given local bin numbers
    ///
    /// @param  [in] localBins local bin indices along each axis
    /// @return center position of bin
    ///
    /// @pre All local bin indices must be a valid index for the corresponding
    ///      axis (excluding the under-/overflow bins for each axis).
    std::array<double, DIM>
    getBinCenter(const std::array<size_t, DIM>& localBins) const
    {
      return grid_bins_helper::getBinCenter(localBins, m_axes);
    }

    /// @brief determine global index for bin containing the given point
    ///
    /// @tparam Point any type with point semantics supporting component access
    ///               through @c operator[]
    ///
    /// @param  [in] point point to look up in the grid
    /// @return global index for bin containing the given point
    ///
    /// @note This could be a under-/overflow bin along one or more axes.
    template <class Point>
    size_t
    getGlobalBinIndex(const Point& point) const
    {
      return global_bin_helper::getGlobalBin(point, m_axes);
    }

    /// @brief determine global bin index from local bin indices along each axis
    ///
    /// @param  [in] localBins local bin indices along each axis
    /// @return global index for bin defined by the local bin indices
    ///
    /// @pre All local bin indices must be a valid index for the corresponding
    ///      axis (including the under-/overflow bin for this axis).
    size_t
    getGlobalBinIndex(const std::array<size_t, DIM>& localBins) const
    {
      return global_bin_helper::getGlobalBin(localBins, m_axes);
    }

    /// @brief determine local bin index for each axis from global bin index
    ///
    /// @param  [in] bin global bin index
    /// @return array with local bin indices along each axis (in same order as
    ///         given @c axes object)
    ///
    /// @note Local bin indices can contain under-/overflow bins along the
    ///       corresponding axis.
    std::array<size_t, DIM>
    getLocalBinIndices(size_t bin) const
    {
      return global_bin_helper::getLocalBinIndices(bin, m_axes);
    }

    /// @brief total number of bins
    ///
    /// @return total number of bins in the grid
    ///
    /// @note This number contains under-and overflow bins along all axes.
    size_t
    size() const
    {
      return grid_bins_helper::getNBins(m_axes);
    }

  private:
    /// set of axis defining the multi-dimensional grid
    std::tuple<Axes...> m_axes;
    /// linear value store for each bin
    std::vector<T> m_values;
  };
}  // namespace detail

}  // namespace Acts
