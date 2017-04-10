// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <tuple>
#include "ACTS/Utilities/detail/get_nbins_helper.hpp"
#include "ACTS/Utilities/detail/global_bin_helper.hpp"

namespace Acts {

namespace Test {
  template <typename T, class... Axes>
  struct GridTester;
}

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
    /// unit test helper class
    friend struct Test::GridTester<T, Axes...>;

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

    /// @brief dimensionality of grid
    ///
    /// @return number of axes spanning the grid
    static constexpr size_t
    dimension()
    {
      return DIM;
    }

    /// @brief total number of bins
    ///
    /// @return total number of bins in the grid
    ///
    /// @note This number contains under-/overflow bins along one or more axes.
    size_t
    size() const
    {
      return get_nbins_helper::getNBins(m_axes);
    }

  private:
    /// @brief determine global bin index in grid defined by a set of axes
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

    /// @brief determine local bin index for each axis from global bin index
    ///
    /// @param  [in] bin global bin index
    /// @return array with local bin indices along each axis (in same order as
    ///         given @c axes object)
    ///
    /// @note Local bin indices could be a under-/overflow bin along this axis.
    std::array<size_t, DIM>
    getLocalBinIndices(size_t bin) const
    {
      return global_bin_helper::getLocalBinIndices(bin, m_axes);
    }

    /// set of axis defining the multi-dimensional grid
    std::tuple<Axes...> m_axes;
    /// linear value store for each bin
    std::vector<T> m_values;
  };
}  // namespace detail

}  // namespace Acts
