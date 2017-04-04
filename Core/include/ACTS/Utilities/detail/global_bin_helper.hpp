// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <tuple>
#include <utility>

namespace Acts {

namespace detail {

  /// @cond
  template <size_t N, size_t MAX>
  struct global_bin_helper_impl;

  template <size_t N, size_t MAX>
  struct global_bin_helper_impl
  {
    template <class Point, class... Axes>
    static void
    getGlobalBin(const Point&               point,
                 const std::tuple<Axes...>& axes,
                 size_t&                    bin,
                 size_t&                    area)
    {
      const auto& thisAxis = std::get<N>(axes);
      bin += area * thisAxis.getBin(point[N]);
      // make sure to account for under-/overflow bins
      area *= (thisAxis.getNBins() + 2);
      global_bin_helper_impl<N + 1, MAX>::getGlobalBin(point, axes, bin, area);
    }
  };

  template <size_t MAX>
  struct global_bin_helper_impl<MAX, MAX>
  {
    template <class Point, class... Axes>
    static void
    getGlobalBin(const Point&               point,
                 const std::tuple<Axes...>& axes,
                 size_t&                    bin,
                 size_t&                    area)
    {
      const auto& thisAxis = std::get<MAX>(axes);
      bin += area * thisAxis.getBin(point[MAX]);
    }
  };
  /// @endcond

  /// @brief determine global bin index in grid defined by a set of axes
  struct global_bin_helper
  {
    /// @brief determine global bin index in grid defined by a set of axes
    ///
    /// @tparam Point any type with point semantics supporting component access
    ///               through @c operator[]
    /// @tparam Axes parameter pack of axis types defining the grid
    ///
    /// @param  [in] point point to look up in the grid
    /// @param  [in] axes  actual axis objects spanning the grid
    /// @return global index for bin containing the given point
    ///
    /// @note This could be a under-/overflow bin along one or more axes.
    template <class Point, class... Axes>
    static size_t
    getGlobalBin(const Point& point, const std::tuple<Axes...>& axes)
    {
      constexpr size_t MAX  = sizeof...(Axes)-1;
      size_t           area = 1;
      size_t           bin  = 0;

      global_bin_helper_impl<0, MAX>::getGlobalBin(point, axes, bin, area);

      return bin;
    }
  };

}  // namespace detail

}  // namespace Acts
