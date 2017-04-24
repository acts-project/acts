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
#include <utility>

namespace Acts {

namespace detail {

  /// @cond
  /// @brief helper struct to calculate number of bins inside a grid
  ///
  /// @tparam N number of axes to consider
  template <size_t N>
  struct grid_bins_helper_impl;

  template <size_t N>
  struct grid_bins_helper_impl
  {
    template <class... Axes>
    static void
    getBinCenter(std::array<double, sizeof...(Axes)>&       center,
                 const std::array<size_t, sizeof...(Axes)>& localIndices,
                 const std::tuple<Axes...>& axes)
    {
      center.at(N) = std::get<N>(axes).getBinCenter(localIndices.at(N));
      grid_bins_helper_impl<N - 1>::getBinCenter(center, localIndices, axes);
    }

    template <class... Axes>
    static size_t
    getNBins(const std::tuple<Axes...>& axes)
    {
      // by convention getNBins does not include under-/overflow bins
      size_t thisAxisNBins = std::get<N>(axes).getNBins() + 2;
      return thisAxisNBins * grid_bins_helper_impl<N - 1>::getNBins(axes);
    }
  };

  template <>
  struct grid_bins_helper_impl<0u>
  {
    template <class... Axes>
    static void
    getBinCenter(std::array<double, sizeof...(Axes)>&       center,
                 const std::array<size_t, sizeof...(Axes)>& localIndices,
                 const std::tuple<Axes...>& axes)
    {
      center.at(0u) = std::get<0u>(axes).getBinCenter(localIndices.at(0u));
    }

    template <class... Axes>
    static size_t
    getNBins(const std::tuple<Axes...>& axes)
    {
      // by convention getNBins does not include under-/overflow bins
      size_t thisAxisNBins = std::get<0u>(axes).getNBins() + 2;
      return thisAxisNBins;
    }
  };
  /// @endcond

  /// @brief calculate total number of bins in a grid defined by a set of axes
  struct grid_bins_helper
  {
    /// @brief retrieve bin center from set of local bin indices
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    /// @param  [in] localIndices local bin indices along each axis
    /// @param  [in] axes         actual axis objects spanning the grid
    /// @return center position of bin
    ///
    /// @pre @c localIndices must only contain valid bin indices (i.e. excluding
    ///      under-/overflow bins).
    /// @pre The dimensions of @c localIndices and @c axes must match.
    template <class... Axes>
    static std::array<double, sizeof...(Axes)>
    getBinCenter(const std::array<size_t, sizeof...(Axes)>& localIndices,
                 const std::tuple<Axes...>& axes)
    {
      std::array<double, sizeof...(Axes)> center;
      constexpr size_t MAX = sizeof...(Axes)-1;
      grid_bins_helper_impl<MAX>::getBinCenter(center, localIndices, axes);

      return center;
    }

    /// @brief calculate total number of bins in a grid defined by a set of axes
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    /// @param  [in] axes actual axis objects spanning the grid
    /// @return total number of bins in the grid
    ///
    /// @note This includes under-/overflow bins along each axis.
    template <class... Axes>
    static size_t
    getNBins(const std::tuple<Axes...>& axes)
    {
      constexpr size_t MAX = sizeof...(Axes)-1;
      return grid_bins_helper_impl<MAX>::getNBins(axes);
    }
  };

}  // namespace detail

}  // namespace Acts
