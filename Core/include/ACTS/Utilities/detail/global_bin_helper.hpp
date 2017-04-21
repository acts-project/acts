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

    template <class... Axes>
    static void
    getGlobalBin(const std::array<size_t, sizeof...(Axes)>& localBins,
                 const std::tuple<Axes...>& axes,
                 size_t&                    bin,
                 size_t&                    area)
    {
      const auto& thisAxis = std::get<N>(axes);
      bin += area * localBins.at(N);
      // make sure to account for under-/overflow bins
      area *= (thisAxis.getNBins() + 2);
      global_bin_helper_impl<N + 1, MAX>::getGlobalBin(
          localBins, axes, bin, area);
    }

    template <class... Axes>
    static void
    getLocalBinIndices(size_t&                    bin,
                       const std::tuple<Axes...>& axes,
                       size_t&                    area,
                       std::array<size_t, sizeof...(Axes)>& indices)
    {
      const auto& thisAxis = std::get<N>(axes);
      // make sure to account for under-/overflow bins
      size_t new_area = area * (thisAxis.getNBins() + 2);
      global_bin_helper_impl<N + 1, MAX>::getLocalBinIndices(
          bin, axes, new_area, indices);
      indices.at(N) = bin / area;
      bin %= area;
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

    template <class... Axes>
    static void
    getGlobalBin(const std::array<size_t, sizeof...(Axes)>& localBins,
                 const std::tuple<Axes...>& axes,
                 size_t&                    bin,
                 size_t&                    area)
    {
      bin += area * localBins.at(MAX);
    }

    template <class... Axes>
    static void
    getLocalBinIndices(size_t&                    bin,
                       const std::tuple<Axes...>& axes,
                       size_t&                    area,
                       std::array<size_t, sizeof...(Axes)>& indices)
    {
      // make sure to account for under-/overflow bins
      indices.at(MAX) = bin / area;
      bin %= area;
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

    /// @brief determine global bin index from local indices along each axis
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    ///
    /// @param  [in] localBins local bin indices along each axis
    /// @param  [in] axes  actual axis objects spanning the grid
    /// @return global index for bin defined by the local bin indices
    ///
    /// @pre All local bin indices must be a valid index for the corresponding
    ///      axis (including the under-/overflow bin for this axis).
    template <class... Axes>
    static size_t
    getGlobalBin(const std::array<size_t, sizeof...(Axes)>& localBins,
                 const std::tuple<Axes...>& axes)
    {
      constexpr size_t MAX  = sizeof...(Axes)-1;
      size_t           area = 1;
      size_t           bin  = 0;

      global_bin_helper_impl<0, MAX>::getGlobalBin(localBins, axes, bin, area);

      return bin;
    }

    /// @brief determine local bin index for each axis from global bin index
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    ///
    /// @param  [in] bin  global bin index
    /// @param  [in] axes actual axis objects spanning the grid
    /// @return array with local bin indices along each axis (in same order as
    ///         given @c axes object)
    ///
    /// @note Local bin indices could be a under-/overflow bin along this axis.
    template <class... Axes>
    static std::array<size_t, sizeof...(Axes)>
    getLocalBinIndices(size_t bin, const std::tuple<Axes...>& axes)
    {
      constexpr size_t MAX  = sizeof...(Axes)-1;
      size_t           area = 1;
      std::array<size_t, sizeof...(Axes)> indices;

      global_bin_helper_impl<0, MAX>::getLocalBinIndices(
          bin, axes, area, indices);

      return indices;
    }
  };

}  // namespace detail

}  // namespace Acts
