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
    static void
    getLowerLeftBinEdge(std::array<double, sizeof...(Axes)>&       llEdge,
                        const std::array<size_t, sizeof...(Axes)>& localIndices,
                        const std::tuple<Axes...>& axes)
    {
      llEdge.at(N) = std::get<N>(axes).getBinLowerBound(localIndices.at(N));
      grid_bins_helper_impl<N - 1>::getLowerLeftBinEdge(
          llEdge, localIndices, axes);
    }

    template <class... Axes>
    static size_t
    getNBins(const std::tuple<Axes...>& axes)
    {
      // by convention getNBins does not include under-/overflow bins
      size_t thisAxisNBins = std::get<N>(axes).getNBins() + 2;
      return thisAxisNBins * grid_bins_helper_impl<N - 1>::getNBins(axes);
    }

    template <class... Axes>
    static void
    getUpperRightBinEdge(
        std::array<double, sizeof...(Axes)>&       urEdge,
        const std::array<size_t, sizeof...(Axes)>& localIndices,
        const std::tuple<Axes...>& axes)
    {
      urEdge.at(N) = std::get<N>(axes).getBinUpperBound(localIndices.at(N));
      grid_bins_helper_impl<N - 1>::getUpperRightBinEdge(
          urEdge, localIndices, axes);
    }

    template <class... Axes>
    static void
    getUpperRightBinIndices(std::array<size_t, sizeof...(Axes)>& localIndices,
                            const std::tuple<Axes...>& axes)
    {
      size_t thisAxisNBins = std::get<N>(axes).getNBins();
      localIndices.at(N)   = std::min(thisAxisNBins + 1, ++localIndices.at(N));
      grid_bins_helper_impl<N - 1>::getUpperRightBinIndices(localIndices, axes);
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
    static void
    getLowerLeftBinEdge(std::array<double, sizeof...(Axes)>&       llEdge,
                        const std::array<size_t, sizeof...(Axes)>& localIndices,
                        const std::tuple<Axes...>& axes)
    {
      llEdge.at(0u) = std::get<0u>(axes).getBinLowerBound(localIndices.at(0u));
    }

    template <class... Axes>
    static size_t
    getNBins(const std::tuple<Axes...>& axes)
    {
      // by convention getNBins does not include under-/overflow bins
      size_t thisAxisNBins = std::get<0u>(axes).getNBins() + 2;
      return thisAxisNBins;
    }

    template <class... Axes>
    static void
    getUpperRightBinEdge(
        std::array<double, sizeof...(Axes)>&       urEdge,
        const std::array<size_t, sizeof...(Axes)>& localIndices,
        const std::tuple<Axes...>& axes)
    {
      urEdge.at(0u) = std::get<0u>(axes).getBinUpperBound(localIndices.at(0u));
    }

    template <class... Axes>
    static void
    getUpperRightBinIndices(std::array<size_t, sizeof...(Axes)>& localIndices,
                            const std::tuple<Axes...>& axes)
    {
      size_t thisAxisNBins = std::get<0u>(axes).getNBins();
      localIndices.at(0u)  = std::min(thisAxisNBins + 1, ++localIndices.at(0u));
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

    /// @brief retrieve lower-left bin edge from set of local bin indices
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    /// @param  [in] localIndices local bin indices along each axis
    /// @param  [in] axes         actual axis objects spanning the grid
    /// @return generalized lower-left bin edge
    ///
    /// @pre @c localIndices must only contain valid bin indices (i.e. excluding
    ///      under-/overflow bins).
    template <class... Axes>
    static std::array<double, sizeof...(Axes)>
    getLowerLeftBinEdge(const std::array<size_t, sizeof...(Axes)>& localIndices,
                        const std::tuple<Axes...>& axes)
    {
      std::array<double, sizeof...(Axes)> llEdge;
      constexpr size_t MAX = sizeof...(Axes)-1;
      grid_bins_helper_impl<MAX>::getLowerLeftBinEdge(
          llEdge, localIndices, axes);

      return llEdge;
    }

    /// @brief get local bin indices for lower-left neighboring bin
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    /// @param  [in] localIndices local bin indices along each axis
    /// @param  [in] axes         actual axis objects spanning the grid
    /// @return array with local bin indices of lower-left neighbor bin
    ///
    /// @pre @c localIndices must only contain valid bin indices (i.e. excluding
    ///      under-/overflow bins).
    ///
    /// This function returns the local bin indices for the generalized
    /// lower-left neighbor which simply means that all local bin indices are
    /// decremented by one.
    template <class... Axes>
    static std::array<size_t, sizeof...(Axes)>
    getLowerLeftBinIndices(
        const std::array<size_t, sizeof...(Axes)>& localIndices,
        const std::tuple<Axes...>&)
    {
      auto llIndices = localIndices;
      for (size_t i = 0; i < sizeof...(Axes); ++i) --llIndices.at(i);

      return llIndices;
    }

    /// @brief calculate total number of bins in a grid defined by a set of
    /// axes
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

    /// @brief retrieve upper-right bin edge from set of local bin indices
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    /// @param  [in] localIndices local bin indices along each axis
    /// @param  [in] axes         actual axis objects spanning the grid
    /// @return generalized upper-right bin edge
    ///
    /// @pre @c localIndices must only contain valid bin indices (i.e. excluding
    ///      under-/overflow bins).
    template <class... Axes>
    static std::array<double, sizeof...(Axes)>
    getUpperRightBinEdge(
        const std::array<size_t, sizeof...(Axes)>& localIndices,
        const std::tuple<Axes...>& axes)
    {
      std::array<double, sizeof...(Axes)> urEdge;
      constexpr size_t MAX = sizeof...(Axes)-1;
      grid_bins_helper_impl<MAX>::getUpperRightBinEdge(
          urEdge, localIndices, axes);

      return urEdge;
    }

    /// @brief get local bin indices for upper-right neighboring bin
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    /// @param  [in] localIndices local bin indices along each axis
    /// @param  [in] axes         actual axis objects spanning the grid
    /// @return array with local bin indices of upper-right neighbor bin
    ///
    /// @pre @c localIndices must only contain valid bin indices (i.e. excluding
    ///      under-/overflow bins).
    ///
    /// This function returns the local bin indices for the generalized
    /// upper-right neighbor which simply means that all local bin indices are
    /// incremented by one.
    template <class... Axes>
    static std::array<size_t, sizeof...(Axes)>
    getUpperRightBinIndices(
        const std::array<size_t, sizeof...(Axes)>& localIndices,
        const std::tuple<Axes...>& axes)
    {
      constexpr size_t MAX       = sizeof...(Axes)-1;
      auto             urIndices = localIndices;
      grid_bins_helper_impl<MAX>::getUpperRightBinIndices(urIndices, axes);

      return urIndices;
    }
  };

}  // namespace detail

}  // namespace Acts
