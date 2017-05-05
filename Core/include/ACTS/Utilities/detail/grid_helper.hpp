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
  struct grid_helper_impl;

  template <size_t N>
  struct grid_helper_impl
  {
    template <class... Axes>
    static void
    getBinCenter(std::array<double, sizeof...(Axes)>&       center,
                 const std::array<size_t, sizeof...(Axes)>& localIndices,
                 const std::tuple<Axes...>& axes)
    {
      center.at(N) = std::get<N>(axes).getBinCenter(localIndices.at(N));
      grid_helper_impl<N - 1>::getBinCenter(center, localIndices, axes);
    }

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
      grid_helper_impl<N - 1>::getGlobalBin(point, axes, bin, area);
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
      grid_helper_impl<N - 1>::getGlobalBin(localBins, axes, bin, area);
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
      grid_helper_impl<N - 1>::getLocalBinIndices(bin, axes, new_area, indices);
      indices.at(N) = bin / area;
      bin %= area;
    }

    template <class... Axes>
    static void
    getLowerLeftBinEdge(std::array<double, sizeof...(Axes)>&       llEdge,
                        const std::array<size_t, sizeof...(Axes)>& localIndices,
                        const std::tuple<Axes...>& axes)
    {
      llEdge.at(N) = std::get<N>(axes).getBinLowerBound(localIndices.at(N));
      grid_helper_impl<N - 1>::getLowerLeftBinEdge(llEdge, localIndices, axes);
    }

    template <class... Axes>
    static size_t
    getNBins(const std::tuple<Axes...>& axes)
    {
      // by convention getNBins does not include under-/overflow bins
      size_t thisAxisNBins = std::get<N>(axes).getNBins() + 2;
      return thisAxisNBins * grid_helper_impl<N - 1>::getNBins(axes);
    }

    template <class... Axes>
    static void
    getUpperRightBinEdge(
        std::array<double, sizeof...(Axes)>&       urEdge,
        const std::array<size_t, sizeof...(Axes)>& localIndices,
        const std::tuple<Axes...>& axes)
    {
      urEdge.at(N) = std::get<N>(axes).getBinUpperBound(localIndices.at(N));
      grid_helper_impl<N - 1>::getUpperRightBinEdge(urEdge, localIndices, axes);
    }

    template <class... Axes>
    static void
    getUpperRightBinIndices(std::array<size_t, sizeof...(Axes)>& localIndices,
                            const std::tuple<Axes...>& axes)
    {
      size_t thisAxisNBins = std::get<N>(axes).getNBins();
      localIndices.at(N)   = std::min(thisAxisNBins + 1, ++localIndices.at(N));
      grid_helper_impl<N - 1>::getUpperRightBinIndices(localIndices, axes);
    }

    template <class Point, class... Axes>
    static bool
    isInside(const Point& position, const std::tuple<Axes...>& axes)
    {
      bool insideThisAxis = std::get<N>(axes).isInside(position[N]);
      return insideThisAxis
          && grid_helper_impl<N - 1>::isInside(position, axes);
    }
  };

  template <>
  struct grid_helper_impl<0u>
  {
    template <class... Axes>
    static void
    getBinCenter(std::array<double, sizeof...(Axes)>&       center,
                 const std::array<size_t, sizeof...(Axes)>& localIndices,
                 const std::tuple<Axes...>& axes)
    {
      center.at(0u) = std::get<0u>(axes).getBinCenter(localIndices.at(0u));
    }

    template <class Point, class... Axes>
    static void
    getGlobalBin(const Point&               point,
                 const std::tuple<Axes...>& axes,
                 size_t&                    bin,
                 size_t&                    area)
    {
      const auto& thisAxis = std::get<0u>(axes);
      bin += area * thisAxis.getBin(point[0u]);
    }

    template <class... Axes>
    static void
    getGlobalBin(const std::array<size_t, sizeof...(Axes)>& localBins,
                 const std::tuple<Axes...>& axes,
                 size_t&                    bin,
                 size_t&                    area)
    {
      bin += area * localBins.at(0u);
    }

    template <class... Axes>
    static void
    getLocalBinIndices(size_t&                    bin,
                       const std::tuple<Axes...>& axes,
                       size_t&                    area,
                       std::array<size_t, sizeof...(Axes)>& indices)
    {
      // make sure to account for under-/overflow bins
      indices.at(0u) = bin / area;
      bin %= area;
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

    template <class Point, class... Axes>
    static bool
    isInside(const Point& position, const std::tuple<Axes...>& axes)
    {
      return std::get<0u>(axes).isInside(position[0u]);
    }
  };
  /// @endcond

  /// @brief helper functions for grid-related operations
  struct grid_helper
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
      grid_helper_impl<MAX>::getBinCenter(center, localIndices, axes);

      return center;
    }

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
    /// @pre The given @c Point type must represent a point in d (or higher)
    ///      dimensions where d is the number of axis objects in the tuple.
    /// @note This could be a under-/overflow bin along one or more axes.
    template <class Point, class... Axes>
    static size_t
    getGlobalBin(const Point& point, const std::tuple<Axes...>& axes)
    {
      constexpr size_t MAX  = sizeof...(Axes)-1;
      size_t           area = 1;
      size_t           bin  = 0;

      grid_helper_impl<MAX>::getGlobalBin(point, axes, bin, area);

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

      grid_helper_impl<MAX>::getGlobalBin(localBins, axes, bin, area);

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
    /// @note Local bin indices can contain under-/overflow bins along any axis.
    template <class... Axes>
    static std::array<size_t, sizeof...(Axes)>
    getLocalBinIndices(size_t bin, const std::tuple<Axes...>& axes)
    {
      constexpr size_t MAX  = sizeof...(Axes)-1;
      size_t           area = 1;
      std::array<size_t, sizeof...(Axes)> indices;

      grid_helper_impl<MAX>::getLocalBinIndices(bin, axes, area, indices);

      return indices;
    }

    /// @brief retrieve lower-left bin edge from set of local bin indices
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    /// @param  [in] localIndices local bin indices along each axis
    /// @param  [in] axes         actual axis objects spanning the grid
    /// @return generalized lower-left bin edge
    ///
    /// @pre @c localIndices must only contain valid bin indices (excluding
    ///      underflow bins).
    template <class... Axes>
    static std::array<double, sizeof...(Axes)>
    getLowerLeftBinEdge(const std::array<size_t, sizeof...(Axes)>& localIndices,
                        const std::tuple<Axes...>& axes)
    {
      std::array<double, sizeof...(Axes)> llEdge;
      constexpr size_t MAX = sizeof...(Axes)-1;
      grid_helper_impl<MAX>::getLowerLeftBinEdge(llEdge, localIndices, axes);

      return llEdge;
    }

    /// @brief get local bin indices for lower-left neighboring bin
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    /// @param  [in] localIndices local bin indices along each axis
    /// @param  [in] axes         actual axis objects spanning the grid
    /// @return array with local bin indices of lower-left neighbor bin
    ///
    /// @pre @c localIndices must only contain valid bin indices (excluding
    ///      underflow bins).
    ///
    /// This function returns the local bin indices for the generalized
    /// lower-left neighbor which simply means that all local bin indices are
    /// decremented by one.
    template <class... Axes>
    static std::array<size_t, sizeof...(Axes)>
    getLowerLeftBinIndices(
        const std::array<size_t, sizeof...(Axes)>& localIndices,
        const std::tuple<Axes...>& axes)
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
      return grid_helper_impl<MAX>::getNBins(axes);
    }

    /// @brief retrieve upper-right bin edge from set of local bin indices
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    /// @param  [in] localIndices local bin indices along each axis
    /// @param  [in] axes         actual axis objects spanning the grid
    /// @return generalized upper-right bin edge
    ///
    /// @pre @c localIndices must only contain valid bin indices (excluding
    ///      overflow bins).
    template <class... Axes>
    static std::array<double, sizeof...(Axes)>
    getUpperRightBinEdge(
        const std::array<size_t, sizeof...(Axes)>& localIndices,
        const std::tuple<Axes...>& axes)
    {
      std::array<double, sizeof...(Axes)> urEdge;
      constexpr size_t MAX = sizeof...(Axes)-1;
      grid_helper_impl<MAX>::getUpperRightBinEdge(urEdge, localIndices, axes);

      return urEdge;
    }

    /// @brief get local bin indices for upper-right neighboring bin
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    /// @param  [in] localIndices local bin indices along each axis
    /// @param  [in] axes         actual axis objects spanning the grid
    /// @return array with local bin indices of upper-right neighbor bin
    ///
    /// @pre @c localIndices must only contain valid bin indices (excluding
    ///      overflow bins).
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
      grid_helper_impl<MAX>::getUpperRightBinIndices(urIndices, axes);

      return urIndices;
    }

    /// @brief check whether given point is inside axes limits
    ///
    /// @tparam Point any type with point semantics supporting component access
    ///               through @c operator[]
    /// @tparam Axes parameter pack of axis types defining the grid
    ///
    /// @param  [in] position point to look up in the grid
    /// @param  [in] axes     actual axis objects spanning the grid
    /// @return @c true if \f$\text{xmin_i} \le x_i < \text{xmax}_i \forall i=0,
    ///         \dots, d-1\f$, otherwise @c false
    ///
    /// @pre The given @c Point type must represent a point in d (or higher)
    ///      dimensions where d is the number of axis objects in the tuple.
    template <class Point, class... Axes>
    static bool
    isInside(const Point& position, const std::tuple<Axes...>& axes)
    {
      constexpr size_t MAX = sizeof...(Axes)-1;
      return grid_helper_impl<MAX>::isInside(position, axes);
    }
  };

}  // namespace detail

}  // namespace Acts
