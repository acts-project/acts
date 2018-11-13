// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <tuple>
#include <utility>
#include "Acts/Utilities/IAxis.hpp"

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
    static void
    getLowerLeftBinIndices(std::array<size_t, sizeof...(Axes)>& localIndices,
                           const std::tuple<Axes...>& axes)
    {
      localIndices.at(N) = std::get<N>(axes).wrapBin(localIndices.at(N) - 1);
      grid_helper_impl<N - 1>::getLowerLeftBinIndices(localIndices, axes);
    }

    template <class... Axes>
    static void
    getNBins(const std::tuple<Axes...>& axes,
             std::array<size_t, sizeof...(Axes)>& nBinsArray)
    {
      // by convention getNBins does not include under-/overflow bins
      nBinsArray[N] = std::get<N>(axes).getNBins();
      grid_helper_impl<N - 1>::getNBins(axes, nBinsArray);
    }

    template <class... Axes>
    static void
    getAxes(const std::tuple<Axes...>& axes,
            std::array<const IAxis*, sizeof...(Axes)>& axesArr)
    {
      axesArr[N] = static_cast<const IAxis*>(&std::get<N>(axes));
      grid_helper_impl<N - 1>::getAxes(axes, axesArr);
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
      localIndices.at(N) = std::get<N>(axes).wrapBin(localIndices.at(N) + 1);
      grid_helper_impl<N - 1>::getUpperRightBinIndices(localIndices, axes);
    }

    template <class... Axes>
    static void
    getMin(const std::tuple<Axes...>& axes,
           std::array<double, sizeof...(Axes)>& minArray)
    {
      minArray[N] = std::get<N>(axes).getMin();
      grid_helper_impl<N - 1>::getMin(axes, minArray);
    }

    template <class... Axes>
    static void
    getMax(const std::tuple<Axes...>& axes,
           std::array<double, sizeof...(Axes)>& maxArray)
    {
      maxArray[N] = std::get<N>(axes).getMax();
      grid_helper_impl<N - 1>::getMax(axes, maxArray);
    }

    template <class Point, class... Axes>
    static bool
    isInside(const Point& position, const std::tuple<Axes...>& axes)
    {
      bool insideThisAxis = std::get<N>(axes).isInside(position[N]);
      return insideThisAxis
          && grid_helper_impl<N - 1>::isInside(position, axes);
    }

    template <class... Axes>
    static void
    neighborHoodIndices(
        const std::array<size_t, sizeof...(Axes)>& localIndices,
        std::pair<size_t, size_t>                  sizes,
        const std::tuple<Axes...>& axes,
        std::array<std::set<size_t>, sizeof...(Axes)>& neighborIndices)
    {

      // ask n-th axis
      size_t           locIdx = localIndices.at(N);
      std::set<size_t> locNeighbors
          = std::get<N>(axes).neighborHoodIndices(locIdx, sizes);
      neighborIndices.at(N) = locNeighbors;

      grid_helper_impl<N - 1>::neighborHoodIndices(
          localIndices, sizes, axes, neighborIndices);
    }

    template <class... Axes>
    static void
    combineNeighborHoodIndices(
        std::array<size_t, sizeof...(Axes)>&           idx,
        std::array<std::set<size_t>, sizeof...(Axes)>& neighborIndices,
        std::set<size_t>&          combinations,
        const std::tuple<Axes...>& axes)
    {
      // iterate over this axis' neighbors
      std::set<size_t>& neighbors = neighborIndices.at(N);
      for (const auto& i : neighbors) {
        idx.at(N) = i;
        // vary other axes recursively
        grid_helper_impl<N - 1>::combineNeighborHoodIndices(
            idx, neighborIndices, combinations, axes);
      }
    }

    template <class... Axes>
    static void
    exteriorBinIndices(std::array<size_t, sizeof...(Axes)>& idx,
                       std::array<bool, sizeof...(Axes)>    isExterior,
                       std::set<size_t>&          combinations,
                       const std::tuple<Axes...>& axes)
    {
      // iterate over this axis' bins, remembering which bins are exterior
      for (size_t i = 0; i < std::get<N>(axes).getNBins() + 2; ++i) {
        idx.at(N)        = i;
        isExterior.at(N) = (i == 0) || (i == std::get<N>(axes).getNBins() + 1);
        // vary other axes recursively
        grid_helper_impl<N - 1>::exteriorBinIndices(
            idx, isExterior, combinations, axes);
      }
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
                 const std::tuple<Axes...>& /*axes*/,
                 size_t& bin,
                 size_t& area)
    {
      bin += area * localBins.at(0u);
    }

    template <class... Axes>
    static void
    getLocalBinIndices(size_t& bin,
                       const std::tuple<Axes...>& /*axes*/,
                       size_t& area,
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
    static void
    getLowerLeftBinIndices(std::array<size_t, sizeof...(Axes)>& localIndices,
                           const std::tuple<Axes...>& axes)
    {
      localIndices.at(0u) = std::get<0u>(axes).wrapBin(localIndices.at(0u) - 1);
    }

    template <class... Axes>
    static void
    getNBins(const std::tuple<Axes...>& axes,
             std::array<size_t, sizeof...(Axes)>& nBinsArray)
    {
      // by convention getNBins does not include under-/overflow bins
      nBinsArray[0u] = std::get<0u>(axes).getNBins();
    }

    template <class... Axes>
    static void
    getAxes(const std::tuple<Axes...>& axes,
            std::array<const IAxis*, sizeof...(Axes)>& axesArr)
    {
      axesArr[0u] = static_cast<const IAxis*>(&std::get<0u>(axes));
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
      localIndices.at(0u) = std::get<0u>(axes).wrapBin(localIndices.at(0u) + 1);
    }

    template <class... Axes>
    static void
    getMin(const std::tuple<Axes...>& axes,
           std::array<double, sizeof...(Axes)>& minArray)
    {
      minArray[0u] = std::get<0u>(axes).getMin();
    }

    template <class... Axes>
    static void
    getMax(const std::tuple<Axes...>& axes,
           std::array<double, sizeof...(Axes)>& maxArray)
    {
      maxArray[0u] = std::get<0u>(axes).getMax();
    }

    template <class Point, class... Axes>
    static bool
    isInside(const Point& position, const std::tuple<Axes...>& axes)
    {
      return std::get<0u>(axes).isInside(position[0u]);
    }

    template <class... Axes>
    static void
    neighborHoodIndices(
        const std::array<size_t, sizeof...(Axes)>& localIndices,
        std::pair<size_t, size_t>                  sizes,
        const std::tuple<Axes...>& axes,
        std::array<std::set<size_t>, sizeof...(Axes)>& neighborIndices)
    {

      // ask 0-th axis
      size_t           locIdx = localIndices.at(0u);
      std::set<size_t> locNeighbors
          = std::get<0u>(axes).neighborHoodIndices(locIdx, sizes);
      neighborIndices.at(0u) = locNeighbors;
    }

    template <class... Axes>
    static void
    combineNeighborHoodIndices(
        std::array<size_t, sizeof...(Axes)>&           idx,
        std::array<std::set<size_t>, sizeof...(Axes)>& neighborIndices,
        std::set<size_t>&          combinations,
        const std::tuple<Axes...>& axes)
    {
      // iterate over this axis' neighbors
      std::set<size_t>& neighbors = neighborIndices.at(0u);
      for (const auto& i : neighbors) {
        idx.at(0u) = i;
        // at this point, combinations are complete, so fill
        size_t bin = 0, area = 1;
        grid_helper_impl<sizeof...(Axes) - 1>::getGlobalBin(
            idx, axes, bin, area);
        combinations.insert(bin);
      }
    }

    template <class... Axes>
    static void
    exteriorBinIndices(std::array<size_t, sizeof...(Axes)>& idx,
                       std::array<bool, sizeof...(Axes)>    isExterior,
                       std::set<size_t>&          combinations,
                       const std::tuple<Axes...>& axes)
    {
      // For each exterior bin on this axis, we will do this
      auto recordExteriorBin = [&](size_t i) {
        idx.at(0u) = i;
        // at this point, combinations are complete: save the global bin
        size_t bin = 0, area = 1;
        grid_helper_impl<sizeof...(Axes) - 1>::getGlobalBin(
            idx, axes, bin, area);
        combinations.insert(bin);
      };

      // The first and last bins on this axis are exterior by definition
      for (size_t i :
           {static_cast<size_t>(0), std::get<0u>(axes).getNBins() + 1}) {
        recordExteriorBin(i);
      }

      // If no other axis is on an exterior index, stop here
      bool otherAxisExterior = false;
      for (size_t N = 1; N < sizeof...(Axes); ++N) {
        otherAxisExterior = otherAxisExterior | isExterior[N];
      }
      if (!otherAxisExterior) {
        return;
      }

      // Otherwise, we're on a grid border: iterate over all the other indices
      for (size_t i = 1; i <= std::get<0u>(axes).getNBins(); ++i) {
        recordExteriorBin(i);
      }
    }
  };
  /// @endcond

  /// @brief helper functions for grid-related operations
  struct grid_helper
  {
    /// @brief get the global indices for closest points on grid
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    /// @param  [in] bin  global bin index for bin of interest
    /// @param  [in] axes actual axis objects spanning the grid
    /// @return set of global bin indices for bins whose lower-left corners are
    ///         the closest points on the grid to every point in the given bin
    ///
    /// @note @c bin must be a valid bin index (excluding under-/overflow bins
    ///       along any axis).
    template <class... Axes>
    static std::set<size_t>
    closestPointsIndices(size_t bin, const std::tuple<Axes...>& axes)
    {
      std::array<size_t, sizeof...(Axes)> localIndices
          = getLocalBinIndices(bin, axes);
      // get neighboring bins, but only increment.
      std::array<std::set<size_t>, sizeof...(Axes)> neighborIndices;
      std::set<size_t> comb
          = neighborHoodIndices(localIndices, std::make_pair(0, 1), axes);
      return comb;
    }

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
      constexpr size_t MAX = sizeof...(Axes) - 1;
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
      constexpr size_t MAX  = sizeof...(Axes) - 1;
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
      constexpr size_t MAX  = sizeof...(Axes) - 1;
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
      constexpr size_t MAX  = sizeof...(Axes) - 1;
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
      constexpr size_t MAX = sizeof...(Axes) - 1;
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
      constexpr size_t MAX       = sizeof...(Axes) - 1;
      auto             llIndices = localIndices;
      grid_helper_impl<MAX>::getLowerLeftBinIndices(llIndices, axes);

      return llIndices;
    }

    /// @brief calculate number of bins in a grid defined by a set of
    /// axes for each axis
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    /// @param  [in] axes actual axis objects spanning the grid
    /// @return array of number of bins for each axis of the grid
    ///
    /// @note This does not include under-/overflow bins along each axis.
    template <class... Axes>
    static std::array<size_t, sizeof...(Axes)>
    getNBins(const std::tuple<Axes...>& axes)
    {
      std::array<size_t, sizeof...(Axes)> nBinsArray;
      grid_helper_impl<sizeof...(Axes) - 1>::getNBins(axes, nBinsArray);
      return nBinsArray;
    }

    /// @brief return an array with copies of the axes, converted
    /// to type AnyAxis
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    /// @param [in] axes actual axis objects spanning the grid
    /// @return array with copies of the axis
    template <class... Axes>
    static std::array<const IAxis*, sizeof...(Axes)>
    getAxes(const std::tuple<Axes...>& axes)
    {
      std::array<const IAxis*, sizeof...(Axes)> arr;
      grid_helper_impl<sizeof...(Axes) - 1>::getAxes(axes, arr);
      return arr;
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
      constexpr size_t MAX = sizeof...(Axes) - 1;
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
      constexpr size_t MAX       = sizeof...(Axes) - 1;
      auto             urIndices = localIndices;
      grid_helper_impl<MAX>::getUpperRightBinIndices(urIndices, axes);

      return urIndices;
    }

    /// @brief get the minimum value of all axes of one grid
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    /// @param  [in] axes actual axis objects spanning the grid
    /// @return array returning the minima of all given axes
    template <class... Axes>
    static std::array<double, sizeof...(Axes)>
    getMin(const std::tuple<Axes...>& axes)
    {
      std::array<double, sizeof...(Axes)> minArray;
      grid_helper_impl<sizeof...(Axes) - 1>::getMin(axes, minArray);
      return minArray;
    }

    /// @brief get the maximum value of all axes of one grid
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    /// @param  [in] axes actual axis objects spanning the grid
    /// @return array returning the maxima of all given axes
    template <class... Axes>
    static std::array<double, sizeof...(Axes)>
    getMax(const std::tuple<Axes...>& axes)
    {
      std::array<double, sizeof...(Axes)> maxArray;
      grid_helper_impl<sizeof...(Axes) - 1>::getMax(axes, maxArray);
      return maxArray;
    }

    /// @brief get global bin indices for bins in specified neighborhood
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    /// @param  [in] localIndices local bin indices along each axis
    /// @param  [in] size         size of neighborhood determining how many
    ///                           adjacent bins along each axis are considered
    /// @param  [in] axes         actual axis objects spanning the grid
    /// @return set of global bin indices for all bins in neighborhood
    ///
    /// @note Over-/underflow bins are included in the neighborhood.
    /// @note The @c size parameter sets the range by how many units each local
    ///       bin index is allowed to be varied. All local bin indices are
    ///       varied independently, that is diagonal neighbors are included.
    ///       Ignoring the truncation of the neighborhood size reaching beyond
    ///       over-/underflow bins, the neighborhood is of size \f$2 \times
    ///       \text{size}+1\f$ along each dimension.
    /// @note The concrete bins which are returned depend on the WrappingTypes
    ///       of the contained axes
    ///
    /// @todo At the moment, this function relies on local-to-global bin index
    ///       conversions for each neighboring bin. In principle, all
    ///       information should be available to work directly on global bin
    ///       indices. The problematic part is the check when going beyond
    ///       under-/overflow bins.
    template <class... Axes>
    static std::set<size_t>
    neighborHoodIndices(const std::array<size_t, sizeof...(Axes)>& localIndices,
                        std::pair<size_t, size_t>                  sizes,
                        const std::tuple<Axes...>& axes)
    {
      constexpr size_t MAX = sizeof...(Axes) - 1;

      // length N array which contains local neighbors based on size par
      std::array<std::set<size_t>, sizeof...(Axes)> neighborIndices;
      // get local bin indices for neighboring bins
      grid_helper_impl<MAX>::neighborHoodIndices(
          localIndices, sizes, axes, neighborIndices);

      std::array<size_t, sizeof...(Axes)> idx;
      std::set<size_t> combinations;
      grid_helper_impl<MAX>::combineNeighborHoodIndices(
          idx, neighborIndices, combinations, axes);

      return combinations;
    }

    template <class... Axes>
    static std::set<size_t>
    neighborHoodIndices(const std::array<size_t, sizeof...(Axes)>& localIndices,
                        size_t                     size,
                        const std::tuple<Axes...>& axes)
    {
      return neighborHoodIndices(
          localIndices, std::make_pair(size, size), axes);
    }

    /// @brief get bin indices of all overflow and underflow bins
    ///
    /// @tparam Axes parameter pack of axis types defining the grid
    /// @param  [in] axes         actual axis objects spanning the grid
    /// @return set of global bin indices for all over- and underflow bins
    template <class... Axes>
    static std::set<size_t>
    exteriorBinIndices(const std::tuple<Axes...>& axes)
    {
      constexpr size_t MAX = sizeof...(Axes) - 1;

      std::array<size_t, sizeof...(Axes)> idx;
      std::array<bool, sizeof...(Axes)>   isExterior;
      std::set<size_t> combinations;
      grid_helper_impl<MAX>::exteriorBinIndices(
          idx, isExterior, combinations, axes);

      return combinations;
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
      constexpr size_t MAX = sizeof...(Axes) - 1;
      return grid_helper_impl<MAX>::isInside(position, axes);
    }
  };

}  // namespace detail

}  // namespace Acts
