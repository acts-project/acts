// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/NeighborHoodIndices.hpp"

#include <array>
#include <set>
#include <tuple>
#include <utility>

#include <boost/container/small_vector.hpp>

namespace Acts::detail {

/// This object can be iterated to produce the (ordered) set of flat indices
/// associated with a neighborhood around a certain point on a grid.
///
/// The goal is to emulate the effect of enumerating the flat indices into
/// an std::set (or into an std::vector that gets subsequently sorted), without
/// paying the price of dynamic memory allocation in hot magnetic field
/// interpolation code.
template <std::size_t DIM>
class FlatNeighborHoodIndices {
 public:
  /// You can get the neighbor multi indices from
  /// MultiAxisHelperImpl<DIM>::neighborHoodIndices and the number of bins in
  /// each direction from MultiAxisHelperImpl<DIM>::getNBins.
  FlatNeighborHoodIndices(std::array<NeighborHoodIndices, DIM>& neighborIndices,
                          const std::array<std::size_t, DIM>& nBinsArray)
      : m_multiIndices(neighborIndices) {
    if (DIM == 1) {
      return;
    }
    std::size_t flatStride = 1;
    for (long i = DIM - 2; i >= 0; --i) {
      flatStride *= (nBinsArray[i + 1] + 2);
      m_flatStride[i] = flatStride;
    }
  }

  class iterator {
   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using pointer = std::size_t*;
    using reference = std::size_t&;

    iterator() = default;

    iterator(const FlatNeighborHoodIndices& parent,
             std::array<NeighborHoodIndices::iterator, DIM>&& multiIndicesIter)
        : m_multiIndicesIter(std::move(multiIndicesIter)), m_parent(&parent) {}

    std::size_t operator*() const {
      std::size_t flatIndex = *m_multiIndicesIter[DIM - 1];
      if (DIM == 1) {
        return flatIndex;
      }
      for (std::size_t i = 0; i < DIM - 1; ++i) {
        flatIndex += m_parent->m_flatStride[i] * (*m_multiIndicesIter[i]);
      }
      return flatIndex;
    }

    iterator& operator++() {
      const auto& multiIndices = m_parent->m_multiIndices;

      // Go to the next flat index via a lexicographic increment:
      // - Start by incrementing the last multi index
      // - If it reaches the end, reset it and increment the previous one...
      for (long i = DIM - 1; i > 0; --i) {
        ++m_multiIndicesIter[i];
        if (m_multiIndicesIter[i] != multiIndices[i].end()) {
          return *this;
        }
        m_multiIndicesIter[i] = multiIndices[i].begin();
      }

      // The first index should stay at the end value when it reaches it, so
      // that we know when we've reached the end of iteration.
      ++m_multiIndicesIter[0];
      return *this;
    }

    iterator operator++(int) {
      iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool isEqual(const iterator& b) const {
      if (b.m_parent == nullptr) {
        return m_multiIndicesIter[0] == m_parent->m_multiIndices[0].end();
      } else {
        return m_multiIndicesIter == b.m_multiIndicesIter;
      }
    }

    friend bool operator==(const iterator& a, const iterator& b) {
      return a.isEqual(b);
    }

   private:
    std::array<NeighborHoodIndices::iterator, DIM> m_multiIndicesIter;
    const FlatNeighborHoodIndices* m_parent = nullptr;
  };

  iterator begin() const {
    std::array<NeighborHoodIndices::iterator, DIM> multiIndicesIter{};
    for (std::size_t i = 0; i < DIM; ++i) {
      multiIndicesIter[i] = m_multiIndices[i].begin();
    }
    return iterator(*this, std::move(multiIndicesIter));
  }

  iterator end() const { return iterator(); }

  /// Number of indices that will be produced if this sequence is iterated
  std::size_t size() const {
    std::size_t result = m_multiIndices[0].size();
    for (std::size_t i = 1; i < DIM; ++i) {
      result *= m_multiIndices[i].size();
    }
    return result;
  }

  /// Collect the sequence of indices into a vector
  auto collect() const {
    boost::container::small_vector<std::size_t, ipow(3, DIM)> result;
    result.reserve(this->size());
    for (std::size_t idx : *this) {
      result.push_back(idx);
    }
    return result;
  }

  std::vector<std::size_t> collectVector() const {
    auto result = collect();
    return {result.begin(), result.end()};
  }

 private:
  std::array<NeighborHoodIndices, DIM> m_multiIndices{};
  std::array<std::size_t, DIM - 1> m_flatStride{};
};

/// @cond
/// helper struct to calculate number of bins inside a grid
///
/// @tparam N number of axes to consider
template <std::size_t N>
struct MultiAxisHelperImpl;

template <std::size_t N>
struct MultiAxisHelperImpl {
  template <class... Axes>
  static void getBinCenter(
      std::array<double, sizeof...(Axes)>& center,
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes) {
    center.at(N) = std::get<N>(axes).getBinCenter(multiIndex.at(N));
    MultiAxisHelperImpl<N - 1>::getBinCenter(center, multiIndex, axes);
  }

  template <class... Axes>
  static void getFlatIndexFromMultiIndex(
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes, std::size_t& bin, std::size_t area) {
    const auto& thisAxis = std::get<N>(axes);
    bin += area * multiIndex.at(N);
    // make sure to account for under-/overflow bins
    area *= (thisAxis.getNBins() + 2);
    MultiAxisHelperImpl<N - 1>::getFlatIndexFromMultiIndex(multiIndex, axes,
                                                           bin, area);
  }

  template <class Point, class... Axes>
  static void getMultiIndexFromPoint(
      const Point& point, const std::tuple<Axes...>& axes,
      std::array<std::size_t, sizeof...(Axes)>& multiIndex) {
    const auto& thisAxis = std::get<N>(axes);
    multiIndex.at(N) = static_cast<std::size_t>(thisAxis.getBin(point[N]));
    MultiAxisHelperImpl<N - 1>::getMultiIndexFromPoint(point, axes, multiIndex);
  }

  template <class... Axes>
  static void getMultiIndexFromFlatIndex(
      std::size_t& bin, const std::tuple<Axes...>& axes, const std::size_t area,
      std::array<std::size_t, sizeof...(Axes)>& multiIndex) {
    const auto& thisAxis = std::get<N>(axes);
    // make sure to account for under-/overflow bins
    std::size_t new_area = area * (thisAxis.getNBins() + 2);
    MultiAxisHelperImpl<N - 1>::getMultiIndexFromFlatIndex(bin, axes, new_area,
                                                           multiIndex);
    multiIndex.at(N) = bin / area;
    bin %= area;
  }

  template <class... Axes>
  static void getLowerLeftBinCorner(
      std::array<double, sizeof...(Axes)>& llEdge,
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes) {
    llEdge.at(N) = std::get<N>(axes).getBinLowerBound(multiIndex.at(N));
    MultiAxisHelperImpl<N - 1>::getLowerLeftBinCorner(llEdge, multiIndex, axes);
  }

  template <class... Axes>
  static void getLowerLeftBinIndices(
      std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes) {
    multiIndex.at(N) = std::get<N>(axes).wrapBin(multiIndex.at(N) - 1);
    MultiAxisHelperImpl<N - 1>::getLowerLeftBinIndices(multiIndex, axes);
  }

  template <class... Axes>
  static void getNBins(const std::tuple<Axes...>& axes,
                       std::array<std::size_t, sizeof...(Axes)>& nBinsArray) {
    // by convention getNBins does not include under-/overflow bins
    nBinsArray[N] = std::get<N>(axes).getNBins();
    MultiAxisHelperImpl<N - 1>::getNBins(axes, nBinsArray);
  }

  template <class... Axes>
  static void getAxes(const std::tuple<Axes...>& axes,
                      std::array<const IAxis*, sizeof...(Axes)>& axesArr) {
    axesArr[N] = static_cast<const IAxis*>(&std::get<N>(axes));
    MultiAxisHelperImpl<N - 1>::getAxes(axes, axesArr);
  }

  template <class... Axes>
  static void getUpperRightBinCorner(
      std::array<double, sizeof...(Axes)>& urEdge,
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes) {
    urEdge.at(N) = std::get<N>(axes).getBinUpperBound(multiIndex.at(N));
    MultiAxisHelperImpl<N - 1>::getUpperRightBinCorner(urEdge, multiIndex,
                                                       axes);
  }

  template <class... Axes>
  static void getUpperRightBinIndices(
      std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes) {
    multiIndex.at(N) = std::get<N>(axes).wrapBin(multiIndex.at(N) + 1);
    MultiAxisHelperImpl<N - 1>::getUpperRightBinIndices(multiIndex, axes);
  }

  template <class... Axes>
  static void getMin(const std::tuple<Axes...>& axes,
                     std::array<double, sizeof...(Axes)>& minArray) {
    minArray[N] = std::get<N>(axes).getMin();
    MultiAxisHelperImpl<N - 1>::getMin(axes, minArray);
  }

  template <class... Axes>
  static void getMax(const std::tuple<Axes...>& axes,
                     std::array<double, sizeof...(Axes)>& maxArray) {
    maxArray[N] = std::get<N>(axes).getMax();
    MultiAxisHelperImpl<N - 1>::getMax(axes, maxArray);
  }

  template <class... Axes>
  static void getWidth(const std::tuple<Axes...>& axes,
                       std::array<double, sizeof...(Axes)>& widthArray) {
    widthArray[N] = std::get<N>(axes).getBinWidth();
    MultiAxisHelperImpl<N - 1>::getWidth(axes, widthArray);
  }

  template <class Point, class... Axes>
  static bool isInside(const Point& point, const std::tuple<Axes...>& axes) {
    bool insideThisAxis = std::get<N>(axes).isInside(point[N]);
    return insideThisAxis && MultiAxisHelperImpl<N - 1>::isInside(point, axes);
  }

  template <class... Axes>
  static void neighborHoodIndices(
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      std::pair<int, int> sizes, const std::tuple<Axes...>& axes,
      std::array<NeighborHoodIndices, sizeof...(Axes)>& neighborIndices) {
    // ask n-th axis
    std::size_t locIdx = multiIndex.at(N);
    NeighborHoodIndices locNeighbors =
        std::get<N>(axes).neighborHoodIndices(locIdx, sizes);
    neighborIndices.at(N) = locNeighbors;

    MultiAxisHelperImpl<N - 1>::neighborHoodIndices(multiIndex, sizes, axes,
                                                    neighborIndices);
  }

  template <class... Axes>
  static void neighborHoodIndices(
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      std::array<std::pair<int, int>, sizeof...(Axes)> sizes,
      const std::tuple<Axes...>& axes,
      std::array<NeighborHoodIndices, sizeof...(Axes)>& neighborIndices) {
    // ask n-th axis
    std::size_t locIdx = multiIndex.at(N);
    NeighborHoodIndices locNeighbors =
        std::get<N>(axes).neighborHoodIndices(locIdx, sizes.at(N));
    neighborIndices.at(N) = locNeighbors;

    MultiAxisHelperImpl<N - 1>::neighborHoodIndices(multiIndex, sizes, axes,
                                                    neighborIndices);
  }

  template <class... Axes>
  static void exteriorBinIndices(
      std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      std::array<bool, sizeof...(Axes)> isExterior,
      std::set<std::size_t>& combinations, const std::tuple<Axes...>& axes) {
    // iterate over this axis' bins, remembering which bins are exterior
    for (std::size_t i = 0; i < std::get<N>(axes).getNBins() + 2; ++i) {
      multiIndex.at(N) = i;
      isExterior.at(N) = (i == 0) || (i == std::get<N>(axes).getNBins() + 1);
      // vary other axes recursively
      MultiAxisHelperImpl<N - 1>::exteriorBinIndices(multiIndex, isExterior,
                                                     combinations, axes);
    }
  }
};

template <>
struct MultiAxisHelperImpl<0u> {
  template <class... Axes>
  static void getBinCenter(
      std::array<double, sizeof...(Axes)>& center,
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes) {
    center.at(0u) = std::get<0u>(axes).getBinCenter(multiIndex.at(0u));
  }

  template <class... Axes>
  static void getFlatIndexFromMultiIndex(
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& /*axes*/, std::size_t& bin,
      const std::size_t area) {
    bin += area * multiIndex.at(0u);
  }

  template <class Point, class... Axes>
  static void getMultiIndexFromPoint(
      const Point& point, const std::tuple<Axes...>& axes,
      std::array<std::size_t, sizeof...(Axes)>& multiIndex) {
    const auto& thisAxis = std::get<0u>(axes);
    multiIndex.at(0u) = thisAxis.getBin(point[0u]);
  }

  template <class... Axes>
  static void getMultiIndexFromFlatIndex(
      std::size_t& bin, const std::tuple<Axes...>& /*axes*/,
      const std::size_t area,
      std::array<std::size_t, sizeof...(Axes)>& multiIndex) {
    // make sure to account for under-/overflow bins
    multiIndex.at(0u) = bin / area;
    bin %= area;
  }

  template <class... Axes>
  static void getLowerLeftBinCorner(
      std::array<double, sizeof...(Axes)>& llEdge,
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes) {
    llEdge.at(0u) = std::get<0u>(axes).getBinLowerBound(multiIndex.at(0u));
  }

  template <class... Axes>
  static void getLowerLeftBinIndices(
      std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes) {
    multiIndex.at(0u) = std::get<0u>(axes).wrapBin(multiIndex.at(0u) - 1);
  }

  template <class... Axes>
  static void getNBins(const std::tuple<Axes...>& axes,
                       std::array<std::size_t, sizeof...(Axes)>& nBinsArray) {
    // by convention getNBins does not include under-/overflow bins
    nBinsArray[0u] = std::get<0u>(axes).getNBins();
  }

  template <class... Axes>
  static void getAxes(const std::tuple<Axes...>& axes,
                      std::array<const IAxis*, sizeof...(Axes)>& axesArr) {
    axesArr[0u] = static_cast<const IAxis*>(&std::get<0u>(axes));
  }

  template <class... Axes>
  static void getUpperRightBinCorner(
      std::array<double, sizeof...(Axes)>& urEdge,
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes) {
    urEdge.at(0u) = std::get<0u>(axes).getBinUpperBound(multiIndex.at(0u));
  }

  template <class... Axes>
  static void getUpperRightBinIndices(
      std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes) {
    multiIndex.at(0u) = std::get<0u>(axes).wrapBin(multiIndex.at(0u) + 1);
  }

  template <class... Axes>
  static void getMin(const std::tuple<Axes...>& axes,
                     std::array<double, sizeof...(Axes)>& minArray) {
    minArray[0u] = std::get<0u>(axes).getMin();
  }

  template <class... Axes>
  static void getMax(const std::tuple<Axes...>& axes,
                     std::array<double, sizeof...(Axes)>& maxArray) {
    maxArray[0u] = std::get<0u>(axes).getMax();
  }

  template <class... Axes>
  static void getWidth(const std::tuple<Axes...>& axes,
                       std::array<double, sizeof...(Axes)>& widthArray) {
    widthArray[0u] = std::get<0u>(axes).getBinWidth();
  }

  template <class Point, class... Axes>
  static bool isInside(const Point& point, const std::tuple<Axes...>& axes) {
    return std::get<0u>(axes).isInside(point[0u]);
  }

  template <class... Axes>
  static void neighborHoodIndices(
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      std::pair<int, int> sizes, const std::tuple<Axes...>& axes,
      std::array<NeighborHoodIndices, sizeof...(Axes)>& neighborIndices) {
    // ask 0-th axis
    std::size_t locIdx = multiIndex.at(0u);
    NeighborHoodIndices locNeighbors =
        std::get<0u>(axes).neighborHoodIndices(locIdx, sizes);
    neighborIndices.at(0u) = locNeighbors;
  }

  template <class... Axes>
  static void neighborHoodIndices(
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      std::array<std::pair<int, int>, sizeof...(Axes)> sizes,
      const std::tuple<Axes...>& axes,
      std::array<NeighborHoodIndices, sizeof...(Axes)>& neighborIndices) {
    // ask 0-th axis
    std::size_t locIdx = multiIndex.at(0u);
    NeighborHoodIndices locNeighbors =
        std::get<0u>(axes).neighborHoodIndices(locIdx, sizes.at(0u));
    neighborIndices.at(0u) = locNeighbors;
  }

  template <class... Axes>
  static void exteriorBinIndices(
      std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      std::array<bool, sizeof...(Axes)> isExterior,
      std::set<std::size_t>& combinations, const std::tuple<Axes...>& axes) {
    // For each exterior bin on this axis, we will do this
    auto recordExteriorBin = [&](std::size_t i) {
      multiIndex.at(0u) = i;
      // at this point, combinations are complete: save the global bin
      std::size_t bin = 0;
      const std::size_t area = 1;
      MultiAxisHelperImpl<sizeof...(Axes) - 1>::getFlatIndexFromMultiIndex(
          multiIndex, axes, bin, area);
      combinations.insert(bin);
    };

    // The first and last bins on this axis are exterior by definition
    for (std::size_t i :
         {static_cast<std::size_t>(0), std::get<0u>(axes).getNBins() + 1}) {
      recordExteriorBin(i);
    }

    // If no other axis is on an exterior index, stop here
    bool otherAxisExterior = false;
    for (std::size_t N = 1; N < sizeof...(Axes); ++N) {
      otherAxisExterior = otherAxisExterior | isExterior[N];
    }
    if (!otherAxisExterior) {
      return;
    }

    // Otherwise, we're on a grid border: iterate over all the other indices
    for (std::size_t i = 1; i <= std::get<0u>(axes).getNBins(); ++i) {
      recordExteriorBin(i);
    }
  }
};
/// @endcond

/// helper functions for grid-related operations
struct MultiAxisHelper {
  /// get the global indices for closest points on grid
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param bin global bin index for bin of interest
  /// @param axes actual axis objects spanning the grid
  /// @return Sorted collection of global bin indices for bins whose
  /// lower-left corners are the closest points on the grid to every
  /// point in the given bin
  ///
  /// @note @c bin must be a valid bin index (excluding under-/overflow bins
  /// along any axis).
  template <class... Axes>
  static FlatNeighborHoodIndices<sizeof...(Axes)> closestPointsIndices(
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes) {
    // get neighboring bins, but only increment.
    return neighborHoodIndices(multiIndex, std::make_pair(0, 1), axes);
  }

  /// retrieve bin center from set of local bin indices
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param multiIndex local bin indices along each axis
  /// @param axes actual axis objects spanning the grid
  /// @return center position of bin
  ///
  /// @pre @c multiIndex must only contain valid bin indices (i.e. excluding
  /// under-/overflow bins).
  template <class... Axes>
  static std::array<double, sizeof...(Axes)> getBinCenter(
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes) {
    std::array<double, sizeof...(Axes)> center{};
    MultiAxisHelperImpl<sizeof...(Axes) - 1>::getBinCenter(center, multiIndex,
                                                           axes);

    return center;
  }

  /// determine global bin index from local indices along each axis
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  ///
  /// @param localBins local bin indices along each axis
  /// @param axes actual axis objects spanning the grid
  /// @return global index for bin defined by the local bin indices
  ///
  /// @pre All local bin indices must be a valid index for the corresponding
  /// axis (including the under-/overflow bin for this axis).
  template <class... Axes>
  static std::size_t getFlatIndexFromMultiIndex(
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes) {
    std::size_t bin = 0;
    const std::size_t area = 1;

    MultiAxisHelperImpl<sizeof...(Axes) - 1>::getFlatIndexFromMultiIndex(
        multiIndex, axes, bin, area);

    return bin;
  }

  /// determine local bin index for each axis from point
  ///
  /// @tparam Point any type with point semantics supporting component access
  /// through @c operator[]
  /// @tparam Axes parameter pack of axis types defining the grid
  ///
  /// @param point point to look up in the grid
  /// @param axes actual axis objects spanning the grid
  /// @return array with local bin indices along each axis (in same order as
  /// given @c axes object)
  ///
  /// @pre The given @c Point type must represent a point in d (or higher)
  /// dimensions where d is the number of axis objects in the tuple.
  /// @note This could be a under-/overflow bin along one or more axes.
  template <class Point, class... Axes>
  static std::array<std::size_t, sizeof...(Axes)> getMultiIndexFromPoint(
      const Point& point, const std::tuple<Axes...>& axes) {
    std::array<std::size_t, sizeof...(Axes)> multiIndex{};

    MultiAxisHelperImpl<sizeof...(Axes) - 1>::getMultiIndexFromPoint(
        point, axes, multiIndex);

    return multiIndex;
  }

  /// determine local bin index for each axis from global bin index
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  ///
  /// @param bin global bin index
  /// @param axes actual axis objects spanning the grid
  /// @return array with local bin indices along each axis (in same order as
  /// given @c axes object)
  ///
  /// @note Local bin indices can contain under-/overflow bins along any axis.
  template <class... Axes>
  static std::array<std::size_t, sizeof...(Axes)> getMultiIndexFromFlatIndex(
      std::size_t bin, const std::tuple<Axes...>& axes) {
    std::size_t area = 1;
    std::array<std::size_t, sizeof...(Axes)> multiIndex{};

    MultiAxisHelperImpl<sizeof...(Axes) - 1>::getMultiIndexFromFlatIndex(
        bin, axes, area, multiIndex);

    return multiIndex;
  }

  /// retrieve lower-left bin edge from set of local bin indices
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param multiIndex local bin indices along each axis
  /// @param axes actual axis objects spanning the grid
  /// @return generalized lower-left bin edge
  ///
  /// @pre @c multiIndex must only contain valid bin indices (excluding
  /// underflow bins).
  template <class... Axes>
  static std::array<double, sizeof...(Axes)> getLowerLeftBinCorner(
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes) {
    std::array<double, sizeof...(Axes)> llEdge{};
    MultiAxisHelperImpl<sizeof...(Axes) - 1>::getLowerLeftBinCorner(
        llEdge, multiIndex, axes);

    return llEdge;
  }

  /// get local bin indices for lower-left neighboring bin
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param multiIndex local bin indices along each axis
  /// @param axes actual axis objects spanning the grid
  /// @return array with local bin indices of lower-left neighbor bin
  ///
  /// @pre @c multiIndex must only contain valid bin indices (excluding
  /// underflow bins).
  ///
  /// This function returns the local bin indices for the generalized
  /// lower-left neighbor which simply means that all local bin indices are
  /// decremented by one.
  template <class... Axes>
  static std::array<std::size_t, sizeof...(Axes)> getLowerLeftBinIndices(
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes) {
    auto llIndices = multiIndex;
    MultiAxisHelperImpl<sizeof...(Axes) - 1>::getLowerLeftBinIndices(llIndices,
                                                                     axes);

    return llIndices;
  }

  /// calculate number of bins in a grid defined by a set of axes for each axis
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param axes actual axis objects spanning the grid
  /// @return array of number of bins for each axis of the grid
  ///
  /// @note This does not include under-/overflow bins along each axis.
  template <class... Axes>
  static std::array<std::size_t, sizeof...(Axes)> getNBins(
      const std::tuple<Axes...>& axes) {
    std::array<std::size_t, sizeof...(Axes)> nBinsArray{};
    MultiAxisHelperImpl<sizeof...(Axes) - 1>::getNBins(axes, nBinsArray);
    return nBinsArray;
  }

  /// return an array with copies of the axes, converted to type AnyAxis
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param axes actual axis objects spanning the grid
  /// @return array with copies of the axis
  template <class... Axes>
  static std::array<const IAxis*, sizeof...(Axes)> getAxes(
      const std::tuple<Axes...>& axes) {
    std::array<const IAxis*, sizeof...(Axes)> arr{};
    MultiAxisHelperImpl<sizeof...(Axes) - 1>::getAxes(axes, arr);
    return arr;
  }

  /// retrieve upper-right bin edge from set of local bin indices
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param multiIndex local bin indices along each axis
  /// @param axes actual axis objects spanning the grid
  /// @return generalized upper-right bin edge
  ///
  /// @pre @c multiIndex must only contain valid bin indices (excluding
  ///      overflow bins).
  template <class... Axes>
  static std::array<double, sizeof...(Axes)> getUpperRightBinCorner(
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes) {
    std::array<double, sizeof...(Axes)> urEdge{};
    MultiAxisHelperImpl<sizeof...(Axes) - 1>::getUpperRightBinCorner(
        urEdge, multiIndex, axes);

    return urEdge;
  }

  /// get local bin indices for upper-right neighboring bin
  ///
  /// @note This function returns the local bin indices for the generalized
  /// upper-right neighbor which simply means that all local bin indices are
  /// incremented by one.
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param multiIndex local bin indices along each axis
  /// @param axes axis objects spanning the grid
  /// @return array with local bin indices of upper-right neighbor bin
  ///
  /// @pre @c multiIndex must only contain valid bin indices (excluding
  ///      overflow bins).
  template <class... Axes>
  static std::array<std::size_t, sizeof...(Axes)> getUpperRightBinIndices(
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      const std::tuple<Axes...>& axes) {
    auto urIndices = multiIndex;
    MultiAxisHelperImpl<sizeof...(Axes) - 1>::getUpperRightBinIndices(urIndices,
                                                                      axes);

    return urIndices;
  }

  /// get the minimum value of all axes of one grid
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param axes actual axis objects spanning the grid
  /// @return array returning the minima of all given axes
  template <class... Axes>
  static std::array<double, sizeof...(Axes)> getMin(
      const std::tuple<Axes...>& axes) {
    std::array<double, sizeof...(Axes)> minArray{};
    MultiAxisHelperImpl<sizeof...(Axes) - 1>::getMin(axes, minArray);
    return minArray;
  }

  /// get the maximum value of all axes of one grid
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param axes actual axis objects spanning the grid
  /// @return array returning the maxima of all given axes
  template <class... Axes>
  static std::array<double, sizeof...(Axes)> getMax(
      const std::tuple<Axes...>& axes) {
    std::array<double, sizeof...(Axes)> maxArray{};
    MultiAxisHelperImpl<sizeof...(Axes) - 1>::getMax(axes, maxArray);
    return maxArray;
  }

  /// get the bin width of all axes of one grid
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param axes actual axis objects spanning the grid
  /// @return array returning the maxima of all given axes
  template <class... Axes>
  static std::array<double, sizeof...(Axes)> getWidth(
      const std::tuple<Axes...>& axes) {
    std::array<double, sizeof...(Axes)> widthArray{};
    MultiAxisHelperImpl<sizeof...(Axes) - 1>::getWidth(axes, widthArray);
    return widthArray;
  }

  /// get global bin indices for bins in specified neighborhood
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param multiIndex local bin indices along each axis
  /// @param size size of neighborhood determining how many
  /// adjacent bins along each axis are considered
  /// @param axes actual axis objects spanning the grid
  /// @return Sorted collection of global bin indices for all bins in
  /// the neighborhood
  ///
  /// @note Over-/underflow bins are included in the neighborhood.
  /// @note The @c size parameter sets the range by how many units each local
  /// bin index is allowed to be varied. All local bin indices are
  /// varied independently, that is diagonal neighbors are included.
  /// Ignoring the truncation of the neighborhood size reaching beyond
  /// over-/underflow bins, the neighborhood is of size \f$2 \times
  /// \text{size}+1\f$ along each dimension.
  /// @note The concrete bins which are returned depend on the WrappingTypes
  /// of the contained axes
  template <class... Axes>
  static FlatNeighborHoodIndices<sizeof...(Axes)> neighborHoodIndices(
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      std::pair<std::size_t, std::size_t> sizes,
      const std::tuple<Axes...>& axes) {
    // length N array which contains local neighbors based on size par
    std::array<NeighborHoodIndices, sizeof...(Axes)> neighborIndices{};
    // get local bin indices for neighboring bins
    MultiAxisHelperImpl<sizeof...(Axes) - 1>::neighborHoodIndices(
        multiIndex, sizes, axes, neighborIndices);

    // Query the number of bins
    std::array<std::size_t, sizeof...(Axes)> nBinsArray = getNBins(axes);

    // Produce iterator of global indices
    return FlatNeighborHoodIndices(neighborIndices, nBinsArray);
  }

  template <class... Axes>
  static FlatNeighborHoodIndices<sizeof...(Axes)> neighborHoodIndices(
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      std::size_t size, const std::tuple<Axes...>& axes) {
    return neighborHoodIndices(
        multiIndex, std::make_pair(-static_cast<int>(size), size), axes);
  }

  /// get global bin indices for bins in specified neighborhood for each axis
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param multiIndex local bin indices along each axis
  /// @param size size of neighborhood for each axis, which
  /// bins along each axis are considered
  /// @param axes actual axis objects spanning the grid
  /// @return Sorted collection of global bin indices for all bins in
  /// the neighborhood
  ///
  /// @note Over-/underflow bins are included in the neighborhood.
  /// @note The @c size parameter sets the range by how many units each local
  /// bin index is allowed to be varied. All local bin indices are
  /// varied independently, that is diagonal neighbors are included.
  /// Ignoring the truncation of the neighborhood size reaching beyond
  /// over-/underflow bins, the neighborhood is of size \f$2 \times
  /// \text{size}+1\f$ along each dimension.
  /// @note The concrete bins which are returned depend on the WrappingTypes
  /// of the contained axes
  template <class... Axes>
  static FlatNeighborHoodIndices<sizeof...(Axes)> neighborHoodIndices(
      const std::array<std::size_t, sizeof...(Axes)>& multiIndex,
      std::array<std::pair<int, int>, sizeof...(Axes)>& sizes,
      const std::tuple<Axes...>& axes) {
    // length N array which contains local neighbors based on size par
    std::array<NeighborHoodIndices, sizeof...(Axes)> neighborIndices{};
    // get local bin indices for neighboring bins
    MultiAxisHelperImpl<sizeof...(Axes) - 1>::neighborHoodIndices(
        multiIndex, sizes, axes, neighborIndices);

    // Query the number of bins
    std::array<std::size_t, sizeof...(Axes)> nBinsArray = getNBins(axes);

    // Produce iterator of global indices
    return FlatNeighborHoodIndices(neighborIndices, nBinsArray);
  }

  /// get bin indices of all overflow and underflow bins
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param axes actual axis objects spanning the grid
  /// @return set of global bin indices for all over- and underflow bins
  template <class... Axes>
  static std::set<std::size_t> exteriorBinIndices(
      const std::tuple<Axes...>& axes) {
    std::array<std::size_t, sizeof...(Axes)> multiIndex{};
    std::array<bool, sizeof...(Axes)> isExterior{};
    std::set<std::size_t> combinations;
    MultiAxisHelperImpl<sizeof...(Axes) - 1>::exteriorBinIndices(
        multiIndex, isExterior, combinations, axes);

    return combinations;
  }

  /// check whether given point is inside axes limits
  ///
  /// @tparam Point any type with point semantics supporting component access
  /// through @c operator[]
  /// @tparam Axes parameter pack of axis types defining the grid
  ///
  /// @param point point to look up in the grid
  /// @param axes actual axis objects spanning the grid
  /// @return @c true if \f$\text{xmin_i} \le x_i < \text{xmax}_i \forall i=0,
  /// \dots, d-1\f$, otherwise @c false
  ///
  /// @pre The given @c Point type must represent a point in d (or higher)
  /// dimensions where d is the number of axis objects in the tuple.
  template <class Point, class... Axes>
  static bool isInside(const Point& point, const std::tuple<Axes...>& axes) {
    constexpr std::size_t MAX = sizeof...(Axes) - 1;
    return MultiAxisHelperImpl<MAX>::isInside(point, axes);
  }
};

}  // namespace Acts::detail
