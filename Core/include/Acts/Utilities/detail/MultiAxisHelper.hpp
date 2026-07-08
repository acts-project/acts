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
  using size_type = std::size_t;

  /// You can get the neighbor multi indices from
  /// MultiAxisHelper::neighborHoodIndices and the number of bins in
  /// each direction from MultiAxisHelper::getNBins.
  FlatNeighborHoodIndices(std::array<NeighborHoodIndices, DIM>& neighborIndices,
                          const std::array<std::size_t, DIM>& nBinsArray)
      : m_localBins(neighborIndices) {
    if constexpr (DIM == 1) {
      return;
    } else {
      std::size_t flatStride = 1;
      for (long i = DIM - 2; i >= 0; --i) {
        flatStride *= (nBinsArray[i + 1] + 2);
        m_flatStride[i] = flatStride;
      }
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
        : m_localBinsIter(std::move(multiIndicesIter)), m_parent(&parent) {}

    std::size_t operator*() const {
      std::size_t globalBin = *m_localBinsIter[DIM - 1];
      if constexpr (DIM == 1) {
        return globalBin;
      } else {
        for (std::size_t i = 0; i < DIM - 1; ++i) {
          globalBin += m_parent->m_flatStride[i] * (*m_localBinsIter[i]);
        }
        return globalBin;
      }
    }

    iterator& operator++() {
      const auto& multiIndices = m_parent->m_localBins;

      // Go to the next flat index via a lexicographic increment:
      // - Start by incrementing the last multi index
      // - If it reaches the end, reset it and increment the previous one...
      for (long i = DIM - 1; i > 0; --i) {
        ++m_localBinsIter[i];
        if (m_localBinsIter[i] != multiIndices[i].end()) {
          return *this;
        }
        m_localBinsIter[i] = multiIndices[i].begin();
      }

      // The first index should stay at the end value when it reaches it, so
      // that we know when we've reached the end of iteration.
      ++m_localBinsIter[0];
      return *this;
    }

    iterator operator++(int) {
      iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool isEqual(const iterator& b) const {
      if (b.m_parent == nullptr) {
        return m_localBinsIter[0] == m_parent->m_localBins[0].end();
      } else {
        return m_localBinsIter == b.m_localBinsIter;
      }
    }

    friend bool operator==(const iterator& a, const iterator& b) {
      return a.isEqual(b);
    }

   private:
    std::array<NeighborHoodIndices::iterator, DIM> m_localBinsIter;
    const FlatNeighborHoodIndices* m_parent = nullptr;
  };

  iterator begin() const {
    std::array<NeighborHoodIndices::iterator, DIM> multiIndicesIter{};
    for (std::size_t i = 0; i < DIM; ++i) {
      multiIndicesIter[i] = m_localBins[i].begin();
    }
    return iterator(*this, std::move(multiIndicesIter));
  }

  iterator end() const { return iterator(); }

  /// Number of indices that will be produced if this sequence is iterated
  std::size_t size() const {
    std::size_t result = m_localBins[0].size();
    for (std::size_t i = 1; i < DIM; ++i) {
      result *= m_localBins[i].size();
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
  std::array<NeighborHoodIndices, DIM> m_localBins{};
  std::array<std::size_t, DIM - 1> m_flatStride{};
};

/// helper functions for grid-related operations
struct MultiAxisHelper {
 private:
  /// Invoke @p f with a compile-time axis index (as an explicit template
  /// argument) for every axis @c 0, 1, ..., sizeof...(Is)-1 in order.
  template <std::size_t... Is, class F>
  static constexpr void forEachAxis(std::index_sequence<Is...> /*seq*/, F&& f) {
    (f.template operator()<Is>(), ...);
  }

  /// Same as @c forEachAxis but visiting the axes in reverse order, i.e.
  /// sizeof...(Is)-1, ..., 1, 0. Used where a running product over the trailing
  /// axes has to be accumulated.
  template <std::size_t... Is, class F>
  static constexpr void forEachAxisReverse(std::index_sequence<Is...> /*seq*/,
                                           F&& f) {
    constexpr std::size_t N = sizeof...(Is);
    (f.template operator()<N - 1 - Is>(), ...);
  }

 public:
  /// retrieve bin center from set of local bin indices
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param localBins local bin indices along each axis
  /// @param axes actual axis objects spanning the grid
  /// @return center position of bin
  ///
  /// @pre @c localBins must only contain valid bin indices (i.e. excluding
  /// under-/overflow bins).
  template <class... Axes>
  static std::array<double, sizeof...(Axes)> getBinCenter(
      const std::array<std::size_t, sizeof...(Axes)>& localBins,
      const std::tuple<Axes...>& axes) {
    std::array<double, sizeof...(Axes)> center{};
    forEachAxis(std::index_sequence_for<Axes...>{}, [&]<std::size_t i>() {
      center[i] = std::get<i>(axes).getBinCenter(localBins[i]);
    });
    return center;
  }

  /// get the bin width of all axes of one grid
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param axes actual axis objects spanning the grid
  /// @return array returning the maxima of all given axes
  template <class... Axes>
  static std::array<double, sizeof...(Axes)> getBinWidth(
      const std::array<std::size_t, sizeof...(Axes)>& localBins,
      const std::tuple<Axes...>& axes) {
    std::array<double, sizeof...(Axes)> widthArray{};
    forEachAxis(std::index_sequence_for<Axes...>{}, [&]<std::size_t i>() {
      widthArray[i] = std::get<i>(axes).getBinWidth(localBins[i]);
    });
    return widthArray;
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
  static std::size_t getGlobalBinFromLocalBins(
      const std::array<std::size_t, sizeof...(Axes)>& localBins,
      const std::tuple<Axes...>& axes) {
    // The trailing axis is the fastest running index. The stride for axis i is
    // the product of (nBins + 2) over all axes following it (the +2 accounts
    // for the under-/overflow bins), accumulated by walking the axes in
    // reverse.
    std::size_t bin = 0;
    std::size_t area = 1;
    forEachAxisReverse(std::index_sequence_for<Axes...>{},
                       [&]<std::size_t i>() {
                         bin += area * localBins[i];
                         area *= std::get<i>(axes).getNBins() + 2;
                       });
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
  static std::array<std::size_t, sizeof...(Axes)> getLocalBinsFromPoint(
      const Point& point, const std::tuple<Axes...>& axes) {
    std::array<std::size_t, sizeof...(Axes)> localBins{};
    forEachAxis(std::index_sequence_for<Axes...>{}, [&]<std::size_t i>() {
      localBins[i] = std::get<i>(axes).getBin(point[i]);
    });
    return localBins;
  }

  /// determine local bin index of the bin with the lower left edge closest to
  /// the given point for each axis
  ///
  /// @tparam Point any type with point semantics supporting component access
  /// through @c operator[]
  /// @tparam Axes parameter pack of axis types defining the grid
  ///
  /// @param point point to look up in the grid
  /// @return array with local bin indices along each axis (in same order as
  /// given @c axes object)
  ///
  /// @pre The given @c Point type must represent a point in d (or higher)
  /// dimensions where d is dimensionality of the grid.
  /// @note This could be a under-/overflow bin along one or more axes.
  template <class Point, class... Axes>
  static std::array<std::size_t, sizeof...(Axes)> getLocalBinsFromLowerLeftEdge(
      const Point& point, const std::tuple<Axes...>& axes) {
    constexpr std::size_t DIM = sizeof...(Axes);

    const auto localBins =
        detail::MultiAxisHelper::getLocalBinsFromPoint(point, axes);

    Point shiftedPoint;
    const auto width = detail::MultiAxisHelper::getBinWidth(localBins, axes);
    for (std::size_t i = 0; i < DIM; i++) {
      shiftedPoint[i] = point[i] + width[i] / 2;
    }

    return detail::MultiAxisHelper::getLocalBinsFromPoint(shiftedPoint, axes);
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
  static std::array<std::size_t, sizeof...(Axes)> getLocalBinsFromGlobalBin(
      std::size_t bin, const std::tuple<Axes...>& axes) {
    using Seq = std::index_sequence_for<Axes...>;

    // Compute the per-axis stride (same layout as getGlobalBinFromLocalBins),
    // then peel off the local bin indices from the most significant axis down.
    std::array<std::size_t, sizeof...(Axes)> strides{};
    std::size_t area = 1;
    forEachAxisReverse(Seq{}, [&]<std::size_t i>() {
      strides[i] = area;
      area *= std::get<i>(axes).getNBins() + 2;
    });

    std::array<std::size_t, sizeof...(Axes)> localBins{};
    forEachAxis(Seq{}, [&]<std::size_t i>() {
      localBins[i] = bin / strides[i];
      bin %= strides[i];
    });
    return localBins;
  }

  /// retrieve lower-left bin edge from set of local bin indices
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param localBins local bin indices along each axis
  /// @param axes actual axis objects spanning the grid
  /// @return generalized lower-left bin edge
  ///
  /// @pre @c localBins must only contain valid bin indices (excluding
  /// underflow bins).
  template <class... Axes>
  static std::array<double, sizeof...(Axes)> getLowerLeftBinEdge(
      const std::array<std::size_t, sizeof...(Axes)>& localBins,
      const std::tuple<Axes...>& axes) {
    std::array<double, sizeof...(Axes)> llEdge{};
    forEachAxis(std::index_sequence_for<Axes...>{}, [&]<std::size_t i>() {
      llEdge[i] = std::get<i>(axes).getBinLowerBound(localBins[i]);
    });
    return llEdge;
  }

  /// get local bin indices for lower-left neighboring bin
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param localBins local bin indices along each axis
  /// @param axes actual axis objects spanning the grid
  /// @return array with local bin indices of lower-left neighbor bin
  ///
  /// @pre @c localBins must only contain valid bin indices (excluding
  /// underflow bins).
  ///
  /// This function returns the local bin indices for the generalized
  /// lower-left neighbor which simply means that all local bin indices are
  /// decremented by one.
  template <class... Axes>
  static std::array<std::size_t, sizeof...(Axes)> getLowerLeftBinIndices(
      const std::array<std::size_t, sizeof...(Axes)>& localBins,
      const std::tuple<Axes...>& axes) {
    auto llIndices = localBins;
    forEachAxis(std::index_sequence_for<Axes...>{}, [&]<std::size_t i>() {
      llIndices[i] = std::get<i>(axes).wrapBin(llIndices[i] - 1);
    });
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
    // by convention getNBins does not include under-/overflow bins
    forEachAxis(std::index_sequence_for<Axes...>{}, [&]<std::size_t i>() {
      nBinsArray[i] = std::get<i>(axes).getNBins();
    });
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
    forEachAxis(std::index_sequence_for<Axes...>{}, [&]<std::size_t i>() {
      arr[i] = static_cast<const IAxis*>(&std::get<i>(axes));
    });
    return arr;
  }

  /// retrieve upper-right bin edge from set of local bin indices
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param localBins local bin indices along each axis
  /// @param axes actual axis objects spanning the grid
  /// @return generalized upper-right bin edge
  ///
  /// @pre @c localBins must only contain valid bin indices (excluding
  ///      overflow bins).
  template <class... Axes>
  static std::array<double, sizeof...(Axes)> getUpperRightBinEdge(
      const std::array<std::size_t, sizeof...(Axes)>& localBins,
      const std::tuple<Axes...>& axes) {
    std::array<double, sizeof...(Axes)> urEdge{};
    forEachAxis(std::index_sequence_for<Axes...>{}, [&]<std::size_t i>() {
      urEdge[i] = std::get<i>(axes).getBinUpperBound(localBins[i]);
    });
    return urEdge;
  }

  /// get local bin indices for upper-right neighboring bin
  ///
  /// @note This function returns the local bin indices for the generalized
  /// upper-right neighbor which simply means that all local bin indices are
  /// incremented by one.
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param localBins local bin indices along each axis
  /// @param axes axis objects spanning the grid
  /// @return array with local bin indices of upper-right neighbor bin
  ///
  /// @pre @c localBins must only contain valid bin indices (excluding
  ///      overflow bins).
  template <class... Axes>
  static std::array<std::size_t, sizeof...(Axes)> getUpperRightBinIndices(
      const std::array<std::size_t, sizeof...(Axes)>& localBins,
      const std::tuple<Axes...>& axes) {
    auto urIndices = localBins;
    forEachAxis(std::index_sequence_for<Axes...>{}, [&]<std::size_t i>() {
      urIndices[i] = std::get<i>(axes).wrapBin(urIndices[i] + 1);
    });
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
    forEachAxis(std::index_sequence_for<Axes...>{}, [&]<std::size_t i>() {
      minArray[i] = std::get<i>(axes).getMin();
    });
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
    forEachAxis(std::index_sequence_for<Axes...>{}, [&]<std::size_t i>() {
      maxArray[i] = std::get<i>(axes).getMax();
    });
    return maxArray;
  }

  /// get global bin indices for bins in specified neighborhood
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param localBins local bin indices along each axis
  /// @param sizes size of neighborhood determining how many
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
      const std::array<std::size_t, sizeof...(Axes)>& localBins,
      std::pair<std::size_t, std::size_t> sizes,
      const std::tuple<Axes...>& axes) {
    // length N array which contains local neighbors based on size par
    std::array<NeighborHoodIndices, sizeof...(Axes)> neighborIndices{};
    // get local bin indices for neighboring bins (same size on every axis)
    forEachAxis(std::index_sequence_for<Axes...>{}, [&]<std::size_t i>() {
      neighborIndices[i] =
          std::get<i>(axes).neighborHoodIndices(localBins[i], sizes);
    });

    // Query the number of bins
    std::array<std::size_t, sizeof...(Axes)> nBinsArray = getNBins(axes);

    // Produce iterator of global indices
    return FlatNeighborHoodIndices(neighborIndices, nBinsArray);
  }

  template <class... Axes>
  static FlatNeighborHoodIndices<sizeof...(Axes)> neighborHoodIndices(
      const std::array<std::size_t, sizeof...(Axes)>& localBins,
      std::size_t size, const std::tuple<Axes...>& axes) {
    return neighborHoodIndices(
        localBins, std::make_pair(-static_cast<int>(size), size), axes);
  }

  /// get global bin indices for bins in specified neighborhood for each axis
  ///
  /// @tparam Axes parameter pack of axis types defining the grid
  /// @param localBins local bin indices along each axis
  /// @param sizes size of neighborhood for each axis, which
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
      const std::array<std::size_t, sizeof...(Axes)>& localBins,
      const std::array<std::pair<int, int>, sizeof...(Axes)>& sizes,
      const std::tuple<Axes...>& axes) {
    // length N array which contains local neighbors based on size par
    std::array<NeighborHoodIndices, sizeof...(Axes)> neighborIndices{};
    // get local bin indices for neighboring bins (per-axis size)
    forEachAxis(std::index_sequence_for<Axes...>{}, [&]<std::size_t i>() {
      neighborIndices[i] =
          std::get<i>(axes).neighborHoodIndices(localBins[i], sizes[i]);
    });

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
    constexpr std::size_t DIM = sizeof...(Axes);
    // getNBins excludes the under-/overflow bins; valid indices run 0..nBins+1.
    const std::array<std::size_t, DIM> nBins = getNBins(axes);

    std::set<std::size_t> combinations;
    std::array<std::size_t, DIM> localBins{};
    const auto record = [&](std::size_t i0) {
      localBins[0] = i0;
      combinations.insert(getGlobalBinFromLocalBins(localBins, axes));
    };

    // The under-/overflow bins of axis 0 are exterior for any combination of
    // the remaining axes. We only need to walk the interior bins of axis 0 when
    // some other axis already sits on an exterior index. Run an odometer over
    // axes 1..DIM-1 (over their full 0..nBins+1 range) and handle axis 0
    // explicitly for each combination.
    while (true) {
      record(0);
      record(nBins[0] + 1);

      bool otherExterior = false;
      for (std::size_t d = 1; d < DIM; ++d) {
        if (localBins[d] == 0 || localBins[d] == nBins[d] + 1) {
          otherExterior = true;
          break;
        }
      }
      if (otherExterior) {
        for (std::size_t i = 1; i <= nBins[0]; ++i) {
          record(i);
        }
      }

      // Increment the odometer over axes 1..DIM-1; stop once it wraps around.
      std::size_t d = 1;
      for (; d < DIM; ++d) {
        if (++localBins[d] <= nBins[d] + 1) {
          break;
        }
        localBins[d] = 0;
      }
      if (d == DIM) {
        break;
      }
    }

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
    return [&]<std::size_t... Is>(std::index_sequence<Is...> /*seq*/) {
      return (std::get<Is>(axes).isInside(point[Is]) && ...);
    }(std::index_sequence_for<Axes...>{});
  }
};

}  // namespace Acts::detail
