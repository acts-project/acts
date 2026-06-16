// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/IMultiAxis.hpp"

#include <algorithm>

namespace Acts {

/// @brief Multi-dimensional binning defined by a product of one-dimensional
/// axes
///
/// This class stores a fixed set of concrete @c Axis objects in a tuple and
/// implements the @c IMultiAxisXD interface for the resulting grid. The grid
/// dimension and the concrete axis types (binning and boundary types) are
/// fixed at compile time, while the @c IMultiAxis base allows handling
/// different multi-axes through a common pointer. The grid index conventions
/// (per-axis local bin indices starting at @c 1, under-/overflow bins at @c 0
/// and <tt>nBins + 1</tt>, and flattened global bin indices including those
/// under-/overflow bins) are described on @c IMultiAxis.
///
/// @tparam Axes parameter pack of concrete @c Axis types spanning the grid
template <AxisConcept... Axes>
class MultiAxis final : public IMultiAxisXD<sizeof...(Axes)> {
 public:
  /// Base interface for this multi-axis' dimension
  using Base = IMultiAxisXD<sizeof...(Axes)>;

  /// Dimension of the grid (number of axes)
  static constexpr std::size_t DIM = Base::DIM;
  /// Flattened global bin index across all axes
  using FlatIndex = typename Base::FlatIndex;
  /// Statically sized multi-index holding one local bin index per axis
  using MultiIndex = typename Base::MultiIndex;
  /// Statically sized point holding one coordinate per axis
  using Point = typename Base::Point;

  /// Tuple type holding the concrete axes
  using AxesTuple = std::tuple<Axes...>;

  /// Construct from a tuple of axes (copy)
  /// @param axes tuple of axes spanning the grid
  explicit MultiAxis(const std::tuple<Axes...>& axes) : m_axes(axes) {}

  /// Construct from a tuple of axes (move)
  /// @param axes tuple of axes spanning the grid
  explicit MultiAxis(std::tuple<Axes...>&& axes) : m_axes(std::move(axes)) {}

  /// Construct from individual axes (forwarding)
  /// @param axes axes spanning the grid
  explicit MultiAxis(Axes&&... axes) : m_axes(std::forward_as_tuple(axes...)) {}

  /// Construct from individual axes (copy)
  /// @param axes axes spanning the grid
  explicit MultiAxis(const Axes&... axes) : m_axes(std::tuple(axes...)) {}

  /// Get the axis at the given dimension
  /// @param i index of the axis
  /// @return const reference to the requested axis
  const IAxis& getAxis(std::size_t i) const override {
    return template_switch_lambda<0, DIM - 1>(
        i, [this]<typename T>(T) -> const IAxis& {
          constexpr std::size_t iValue = T::value;
          return std::get<iValue>(m_axes);
        });
  }

  /// Get the tuple of concrete axes
  /// @return const reference to the stored axes tuple
  const AxesTuple& getAxesTuple() const { return m_axes; }

  /// Get the number of bins along each axis
  /// @return per-axis number of bins (excluding under-/overflow bins)
  MultiIndex getNBins() const override {
    return detail::MultiAxisHelper::getNBins(m_axes);
  }

  /// Get the lower boundary of the grid range along each axis
  /// @return point holding the minimum of each axis
  Point getMinPoint() const override {
    return detail::MultiAxisHelper::getMin(m_axes);
  }

  /// Get the upper boundary of the grid range along each axis
  /// @return point holding the maximum of each axis
  Point getMaxPoint() const override {
    return detail::MultiAxisHelper::getMax(m_axes);
  }

  /// Check whether a point lies inside the grid limits
  /// @param point coordinates to check, one per axis
  /// @return @c true if the point is within range along every axis
  bool isInside(const Point& point) const override {
    return detail::MultiAxisHelper::isInside(point, m_axes);
  }

  /// Get the lower-left corner of the bin given by a multi-index
  /// @param multiIndex local bin indices along each axis
  /// @return point holding the lower bin boundary of each axis
  Point getLowerLeftBinCorner(const MultiIndex& multiIndex) const override {
    return detail::MultiAxisHelper::getLowerLeftBinCorner(multiIndex, m_axes);
  }

  /// Get the upper-right corner of the bin given by a multi-index
  /// @param multiIndex local bin indices along each axis
  /// @return point holding the upper bin boundary of each axis
  Point getUpperRightBinCorner(const MultiIndex& multiIndex) const override {
    return detail::MultiAxisHelper::getUpperRightBinCorner(multiIndex, m_axes);
  }

  /// Get the center of the bin given by a multi-index
  /// @param multiIndex local bin indices along each axis
  /// @return point holding the bin center coordinate of each axis
  Point getBinCenter(const MultiIndex& multiIndex) const override {
    return detail::MultiAxisHelper::getBinCenter(multiIndex, m_axes);
  }

  /// Determine the flattened global bin index for a given point
  /// @param point coordinates to look up, one per axis
  /// @return global bin index of the bin containing the point
  FlatIndex getFlatIndexFromPoint(const Point& point) const override {
    return getFlatIndexFromMultiIndex(getMultiIndexFromPoint(point));
  }

  /// Determine the flattened global bin index from a multi-index
  /// @param multiIndex local bin indices along each axis (under-/overflow bins
  ///        are allowed)
  /// @return global bin index of the bin
  FlatIndex getFlatIndexFromMultiIndex(
      const MultiIndex& multiIndex) const override {
    return detail::MultiAxisHelper::getFlatIndexFromMultiIndex(multiIndex,
                                                               m_axes);
  }

  /// Determine the multi-index of local bin indices for a given point
  /// @param point coordinates to look up, one per axis
  /// @return local bin indices along each axis (may be under-/overflow bins)
  MultiIndex getMultiIndexFromPoint(const Point& point) const override {
    return detail::MultiAxisHelper::getMultiIndexFromPoint(point, m_axes);
  }

  /// Determine the multi-index of local bin indices from a flattened global
  /// bin index
  /// @param flatIndex global bin index
  /// @return local bin indices along each axis (may be under-/overflow bins)
  MultiIndex getMultiIndexFromFlatIndex(FlatIndex flatIndex) const override {
    return detail::MultiAxisHelper::getMultiIndexFromFlatIndex(flatIndex,
                                                               m_axes);
  }

  /// Get the global bin indices of the bins in the neighborhood of a bin
  /// @param multiIndex local bin indices of the bin of interest
  /// @param size number of adjacent bins to include along each axis (symmetric)
  /// @return sorted collection of global bin indices in the neighborhood
  detail::FlatNeighborHoodIndices<DIM> getNeighborHoodIndices(
      const MultiIndex& multiIndex, std::size_t size = 1u) const override {
    return detail::MultiAxisHelper::neighborHoodIndices(multiIndex, size,
                                                        m_axes);
  }

  /// Get the global bin indices of the bins in the neighborhood of a bin, with
  /// a separate neighborhood size per axis
  /// @param multiIndex local bin indices of the bin of interest
  /// @param sizePerAxis per-axis lower/upper number of adjacent bins to include
  /// @return sorted collection of global bin indices in the neighborhood
  detail::FlatNeighborHoodIndices<DIM> getNeighborHoodIndices(
      const MultiIndex& multiIndex,
      std::array<std::pair<int, int>, DIM>& sizePerAxis) const override {
    return detail::MultiAxisHelper::neighborHoodIndices(multiIndex, sizePerAxis,
                                                        m_axes);
  }

  /// Get the global bin indices of the grid points closest to the given bin
  /// @param multiIndex local bin indices of the bin of interest
  /// @return sorted collection of global bin indices whose lower-left corners
  ///         are the closest grid points to every point in the given bin
  detail::FlatNeighborHoodIndices<DIM> getClosestPointsIndices(
      const MultiIndex& multiIndex) const override {
    return detail::MultiAxisHelper::closestPointsIndices(multiIndex, m_axes);
  }

  using Base::getClosestPointsIndices;

 private:
  /// tuple of concrete axes spanning the grid
  std::tuple<Axes...> m_axes;
};

}  // namespace Acts
