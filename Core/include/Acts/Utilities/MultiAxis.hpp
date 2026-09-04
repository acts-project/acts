// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
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
class MultiAxis : public IMultiAxisXD<sizeof...(Axes)> {
 public:
  /// Base interface for this multi-axis' dimension
  using Base = IMultiAxisXD<sizeof...(Axes)>;

  /// Dimension of the grid (number of axes)
  static constexpr std::size_t DIM = Base::DIM;
  /// Flattened global bin index across all axes
  using GlobalBin = typename Base::GlobalBin;
  /// Statically sized multi-index holding one local bin index per axis
  using LocalBins = typename Base::LocalBins;
  /// Statically sized point holding one coordinate per axis
  using Point = typename Base::Point;
  /// Point holding one coordinate per axis as an Eigen (algebra) vector
  using VectorPoint = Vector<DIM>;

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

  /// @copydoc IMultiAxisXD::getAxis(std::size_t) const
  const IAxis& getAxis(std::size_t i) const final {
    return template_switch_lambda<0, DIM - 1>(
        i, [this]<typename T>(T) -> const IAxis& {
          constexpr std::size_t iValue = T::value;
          return std::get<iValue>(m_axes);
        });
  }

  /// Get the tuple of concrete axes
  /// @return const reference to the stored axes tuple
  const AxesTuple& getAxesTuple() const { return m_axes; }

  /// @copydoc IMultiAxisXD::getNBins() const
  LocalBins getNBins() const final {
    return detail::MultiAxisHelper::getNBins(m_axes);
  }

  /// @copydoc IMultiAxisXD::getMinPoint() const
  Point getMinPoint() const final {
    return detail::MultiAxisHelper::getMin(m_axes);
  }

  /// @copydoc IMultiAxisXD::getMaxPoint() const
  Point getMaxPoint() const final {
    return detail::MultiAxisHelper::getMax(m_axes);
  }

  /// @copydoc IMultiAxisXD::isInside(const Point&) const
  bool isInside(const Point& point) const final {
    return detail::MultiAxisHelper::isInside(point, m_axes);
  }

  /// Check whether a point lies inside the grid limits
  /// @tparam Derived Eigen expression type of the point (e.g. @c VectorPoint)
  /// @param point coordinates to check, one per axis
  /// @return @c true if the point is within range along every axis
  template <typename Derived>
  bool isInside(const Eigen::MatrixBase<Derived>& point) const {
    return detail::MultiAxisHelper::isInside(point.derived(), m_axes);
  }

  /// @copydoc IMultiAxisXD::getLowerLeftBinEdge(const LocalBins&) const
  Point getLowerLeftBinEdge(const LocalBins& localBins) const final {
    return detail::MultiAxisHelper::getLowerLeftBinEdge(localBins, m_axes);
  }

  /// @copydoc IMultiAxisXD::getUpperRightBinEdge(const LocalBins&) const
  Point getUpperRightBinEdge(const LocalBins& localBins) const final {
    return detail::MultiAxisHelper::getUpperRightBinEdge(localBins, m_axes);
  }

  /// @copydoc IMultiAxisXD::getBinCenter(const LocalBins&) const
  Point getBinCenter(const LocalBins& localBins) const final {
    return detail::MultiAxisHelper::getBinCenter(localBins, m_axes);
  }

  /// @copydoc IMultiAxisXD::getBinWidth(const LocalBins&) const
  Point getBinWidth(const LocalBins& localBins) const final {
    return detail::MultiAxisHelper::getBinWidth(localBins, m_axes);
  }

  /// @copydoc IMultiAxisXD::getGlobalBinFromPoint(const Point&) const
  GlobalBin getGlobalBinFromPoint(const Point& point) const final {
    return getGlobalBinFromLocalBins(getLocalBinsFromPoint(point));
  }

  /// Determine the flattened global bin index for a given point
  /// @tparam Derived Eigen expression type of the point (e.g. @c VectorPoint)
  /// @param point coordinates to look up, one per axis
  /// @return global bin index of the bin containing the point
  template <typename Derived>
  GlobalBin getGlobalBinFromPoint(
      const Eigen::MatrixBase<Derived>& point) const {
    return getGlobalBinFromLocalBins(getLocalBinsFromPoint(point));
  }

  /// Determine the flattened global bin index of the bin with the lower left
  /// edge closest to the given point for each axis
  /// @param point coordinates to look up, one per axis
  /// @return global bin index of the bin
  GlobalBin getGlobalBinFromLowerLeftEdge(const Point& point) const {
    return getGlobalBinFromLocalBins(getLocalBinsFromLowerLeftEdge(point));
  }

  /// @copydoc getGlobalBinFromLowerLeftEdge(const Point&) const
  /// @tparam Derived Eigen expression type of the point (e.g. @c VectorPoint)
  template <typename Derived>
  GlobalBin getGlobalBinFromLowerLeftEdge(
      const Eigen::MatrixBase<Derived>& point) const {
    return getGlobalBinFromLocalBins(getLocalBinsFromLowerLeftEdge(point));
  }

  /// @copydoc IMultiAxisXD::getGlobalBinFromLocalBins(const LocalBins&) const
  GlobalBin getGlobalBinFromLocalBins(const LocalBins& localBins) const final {
    return detail::MultiAxisHelper::getGlobalBinFromLocalBins(localBins,
                                                              m_axes);
  }

  /// @copydoc IMultiAxisXD::getLocalBinsFromPoint(const Point&) const
  LocalBins getLocalBinsFromPoint(const Point& point) const final {
    return detail::MultiAxisHelper::getLocalBinsFromPoint(point, m_axes);
  }

  /// Determine the multi-index of local bin indices for a given point
  /// @tparam Derived Eigen expression type of the point (e.g. @c VectorPoint)
  /// @param point coordinates to look up, one per axis
  /// @return local bin indices along each axis (may be under-/overflow bins)
  template <typename Derived>
  LocalBins getLocalBinsFromPoint(
      const Eigen::MatrixBase<Derived>& point) const {
    return detail::MultiAxisHelper::getLocalBinsFromPoint(point.derived(),
                                                          m_axes);
  }

  /// Determine the multi-index of local bin indices for a given point, where
  /// the point is interpreted as being shifted by half a bin width along each
  /// axis.
  /// @param point coordinates to look up, one per axis
  /// @return local bin indices along each axis (may be under-/overflow bins)
  LocalBins getLocalBinsFromLowerLeftEdge(const Point& point) const {
    return detail::MultiAxisHelper::getLocalBinsFromLowerLeftEdge(point,
                                                                  m_axes);
  }

  /// @copydoc getLocalBinsFromLowerLeftEdge(const Point&) const
  /// @tparam Derived Eigen expression type of the point (e.g. @c VectorPoint)
  template <typename Derived>
  LocalBins getLocalBinsFromLowerLeftEdge(
      const Eigen::MatrixBase<Derived>& point) const {
    return detail::MultiAxisHelper::getLocalBinsFromLowerLeftEdge(
        point.derived(), m_axes);
  }

  /// @copydoc IMultiAxisXD::getLocalBinsFromGlobalBin(GlobalBin) const
  LocalBins getLocalBinsFromGlobalBin(GlobalBin globalBin) const final {
    return detail::MultiAxisHelper::getLocalBinsFromGlobalBin(globalBin,
                                                              m_axes);
  }

  /// @copydoc IMultiAxisXD::getNeighborHoodIndices(const LocalBins&, std::size_t) const
  detail::FlatNeighborHoodIndices<DIM> getNeighborHoodIndices(
      const LocalBins& localBins, std::size_t size = 1u) const final {
    return detail::MultiAxisHelper::neighborHoodIndices(localBins, size,
                                                        m_axes);
  }

  /// @copydoc IMultiAxisXD::getNeighborHoodIndices(const LocalBins&, const std::pair<int, int>&) const
  detail::FlatNeighborHoodIndices<DIM> getNeighborHoodIndices(
      const LocalBins& localBins, const std::pair<int, int>& size) const final {
    return detail::MultiAxisHelper::neighborHoodIndices(localBins, size,
                                                        m_axes);
  }

  /// @copydoc IMultiAxisXD::getNeighborHoodIndices(const LocalBins&, const std::array<std::pair<int, int>, DIM>&) const
  detail::FlatNeighborHoodIndices<DIM> getNeighborHoodIndices(
      const LocalBins& localBins,
      const std::array<std::pair<int, int>, DIM>& sizePerAxis) const final {
    return detail::MultiAxisHelper::neighborHoodIndices(localBins, sizePerAxis,
                                                        m_axes);
  }

  /// @copydoc IMultiAxisXD::getClosestPointsIndices(const LocalBins&) const
  detail::FlatNeighborHoodIndices<DIM> getClosestPointsIndices(
      const LocalBins& localBins) const override {
    return getNeighborHoodIndices(localBins, {0, 1});
  }

  /// @copydoc IMultiAxisXD::getClosestPointsIndices(const Point&) const
  detail::FlatNeighborHoodIndices<DIM> getClosestPointsIndices(
      const Point& point) const override {
    return getNeighborHoodIndices(getLocalBinsFromPoint(point), {0, 1});
  }

  /// Get the global bin indices of the grid points closest to the given point
  /// @tparam Derived Eigen expression type of the point (e.g. @c VectorPoint)
  /// @param point coordinates to look up, one per axis
  /// @return sorted collection of global bin indices of the closest grid points
  template <typename Derived>
  detail::FlatNeighborHoodIndices<DIM> getClosestPointsIndices(
      const Eigen::MatrixBase<Derived>& point) const {
    return getNeighborHoodIndices(getLocalBinsFromPoint(point), {0, 1});
  }

 private:
  /// tuple of concrete axes spanning the grid
  std::tuple<Axes...> m_axes;
};

}  // namespace Acts
