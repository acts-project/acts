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

template <class... Axes>
class MultiAxis final : public IMultiAxisND<sizeof...(Axes)> {
 public:
  using Base = IMultiAxisND<sizeof...(Axes)>;

  static constexpr std::size_t DIM = Base::DIM;
  using FlatIndex = typename Base::FlatIndex;
  using MultiIndex = typename Base::MultiIndex;
  using Point = typename Base::Point;

  using AxesTuple = std::tuple<Axes...>;

  explicit MultiAxis(const std::tuple<Axes...>& axes) : m_axes(axes) {}

  explicit MultiAxis(std::tuple<Axes...>&& axes) : m_axes(std::move(axes)) {}

  explicit MultiAxis(Axes&&... axes) : m_axes(std::forward_as_tuple(axes...)) {}

  explicit MultiAxis(const Axes&... axes) : m_axes(std::tuple(axes...)) {}

  const IAxis& getAxis(std::size_t i) const override {
    return template_switch_lambda<0, DIM - 1>(
        i, [this](auto iType) -> const IAxis& {
          constexpr std::size_t iValue = decltype(iType)::value;
          return std::get<iValue>(m_axes);
        });
  }

  const AxesTuple& getAxesTuple() const { return m_axes; }

  MultiIndex getNBins() const override {
    return detail::grid_helper::getNBins(m_axes);
  }

  Point getMinPoint() const override {
    return detail::grid_helper::getMin(m_axes);
  }

  Point getMaxPoint() const override {
    return detail::grid_helper::getMax(m_axes);
  }

  bool isInside(const Point& point) const override {
    return detail::grid_helper::isInside(point, m_axes);
  }

  Point getLowerLeftBinCorner(const MultiIndex& multiIndex) const override {
    return detail::grid_helper::getLowerLeftBinEdge(multiIndex, m_axes);
  }

  Point getUpperRightBinCorner(const MultiIndex& multiIndex) const override {
    return detail::grid_helper::getUpperRightBinEdge(multiIndex, m_axes);
  }

  Point getBinCenter(const MultiIndex& multiIndex) const override {
    return detail::grid_helper::getBinCenter(multiIndex, m_axes);
  }

  FlatIndex getFlatIndexFromPoint(const Point& point) const override {
    return getFlatIndexFromMultiIndex(getMultiIndexFromPoint(point));
  }

  FlatIndex getFlatIndexFromMultiIndex(
      const MultiIndex& multiIndex) const override {
    return detail::grid_helper::getGlobalBin(multiIndex, m_axes);
  }

  MultiIndex getMultiIndexFromPoint(const Point& point) const override {
    return detail::grid_helper::getLocalBinIndices(point, m_axes);
  }

  MultiIndex getMultiIndexFromFlatIndex(FlatIndex flatIndex) const override {
    return detail::grid_helper::getLocalBinIndices(flatIndex, m_axes);
  }

  detail::GlobalNeighborHoodIndices<DIM> getNeighborHoodIndices(
      const MultiIndex& multiIndex, std::size_t size = 1u) const override {
    return detail::grid_helper::neighborHoodIndices(multiIndex, size, m_axes);
  }

  detail::GlobalNeighborHoodIndices<DIM> getNeighborHoodIndices(
      const MultiIndex& multiIndex,
      std::array<std::pair<int, int>, DIM>& sizePerAxis) const override {
    return detail::grid_helper::neighborHoodIndices(multiIndex, sizePerAxis,
                                                    m_axes);
  }

  detail::GlobalNeighborHoodIndices<DIM> getClosestPointsIndices(
      const MultiIndex& multiIndex) const override {
    return detail::grid_helper::closestPointsIndices(multiIndex, m_axes);
  }

  using Base::getClosestPointsIndices;

 private:
  std::tuple<Axes...> m_axes;
};

}  // namespace Acts
