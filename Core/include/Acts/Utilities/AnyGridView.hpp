// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <type_traits>

namespace Acts {

namespace detail {

template <typename T, bool isConst>
class AnyGridViewBase {
 public:
  using GridPointerType = std::conditional_t<isConst, const IGrid*, IGrid*>;
  using AnyIndexType = IGrid::AnyIndexType;
  using AnyPointType = IGrid::AnyPointType;

  explicit AnyGridViewBase(const IGrid& grid)
    requires(isConst)
      : m_grid(&grid) {
    checkType();
  }

  explicit AnyGridViewBase(IGrid& grid) : m_grid(&grid) { checkType(); }

  template <typename... Axes>
  AnyGridViewBase(Grid<T, Axes...>& grid)
    requires(!isConst)
      : m_grid(&grid) {
    checkType();
  }

  template <typename... Axes>
  AnyGridViewBase(const Grid<T, Axes...>& grid)
    requires(isConst)
      : m_grid(&grid) {
    checkType();
  }

  AnyGridViewBase(const AnyGridViewBase& other) = default;
  AnyGridViewBase& operator=(const AnyGridViewBase& other) = default;

  AnyGridViewBase(AnyGridViewBase&&) noexcept = default;
  AnyGridViewBase& operator=(AnyGridViewBase&&) noexcept = default;

  T& atLocalBins(const AnyIndexType& indices)
    requires(!isConst)
  {
    std::any any = m_grid->atLocalBinsAny(indices);
    return *std::any_cast<T*>(any);
  }

  const T& atLocalBins(const AnyIndexType& indices) const {
    std::any any = m_grid->atLocalBinsAny(indices);
    return *std::any_cast<T const*>(any);
  }

  std::size_t dimensions() const { return m_grid->dimensions(); }

  AnyPointType binCenter(const IGrid::AnyIndexType& indices) const {
    return m_grid->binCenterAny(indices);
  }

  AnyPointType lowerLeftBinEdge(const IGrid::AnyIndexType& indices) const {
    return m_grid->lowerLeftBinEdgeAny(indices);
  }

  AnyPointType upperRightBinEdge(const IGrid::AnyIndexType& indices) const {
    return m_grid->upperRightBinEdgeAny(indices);
  }

  AnyIndexType numLocalBins() const { return m_grid->numLocalBinsAny(); }

 private:
  void checkType() {
    if (m_grid->valueType() != typeid(T)) {
      throw std::invalid_argument("Type mismatch between grid and view type");
    }
  }

  GridPointerType m_grid;
};

}  // namespace detail

template <typename T>
class AnyGridView : public detail::AnyGridViewBase<T, false> {
 public:
  using detail::AnyGridViewBase<T, false>::AnyGridViewBase;
};

template <typename T, typename... Axes>
AnyGridView(Grid<T, Axes...>& grid) -> AnyGridView<T>;

template <typename T>
class AnyGridConstView : public detail::AnyGridViewBase<T, true> {
 public:
  using detail::AnyGridViewBase<T, true>::AnyGridViewBase;
};

template <typename T, typename... Axes>
AnyGridConstView(const Grid<T, Axes...>& grid) -> AnyGridConstView<T>;

}  // namespace Acts
