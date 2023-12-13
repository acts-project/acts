// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Holders.hpp"

#include <array>

namespace Acts {

// Using Global iterator, including over/under flow bins
template <typename T, class... Axes>
class GridGlobalIterator {
 public:
  static constexpr std::size_t DIM = sizeof...(Axes);

  using iterator_category = std::random_access_iterator_tag;
  using value_type = T;
  using difference_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;

  GridGlobalIterator() = default;
  GridGlobalIterator(Acts::Grid<T, Axes...>&& grid, std::size_t idx) = delete;
  GridGlobalIterator(const Acts::Grid<T, Axes...>& grid, std::size_t idx = 0ul);

  GridGlobalIterator(const GridGlobalIterator<T, Axes...>& other) = default;
  GridGlobalIterator<T, Axes...>& operator=(
      const GridGlobalIterator<T, Axes...>& other) = default;

  GridGlobalIterator(GridGlobalIterator<T, Axes...>&& other) noexcept;
  GridGlobalIterator<T, Axes...>& operator=(
      GridGlobalIterator<T, Axes...>&& other) noexcept;

  ~GridGlobalIterator() = default;

  bool operator==(const GridGlobalIterator<T, Axes...>& other) const;
  bool operator!=(const GridGlobalIterator<T, Axes...>& other) const;

  bool operator<(const GridGlobalIterator<T, Axes...>& other) const;
  bool operator>(const GridGlobalIterator<T, Axes...>& other) const;
  bool operator<=(const GridGlobalIterator<T, Axes...>& other) const;
  bool operator>=(const GridGlobalIterator<T, Axes...>& other) const;

  GridGlobalIterator<T, Axes...>& operator+=(const std::size_t offset);
  GridGlobalIterator<T, Axes...>& operator-=(const std::size_t offset);
  GridGlobalIterator<T, Axes...> operator+(const std::size_t offset) const;
  GridGlobalIterator<T, Axes...> operator-(const std::size_t offset) const;

  difference_type operator-(const GridGlobalIterator<T, Axes...>& other) const;
  const value_type& operator*() const;

  GridGlobalIterator<T, Axes...>& operator++();
  GridGlobalIterator<T, Axes...> operator++(int);

 private:
  Acts::detail::RefHolder<const Acts::Grid<T, Axes...>> m_grid{nullptr};
  std::size_t m_idx{0ul};
};

// Using Local iterator, excluding over/under flow bins
// Can also allow for custom navigation pattern along axes
template <typename T, class... Axes>
class GridLocalIterator {
 public:
  static constexpr std::size_t DIM = sizeof...(Axes);

  using iterator_category = std::bidirectional_iterator_tag;
  using value_type = T;
  using difference_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;

  GridLocalIterator() = default;
  GridLocalIterator(Acts::Grid<T, Axes...>&& grid,
                    const std::array<std::size_t, DIM>& indexes) = delete;
  GridLocalIterator(Acts::Grid<T, Axes...>&& grid,
                    const std::array<std::size_t, DIM>& indexes,
                    std::array<std::vector<std::size_t>, DIM> navigation) =
      delete;
  GridLocalIterator(const Acts::Grid<T, Axes...>& grid,
                    const std::array<std::size_t, DIM>& indexes);
  GridLocalIterator(const Acts::Grid<T, Axes...>& grid,
                    const std::array<std::size_t, DIM>& indexes,
                    std::array<std::vector<std::size_t>, DIM> navigation);

  GridLocalIterator(const GridLocalIterator<T, Axes...>& other) = default;
  GridLocalIterator<T, Axes...>& operator=(
      const GridLocalIterator<T, Axes...>& other) = default;

  GridLocalIterator(GridLocalIterator<T, Axes...>&& other) noexcept;
  GridLocalIterator<T, Axes...>& operator=(
      GridLocalIterator<T, Axes...>&& other) noexcept;

  ~GridLocalIterator() = default;

  bool operator==(const Acts::GridLocalIterator<T, Axes...>& other) const;
  bool operator!=(const Acts::GridLocalIterator<T, Axes...>& other) const;

  const value_type& operator*() const;

  GridLocalIterator<T, Axes...>& operator++();
  GridLocalIterator<T, Axes...> operator++(int);

  std::array<std::size_t, DIM> localPosition() const;

 private:
  template <std::size_t N>
  void increment();

 private:
  Acts::detail::RefHolder<const Acts::Grid<T, Axes...>> m_grid{nullptr};
  std::array<std::size_t, DIM> m_numLocalBins{};
  std::array<std::size_t, DIM> m_currentIndex{};
  std::array<std::vector<std::size_t>, DIM> m_navigationIndex{};
};

}  // namespace Acts

#include "Acts/Utilities/GridIterator.ipp"
