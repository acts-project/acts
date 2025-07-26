// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <iterator>

namespace Acts {

template <typename Index>
class IndexRange {
 public:
  using index_type = Index;
  using index_range_type = std::pair<index_type, index_type>;

  class Iterator {
   public:
    using value_type = Index;
    using difference_type = std::ptrdiff_t;
    using pointer = Index *;
    using reference = Index &;

    using iterator_category = std::random_access_iterator_tag;
    using iterator_concept = std::random_access_iterator_tag;

    constexpr Iterator() noexcept = default;
    explicit constexpr Iterator(Index index) noexcept : m_index{index} {}

    constexpr value_type operator*() const noexcept { return m_index; }
    constexpr value_type operator[](difference_type n) const noexcept {
      return m_index + n;
    }

    constexpr Iterator &operator++() noexcept {
      ++m_index;
      return *this;
    }
    constexpr Iterator operator++(int) noexcept {
      auto tmp = *this;
      ++(*this);
      return tmp;
    }
    constexpr Iterator &operator--() noexcept {
      --m_index;
      return *this;
    }
    constexpr Iterator operator--(int) noexcept {
      auto tmp = *this;
      --(*this);
      return tmp;
    }

    constexpr Iterator &operator+=(difference_type n) noexcept {
      m_index += n;
      return *this;
    }
    constexpr Iterator &operator-=(difference_type n) noexcept {
      m_index -= n;
      return *this;
    }

   private:
    Index m_index{0};

    friend constexpr Iterator operator+(Iterator it,
                                        difference_type n) noexcept {
      return it += n;
    }

    friend constexpr Iterator operator+(difference_type n,
                                        Iterator it) noexcept {
      return it += n;
    }

    friend constexpr Iterator operator-(Iterator it,
                                        difference_type n) noexcept {
      return it -= n;
    }

    friend constexpr difference_type operator-(const Iterator &lhs,
                                               const Iterator &rhs) noexcept {
      return lhs.m_index - rhs.m_index;
    }

    friend constexpr auto operator<=>(const Iterator &a,
                                      const Iterator &b) noexcept = default;
    friend constexpr bool operator==(const Iterator &a,
                                     const Iterator &b) noexcept = default;
  };
  using iterator = Iterator;

  explicit IndexRange(index_range_type range) noexcept : m_range(range) {}

  std::size_t size() const noexcept { return m_range.second - m_range.first; }
  bool empty() const noexcept { return size() == 0; }

  iterator begin() const noexcept { return iterator(m_range.first); }
  iterator end() const noexcept { return iterator(m_range.second); }

 private:
  index_range_type m_range{};
};

}  // namespace Acts
