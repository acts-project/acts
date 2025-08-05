// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/TypeTraits.hpp"

#include <utility>

namespace Acts {

template <typename Derived, typename DerivedReadOnly, typename Container,
          typename Index, bool ReadOnly>
class ContainerRange {
 public:
  using container_type = const_if_t<ReadOnly, Container>;
  using index_type = Index;
  static constexpr bool read_only = ReadOnly;
  using index_range_type = std::pair<index_type, index_type>;

  constexpr ContainerRange(container_type &container,
                           const index_range_type &range) noexcept
      : m_container(&container), m_range(range) {}
  template <bool OtherReadOnly>
  explicit constexpr ContainerRange(
      const ContainerRange<Derived, DerivedReadOnly, Container, Index,
                           OtherReadOnly> &other) noexcept
    requires(ReadOnly && !OtherReadOnly)
      : m_container(&other.container()), m_range(other.range()) {}

  constexpr DerivedReadOnly asConst() const noexcept
    requires(!ReadOnly)
  {
    return {container(), range()};
  }

  constexpr container_type &container() const noexcept { return *m_container; }
  constexpr const index_range_type &range() const noexcept { return m_range; }

  constexpr std::size_t size() const noexcept {
    return m_range.second - m_range.first;
  }
  constexpr bool empty() const noexcept { return size() == 0; }

  constexpr Derived subrange(index_type offset) const noexcept {
    assert(offset <= m_range.second - m_range.first &&
           "Subrange offset out of bounds");
    return {container(), {m_range.first + offset, m_range.second}};
  }
  constexpr Derived subrange(index_type offset,
                             index_type count) const noexcept {
    assert(offset <= m_range.second - m_range.first &&
           "Subrange offset out of bounds");
    assert(count <= m_range.second - m_range.first - offset &&
           "Subrange count out of bounds");
    return {container(),
            {m_range.first + offset, m_range.first + offset + count}};
  }

  constexpr auto front() const noexcept { return container()[m_range.first]; }
  constexpr auto back() const noexcept {
    return container()[m_range.second - 1];
  }

  constexpr auto begin() const noexcept {
    return container().begin() + m_range.first;
  }
  constexpr auto end() const noexcept {
    return container().begin() + m_range.second;
  }

 private:
  container_type *m_container{};
  index_range_type m_range{};
};

}  // namespace Acts
