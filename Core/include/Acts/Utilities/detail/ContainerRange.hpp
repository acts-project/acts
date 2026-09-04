// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Utilities/detail/ContainerConcepts.hpp"

#include <concepts>
#include <ranges>
#include <stdexcept>
#include <utility>

namespace Acts::detail {

/// CRTP helper class to represent a subrange of a container defined by a start
/// and end index. The subrange provides an iterator that dereferences to the
/// values in the container corresponding to the indices in the index range. The
/// subrange can be read-only or mutable, depending on the template parameter.
///
/// The user is expected to derive from this class using CRTP and provide the
/// appropriate template parameters. The derived class can then be used to
/// create subsets of the container and access the values through the provided
/// iterator interface.
///
/// @tparam derived_t The derived class type, used for CRTP.
/// @tparam derived_read_only_t The read-only version of the derived class, used for
///         const-correctness.
/// @tparam container_t The type of the underlying container.
/// @tparam index_t The type of the indices that define the subrange.
/// @tparam read_only A boolean indicating whether the subset is read-only.
template <typename derived_t, typename derived_read_only_t,
          typename container_t, std::integral index_t, bool read_only>
class ContainerRange {
 public:
  using Container = const_if_t<read_only, container_t>;
  using Index = index_t;
  static constexpr bool ReadOnly = read_only;
  using IndexRange = std::pair<Index, Index>;

  constexpr ContainerRange(Container &container,
                           const IndexRange &range) noexcept
      : m_container(&container), m_range(range) {}

  template <bool OtherReadOnly>
  explicit constexpr ContainerRange(
      const ContainerRange<derived_t, derived_read_only_t, container_t, Index,
                           OtherReadOnly> &other) noexcept
    requires(ReadOnly && !OtherReadOnly)
      : m_container(&other.container()), m_range(other.range()) {}

  constexpr derived_read_only_t asConst() const noexcept
    requires(!ReadOnly)
  {
    return {container(), range()};
  }

  constexpr Container &container() const noexcept { return *m_container; }

  constexpr const IndexRange &range() const noexcept { return m_range; }

  constexpr std::size_t size() const noexcept {
    return m_range.second - m_range.first;
  }

  constexpr bool empty() const noexcept { return size() == 0; }

  constexpr derived_t subrange(Index offset) const noexcept {
    assert(offset <= m_range.second - m_range.first &&
           "Subrange offset out of bounds");
    return {container(), {m_range.first + offset, m_range.second}};
  }

  constexpr derived_t subrange(Index offset, Index count) const noexcept {
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

  constexpr auto operator[](Index index) const noexcept
    requires(ContainerHasArrayAccess<Container>)
  {
    assert(index < size() && "Index out of bounds");
    return (*m_container)[m_range.first + index];
  }

  constexpr auto at(Index index) const
    requires(ContainerHasAt<Container>)
  {
    if (index >= size()) {
      throw std::out_of_range("Index out of bounds");
    }
    return m_container->at(m_range.first + index);
  }

 private:
  Container *m_container{};
  IndexRange m_range{};
};

}  // namespace Acts::detail
