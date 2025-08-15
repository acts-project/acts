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

#include <cstddef>
#include <iterator>

namespace Acts::detail {

template <typename Container, typename Value, typename Index, bool ReadOnly>
class ContainerIterator {
 public:
  using container_type = const_if_t<ReadOnly, Container>;
  using index_type = Index;
  static constexpr bool read_only = ReadOnly;

  using value_type = Value;
  using difference_type = std::ptrdiff_t;
  using pointer = void;
  using reference = void;

  using iterator_category = std::random_access_iterator_tag;
  using iterator_concept = std::random_access_iterator_tag;

  constexpr ContainerIterator() noexcept = default;
  constexpr ContainerIterator(container_type &container,
                              index_type index) noexcept
      : m_container(&container), m_index(index) {}
  template <bool OtherReadOnly>
  explicit constexpr ContainerIterator(
      const ContainerIterator<Container, Value, Index, OtherReadOnly>
          &other) noexcept
    requires(ReadOnly && !OtherReadOnly)
      : m_container(&other.container()), m_index(other.index()) {}

  constexpr ContainerIterator<Container, Value, Index, true> asConst()
      const noexcept
    requires(!ReadOnly)
  {
    return {*m_container, m_index};
  }

  constexpr container_type &container() const noexcept { return *m_container; }
  constexpr index_type index() const noexcept { return m_index; }

  constexpr value_type operator*() const {
    static_assert(
        ContainerHasAt<Container> || ContainerHasArrayAccess<Container>,
        "Container must support at() or operator[] for indexing");
    constexpr bool HasArrayAccess = ContainerHasArrayAccess<Container>;

    if constexpr (HasArrayAccess) {
      return (*m_container)[m_index];
    } else {
      return m_container->at(m_index);
    }
  }
  constexpr value_type operator[](difference_type n) const {
    static_assert(
        ContainerHasAt<Container> || ContainerHasArrayAccess<Container>,
        "Container must support at() or operator[] for indexing");
    constexpr bool HasArrayAccess = ContainerHasArrayAccess<Container>;

    if constexpr (HasArrayAccess) {
      return (*m_container)[m_index + n];
    } else {
      return m_container->at(m_index + n);
    }
  }

  constexpr ContainerIterator &operator++() noexcept {
    ++m_index;
    return *this;
  }
  constexpr ContainerIterator operator++(int) noexcept {
    auto tmp = *this;
    ++(*this);
    return tmp;
  }
  constexpr ContainerIterator &operator--() noexcept {
    --m_index;
    return *this;
  }
  constexpr ContainerIterator operator--(int) noexcept {
    auto tmp = *this;
    --(*this);
    return tmp;
  }

  constexpr ContainerIterator &operator+=(difference_type n) noexcept {
    m_index += n;
    return *this;
  }
  constexpr ContainerIterator &operator-=(difference_type n) noexcept {
    m_index -= n;
    return *this;
  }

 private:
  container_type *m_container{};
  index_type m_index{};

  friend constexpr ContainerIterator operator+(ContainerIterator it,
                                               difference_type n) noexcept {
    return it += n;
  }

  friend constexpr ContainerIterator operator+(difference_type n,
                                               ContainerIterator it) noexcept {
    return it += n;
  }

  friend constexpr ContainerIterator operator-(ContainerIterator it,
                                               difference_type n) noexcept {
    return it -= n;
  }

  friend constexpr difference_type operator-(
      const ContainerIterator &lhs, const ContainerIterator &rhs) noexcept {
    return lhs.m_index - rhs.m_index;
  }

  friend constexpr auto operator<=>(const ContainerIterator &a,
                                    const ContainerIterator &b) noexcept {
    return a.m_index <=> b.m_index;
  }
  friend constexpr bool operator==(const ContainerIterator &a,
                                   const ContainerIterator &b) noexcept {
    return a.m_index == b.m_index;
  }
};

}  // namespace Acts::detail
