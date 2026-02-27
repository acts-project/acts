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

#include <compare>
#include <iterator>
#include <stdexcept>

namespace Acts::detail {

template <typename Derived, typename DerivedReadOnly, typename Container,
          typename Value, typename IndexContainer, bool ReadOnly>
class ContainerSubset {
 public:
  using container_type = const_if_t<ReadOnly, Container>;
  using value_type = Value;
  static constexpr bool read_only = ReadOnly;
  using index_container_type = IndexContainer;
  using index_type = typename index_container_type::value_type;

  class Iterator {
   public:
    using subset_iterator = index_container_type::iterator;

    using value_type = Value;
    using difference_type = std::ptrdiff_t;
    using pointer = void;
    using reference = void;

    using iterator_category = std::random_access_iterator_tag;
    using iterator_concept = std::random_access_iterator_tag;

    constexpr Iterator() noexcept = default;
    constexpr Iterator(container_type &container,
                       subset_iterator iterator) noexcept
        : m_container(&container), m_iterator(iterator) {}

    constexpr value_type operator*() const {
      static_assert(
          ContainerHasAt<Container> || ContainerHasArrayAccess<Container>,
          "Container must support at() or operator[] for indexing");
      constexpr bool HasArrayAccess = ContainerHasArrayAccess<Container>;

      if constexpr (HasArrayAccess) {
        return (*m_container)[*m_iterator];
      } else {
        return m_container->at(*m_iterator);
      }
    }
    constexpr value_type operator[](difference_type n) const {
      static_assert(
          ContainerHasAt<Container> || ContainerHasArrayAccess<Container>,
          "Container must support at() or operator[] for indexing");
      constexpr bool HasArrayAccess = ContainerHasArrayAccess<Container>;

      if constexpr (HasArrayAccess) {
        return (*m_container)[m_iterator[n]];
      } else {
        return m_container->at(m_iterator[n]);
      }
    }

    constexpr Iterator &operator++() noexcept {
      ++m_iterator;
      return *this;
    }
    constexpr Iterator operator++(int) noexcept {
      auto tmp = *this;
      ++(*this);
      return tmp;
    }
    constexpr Iterator &operator--() noexcept {
      --m_iterator;
      return *this;
    }
    constexpr Iterator operator--(int) noexcept {
      auto tmp = *this;
      --(*this);
      return tmp;
    }

    constexpr Iterator &operator+=(difference_type n) noexcept {
      m_iterator += n;
      return *this;
    }
    constexpr Iterator &operator-=(difference_type n) noexcept {
      m_iterator -= n;
      return *this;
    }

   private:
    container_type *m_container{};
    subset_iterator m_iterator{};

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
      return lhs.m_iterator - rhs.m_iterator;
    }

    friend constexpr auto operator<=>(const Iterator &a,
                                      const Iterator &b) noexcept {
      if (a.m_iterator < b.m_iterator) {
        return std::strong_ordering::less;
      }
      if (a.m_iterator > b.m_iterator) {
        return std::strong_ordering::greater;
      }
      return std::strong_ordering::equal;
    }
    friend constexpr bool operator==(const Iterator &a,
                                     const Iterator &b) noexcept {
      return a.m_iterator == b.m_iterator;
    }
  };
  using iterator = Iterator;

  constexpr ContainerSubset(container_type &container,
                            index_container_type subset) noexcept
      : m_container(&container), m_subset(std::move(subset)) {}
  template <bool OtherReadOnly>
  explicit constexpr ContainerSubset(
      const ContainerSubset<Derived, DerivedReadOnly, Container, Value,
                            IndexContainer, OtherReadOnly> &other) noexcept
    requires(ReadOnly && !OtherReadOnly)
      : m_container(&other.container()), m_subset(other.subset()) {}

  constexpr DerivedReadOnly asConst() const noexcept
    requires(!ReadOnly)
  {
    return {*m_container, m_subset};
  }

  constexpr container_type &container() const noexcept { return *m_container; }
  constexpr const index_container_type &subset() const noexcept {
    return m_subset;
  }

  constexpr std::size_t size() const noexcept { return subset().size(); }
  constexpr bool empty() const noexcept { return size() == 0; }

  constexpr Derived subrange(std::size_t offset) const noexcept {
    return {container(), subset().subspan(offset)};
  }
  constexpr Derived subrange(std::size_t offset,
                             std::size_t count) const noexcept {
    return {container(), subset().subspan(offset, count)};
  }

  constexpr auto front() const noexcept {
    return container()[subset().front()];
  }
  constexpr auto back() const noexcept { return container()[subset().back()]; }

  constexpr iterator begin() const noexcept {
    return iterator(container(), subset().begin());
  }
  constexpr iterator end() const noexcept {
    return iterator(container(), subset().end());
  }

  constexpr auto operator[](index_type index) const noexcept
    requires(ContainerHasArrayAccess<Container>)
  {
    assert(index < size() && "Index out of bounds");
    return container()[subset()[index]];
  }
  constexpr auto at(index_type index) const
    requires(ContainerHasAt<Container>)
  {
    if (index >= size()) {
      throw std::out_of_range("Index out of bounds");
    }
    return container().at(subset()[index]);
  }

 private:
  container_type *m_container{};
  index_container_type m_subset{};
};

}  // namespace Acts::detail
