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
#include <ranges>
#include <stdexcept>

namespace Acts::detail {

/// CRTP helper class to represent a subset of a container defined by an index
/// range. The subset provides an iterator that dereferences to the values in
/// the container corresponding to the indices in the index range. The subset
/// can be read-only or mutable, depending on the template parameter.
///
/// The user is expected to derive from this class using CRTP and provide the
/// appropriate template parameters. The derived class can then be used to
/// create subsets of the container and access the values through the provided
/// iterator interface.
///
/// @tparam Derived The derived class type, used for CRTP.
/// @tparam DerivedReadOnly The read-only version of the derived class, used for
///         const-correctness.
/// @tparam Container The type of the underlying container.
/// @tparam Value The type of the values in the container.
/// @tparam IndexRange The type of the index range that defines the subset.
/// @tparam ReadOnly A boolean indicating whether the subset is read-only.
template <typename Derived, typename DerivedReadOnly, typename Container,
          typename Value, std::ranges::range IndexRange, bool ReadOnly>
class ContainerSubset {
 public:
  using container_type = const_if_t<ReadOnly, Container>;
  using value_type = Value;
  static constexpr bool read_only = ReadOnly;
  using index_range_type = IndexRange;
  using index_type = typename index_range_type::value_type;

  template <typename subset_iterator>
  class Iterator {
   private:
    static_assert(std::input_iterator<subset_iterator>,
                  "Subset iterator must satisfy input iterator requirements");
    static constexpr bool IsBidirectional =
        std::bidirectional_iterator<subset_iterator>;
    static constexpr bool IsRandomAccess =
        std::random_access_iterator<subset_iterator>;

   public:
    using value_type = Value;
    using difference_type = index_range_type::difference_type;
    using pointer = void;
    using reference = void;

    using iterator_category = std::conditional_t<
        IsRandomAccess, std::random_access_iterator_tag,
        std::conditional_t<IsBidirectional, std::bidirectional_iterator_tag,
                           std::input_iterator_tag>>;
    using iterator_concept = iterator_category;

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

    constexpr value_type operator[](difference_type n) const
      requires(IsRandomAccess)
    {
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

    constexpr Iterator &operator--() noexcept
      requires(IsBidirectional)
    {
      --m_iterator;
      return *this;
    }

    constexpr Iterator operator--(int) noexcept
      requires(IsBidirectional)
    {
      auto tmp = *this;
      --(*this);
      return tmp;
    }

    constexpr Iterator &operator+=(difference_type n) noexcept
      requires(IsRandomAccess)
    {
      m_iterator += n;
      return *this;
    }

    constexpr Iterator &operator-=(difference_type n) noexcept
      requires(IsRandomAccess)
    {
      m_iterator -= n;
      return *this;
    }

   private:
    container_type *m_container{};
    subset_iterator m_iterator{};

    friend constexpr Iterator operator+(Iterator it, difference_type n) noexcept
      requires(IsRandomAccess)
    {
      return it += n;
    }

    friend constexpr Iterator operator+(difference_type n, Iterator it) noexcept
      requires(IsRandomAccess)
    {
      return it += n;
    }

    friend constexpr Iterator operator-(Iterator it, difference_type n) noexcept
      requires(IsRandomAccess)
    {
      return it -= n;
    }

    friend constexpr difference_type operator-(const Iterator &lhs,
                                               const Iterator &rhs) noexcept
      requires(IsRandomAccess)
    {
      return lhs.m_iterator - rhs.m_iterator;
    }

    friend constexpr auto operator<=>(const Iterator &a,
                                      const Iterator &b) noexcept
      requires(IsRandomAccess)
    {
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
  using iterator =
      Iterator<decltype(std::declval<const index_range_type>().begin())>;

  constexpr ContainerSubset(container_type &container,
                            index_range_type subset) noexcept
      : m_container(&container), m_subset(std::move(subset)) {}
  template <bool OtherReadOnly>
  explicit constexpr ContainerSubset(
      const ContainerSubset<Derived, DerivedReadOnly, Container, Value,
                            IndexRange, OtherReadOnly> &other) noexcept
    requires(read_only && !OtherReadOnly)
      : m_container(&other.container()), m_subset(other.subset()) {}

  constexpr DerivedReadOnly asConst() const noexcept
    requires(!read_only)
  {
    return {*m_container, m_subset};
  }

  constexpr container_type &container() const noexcept { return *m_container; }
  constexpr const index_range_type &subset() const noexcept { return m_subset; }

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

  constexpr auto begin() const noexcept {
    return Iterator(container(), subset().begin());
  }
  constexpr auto end() const noexcept {
    return Iterator(container(), subset().end());
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
  index_range_type m_subset{};
};

}  // namespace Acts::detail
