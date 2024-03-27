// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iterator>
#include <type_traits>
#include <utility>

namespace Acts {

namespace detail {

struct Dereference {
  template <typename input_t>
  constexpr static decltype(auto) apply(input_t&& value) {
    return *value;
  }
};

struct ConstDereference {
  template <typename input_t>
  constexpr static decltype(auto) apply(input_t&& value) {
    return std::as_const(*value);
  }
};

struct DotGet {
  template <typename input_t>
  constexpr static decltype(auto) apply(input_t&& value) {
    return value.get();
  }
};

template <typename Callable, typename iterator_t, bool force_const>
struct TransformRangeIterator {
 private:
  using internal_value_type =
      typename std::iterator_traits<iterator_t>::value_type;

  using raw_value_type = std::remove_reference_t<decltype(Callable::apply(
      std::declval<internal_value_type>()))>;

 public:
  using value_type =
      std::conditional_t<force_const, std::add_const_t<raw_value_type>,
                         raw_value_type>;

  using difference_type =
      typename std::iterator_traits<iterator_t>::difference_type;
  using pointer = std::remove_reference_t<value_type>*;
  using reference = value_type&;
  using iterator_category = std::forward_iterator_tag;

  explicit TransformRangeIterator(iterator_t iterator) : m_iterator(iterator) {}

  reference operator*() { return Callable::apply(*m_iterator); }

  reference operator*() const { return Callable::apply(*m_iterator); }

  TransformRangeIterator& operator++() {
    ++m_iterator;
    return *this;
  }

  bool operator==(const TransformRangeIterator& other) const {
    return m_iterator == other.m_iterator;
  }

 private:
  iterator_t m_iterator;
};

template <typename Callable, typename container_t>
struct TransformRange {
 private:
  using internal_value_type = typename container_t::value_type;

  using raw_value_type = std::remove_reference_t<decltype(Callable::apply(
      std::declval<internal_value_type>()))>;

 public:
  using value_type =
      std::conditional_t<std::is_const_v<container_t>,
                         std::add_const_t<raw_value_type>, raw_value_type>;

  using reference = value_type&;
  using const_reference = const value_type&;

  using iterator = TransformRangeIterator<
      Callable,
      std::conditional_t<std::is_const_v<container_t>,
                         typename container_t::const_iterator,
                         typename container_t::iterator>,
      std::is_const_v<container_t>>;
  using const_iterator =
      TransformRangeIterator<Callable, typename container_t::const_iterator,
                             true>;

  explicit TransformRange(Callable&& /*callable*/, container_t& container)
      : m_container(&container) {}

  reference operator[](std::size_t i) {
    return Callable::apply((*m_container)[i]);
  }

  const_reference operator[](std::size_t i) const {
    return std::as_const(Callable::apply((*m_container)[i]));
  }

  reference at(std::size_t i) { return Callable::apply(m_container->at(i)); }

  const_reference at(std::size_t i) const {
    return std::as_const(Callable::apply(m_container->at(i)));
  }

  iterator begin() { return iterator{m_container->begin()}; }

  iterator end() { return iterator{m_container->end()}; }

  const_iterator begin() const { return const_iterator{m_container->begin()}; }

  const_iterator end() const { return const_iterator{m_container->end()}; }

  std::size_t size() const { return m_container->size(); }
  bool empty() const { return m_container->empty(); }

 private:
  container_t* m_container;
};

}  // namespace detail
}  // namespace Acts
