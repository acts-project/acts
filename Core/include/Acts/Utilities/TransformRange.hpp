// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <type_traits>
#include <utility>

namespace Acts {

namespace detail {

template <typename input_t>
struct Dereference {
  using value_type = decltype(*std::declval<input_t>());
  constexpr static value_type apply(input_t& value) { return *value; }
  constexpr static value_type apply(const input_t& value) { return *value; }
};

template <typename input_t>
struct ConstDereference {
  using value_type = decltype(std::as_const(*std::declval<input_t>()));

  constexpr static value_type apply(input_t& value) { return *value; }
  constexpr static value_type apply(const input_t& value) { return *value; }
};

template <template <typename> typename Callable, typename iterator_t,
          bool force_const>
struct TransformRangeIterator {
 private:
  using internal_value_type = typename iterator_t::value_type;

 public:
  using value_type = std::conditional_t<
      force_const,
      decltype(std::as_const(
          std::declval<typename Callable<internal_value_type>::value_type>())),
      typename Callable<internal_value_type>::value_type>;

  explicit TransformRangeIterator(iterator_t iterator) : m_iterator(iterator) {}

  value_type operator*() {
    return Callable<typename iterator_t::value_type>::apply(*m_iterator);
  }

  value_type operator*() const {
    return Callable<typename iterator_t::value_type>::apply(*m_iterator);
  }

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

template <template <typename> typename Callable, typename container_t>
struct TransformRange {
 private:
  using internal_value_type = typename container_t::value_type;

 public:
  using value_type = std::conditional_t<
      std::is_const_v<container_t>,
      decltype(std::as_const(
          std::declval<typename Callable<internal_value_type>::value_type>())),
      typename Callable<internal_value_type>::value_type>;

  using iterator = TransformRangeIterator<
      Callable,
      std::conditional_t<std::is_const_v<container_t>,
                         typename container_t::const_iterator,
                         typename container_t::iterator>,
      std::is_const_v<container_t>>;
  using const_iterator =
      TransformRangeIterator<Callable, typename container_t::const_iterator,
                             true>;

  explicit TransformRange(container_t& container) : m_container(&container) {}

  value_type operator[](std::size_t i) {
    return Callable<internal_value_type>::apply((*m_container)[i]);
  }

  decltype(auto) operator[](std::size_t i) const {
    return std::as_const(
        Callable<internal_value_type>::apply((*m_container)[i]));
  }

  value_type at(std::size_t i) {
    return Callable<internal_value_type>::apply(m_container->at(i));
  }

  decltype(auto) at(std::size_t i) const {
    return std::as_const(
        Callable<internal_value_type>::apply(m_container->at(i)));
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

template <template <typename> typename Callable, typename container_t>
TransformRange<Callable, container_t> make_transform_range(
    container_t& container) {
  return TransformRange<Callable, container_t>(container);
}

}  // namespace detail
}  // namespace Acts
