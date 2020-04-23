// SPDX-License-Identifier: MIT
// Copyright 2019 Moritz Kiehn
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

/// \file
/// \brief   Minimal associative containers based on continous storage
/// \author  Moritz Kiehn <msmk@cern.ch>
/// \date    2019-02-27, Initial version

#pragma once

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <vector>

namespace dfe {

/// An container adaptor to store a set of elements in a sequential container.
///
/// \tparam T         Stored element type
/// \tparam Compare   Function satisfying the `Compare` named requirement
/// \tparam Container Sequential container
///
/// Supports access by equivalence, iteration over all elements, removing all
/// elements, adding elements while maintaining uniqueness, and membership
/// checks. By using a sequential container, memory allocation is greatly
/// simplified and lookups benefit from higher memory locality at the expense of
/// slower insertion of elements. Should work best for smaller sets with
/// frequent lookups.
///
/// The set elements can not be modified on purpose. With a non-standard
/// `Compare` function modifying a contained object might change its identity
/// and thus its position in the set. This would break the internal sorting.
template<
  typename T, typename Compare = std::less<T>,
  typename Container = std::vector<T>>
class FlatSet {
public:
  using value_type = T;
  using size_type = typename Container::size_type;
  using const_iterator = typename Container::const_iterator;

  /// Access the equivalent element or throw if it does not exists.
  template<typename U>
  const value_type& at(U&& u) const;

  const_iterator begin() const { return m_items.begin(); }
  const_iterator end() const { return m_items.end(); }

  /// Return true if there are no elements in the set.
  bool empty() const { return m_items.empty(); }
  /// Return the number of elements in the set.
  size_type size() const { return m_items.size(); }

  /// Remove all elements from the container.
  void clear() { m_items.clear(); }
  /// Add the element to the set or replace the existing equivalent element.
  ///
  /// Depending on the `Compare` function we might have elements with different
  /// values but the same identity with respect to the chosen `Compare`
  /// function. Only one can be kept and this function replaces the existing
  /// element with the new one in such a case.
  void insert_or_assign(const T& t);

  /// Return an interator to the equivalent element or `.end()` if not found.
  template<typename U>
  const_iterator find(U&& u) const;
  /// Return true if the equivalent element is in the set.
  template<typename U>
  bool contains(U&& u) const;

private:
  Container m_items;
};

/// A key-value map that stores keys and values in sequential containers.
///
/// \tparam Key     Stored element key type
/// \tparam T       Stored element value type
/// \tparam Compare Function satisfying the `Compare` name requirements for keys
///
/// Supports access by key, clearing all elements, adding or replacing the
/// stored value for a given key, and membership checks. Keys and values are
/// stored in separate sequential containers to simplify allocation and benefit
/// from greater memory locality.
template<typename Key, typename T, typename Compare = std::less<Key>>
class FlatMap {
public:
  using key_type = Key;
  using value_type = T;
  using size_type = std::size_t;

  /// Writable access to an element or throw if it does not exists.
  value_type& at(const Key& key) { return m_items[m_keys.at(key).index]; }
  /// Read-only access to an element or throw if it does not exists.
  const value_type& at(const Key& key) const {
    return m_items[m_keys.at(key).index];
  }

  /// Return true if there are no elements in the map.
  bool empty() const { return m_keys.empty(); }
  /// Return the number of elements in the container.
  size_type size() const { return m_keys.size(); }

  /// Remove all elements from the container.
  void clear() { m_keys.clear(), m_items.clear(); }
  /// Add the element under the given key or replace an existing element.
  ///
  /// New elements are constructed or assigned in-place with the parameters
  /// forwarded to a `T(...)` constructor call.
  template<typename... Params>
  void emplace(const Key& key, Params&&... params);

  /// Return true if an element exists for the given key
  bool contains(const Key& key) const { return m_keys.contains(key); }

private:
  struct KeyIndex {
    Key key;
    size_type index;
  };
  struct KeyCompare {
    constexpr bool operator()(const KeyIndex& lhs, const KeyIndex& rhs) const {
      return Compare()(lhs.key, rhs.key);
    }
    constexpr bool operator()(const KeyIndex& lhs, const Key& rhs_key) const {
      return Compare()(lhs.key, rhs_key);
    }
    constexpr bool operator()(const Key& lhs_key, const KeyIndex& rhs) const {
      return Compare()(lhs_key, rhs.key);
    }
  };

  FlatSet<KeyIndex, KeyCompare> m_keys;
  std::vector<T> m_items;
};

// implementation FlatSet

template<typename T, typename Compare, typename Container>
template<typename U>
inline const typename FlatSet<T, Compare, Container>::value_type&
FlatSet<T, Compare, Container>::at(U&& u) const {
  auto pos = find(std::forward<U>(u));
  if (pos == end()) {
    throw std::out_of_range("The requested element does not exists");
  }
  return *pos;
}

template<typename T, typename Compare, typename Container>
inline void
FlatSet<T, Compare, Container>::insert_or_assign(const T& t) {
  auto pos = std::lower_bound(m_items.begin(), m_items.end(), t, Compare());
  if (((pos != m_items.end()) and !Compare()(t, *pos))) {
    *pos = t;
  } else {
    m_items.emplace(pos, t);
  }
}

template<typename T, typename Compare, typename Container>
template<typename U>
inline typename FlatSet<T, Compare, Container>::const_iterator
FlatSet<T, Compare, Container>::find(U&& u) const {
  auto end = m_items.end();
  auto pos =
    std::lower_bound(m_items.begin(), end, std::forward<U>(u), Compare());
  return ((pos != end) and !Compare()(std::forward<U>(u), *pos)) ? pos : end;
}

template<typename T, typename Compare, typename Container>
template<typename U>
inline bool
FlatSet<T, Compare, Container>::contains(U&& u) const {
  return std::binary_search(
    m_items.begin(), m_items.end(), std::forward<U>(u), Compare());
}

// implementation FlatMap

template<typename Key, typename T, typename Compare>
template<typename... Params>
inline void
FlatMap<Key, T, Compare>::emplace(const Key& key, Params&&... params) {
  auto idx = m_keys.find(key);
  if (idx != m_keys.end()) {
    m_items[idx->index] = T(std::forward<Params>(params)...);
  } else {
    m_items.emplace_back(std::forward<Params>(params)...);
    m_keys.insert_or_assign(KeyIndex{key, m_items.size() - 1});
  }
}

} // namespace dfe
