// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <functional>
#include <list>
#include <optional>
#include <stdexcept>
#include <unordered_map>

namespace Acts::detail {

/// @brief Fixed-capacity least-recently-used (LRU) cache.
///
/// Associates keys of type @p Key with values of type @p Value and keeps at
/// most @p capacity entries in memory. When the cache is full and a new entry
/// is inserted the least-recently-used entry is evicted.
///
/// All operations — @c get, @c put — are O(1). "Age" is structural: entries
/// are stored in a doubly-linked list ordered MRU → LRU; a cache hit splices
/// the found node to the front in O(1) without invalidating any iterators.
///
/// Not thread-safe. Callers that need concurrent access must hold an external
/// mutex around every @c get / @c put call.
///
/// @tparam Key   Key type. Must be hashable with @p Hash and comparable with
///               @p Equal.
/// @tparam Value Value type. Must be copy- or move-constructible.
/// @tparam Hash  Hash functor for @p Key (default: @c std::hash<Key>).
/// @tparam Equal Equality functor for @p Key (default: @c std::equal_to<Key>).
template <typename Key, typename Value, typename Hash = std::hash<Key>,
          typename Equal = std::equal_to<Key>>
class LruCache {
 public:
  /// @brief Construct a cache with the given capacity.
  /// @param capacity Maximum number of entries to keep.
  /// @throws std::invalid_argument if @p capacity is zero.
  explicit LruCache(std::size_t capacity) : m_capacity(capacity) {
    if (capacity == 0) {
      throw std::invalid_argument("LruCache: capacity must be >= 1");
    }
  }

  /// @brief Look up @p key; promote it to MRU on hit.
  /// @param key Key to look up.
  /// @return A copy of the cached value, or @c std::nullopt on miss.
  std::optional<Value> get(const Key& key) {
    auto it = m_index.find(key);
    if (it == m_index.end()) {
      return std::nullopt;
    }
    // Splice to front — O(1), does not invalidate other iterators.
    m_list.splice(m_list.begin(), m_list, it->second);
    return it->second->value;
  }

  /// @brief Insert or update @p key → @p value; evict LRU when at capacity.
  /// @param key Key to insert or update.
  /// @param value Value to associate with @p key.
  void put(const Key& key, Value value) {
    auto it = m_index.find(key);
    if (it != m_index.end()) {
      it->second->value = std::move(value);
      m_list.splice(m_list.begin(), m_list, it->second);
      return;
    }
    m_list.push_front(Entry{key, std::move(value)});
    m_index.emplace(key, m_list.begin());
    if (m_list.size() > m_capacity) {
      m_index.erase(m_list.back().key);
      m_list.pop_back();
    }
  }

  /// @return Number of entries currently in the cache.
  std::size_t size() const { return m_list.size(); }

  /// @return Maximum number of entries the cache will hold.
  std::size_t capacity() const { return m_capacity; }

 private:
  struct Entry {
    Key key;
    Value value;
  };

  std::size_t m_capacity;
  std::list<Entry> m_list;  // front = MRU, back = LRU
  std::unordered_map<Key, typename std::list<Entry>::iterator, Hash, Equal>
      m_index;
};

}  // namespace Acts::detail
