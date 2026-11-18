// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/indexing.hpp"
#include "detray/utils/hash_combine.hpp"

// System include(s)
#include <algorithm>
#include <concepts>
#include <cstddef>
#include <functional>
#include <typeindex>
#include <unordered_map>
#include <vector>

namespace detray::detail {

/// @brief Hash based lookup of identical objects in detector data collections
///
/// The cache holds one hash map per detector collection and hash function
/// type. It is filled lazily: Before every lookup, the cache hashes the
/// objects that were added to the collection since the last lookup. This way
/// the cache is always consistent with the collection, also if objects were
/// added without the cache (e.g. when the deduplication was switched off).
///
/// The hash function only selects the candidates. The equality operator of
/// the object type decides if an object is identical. The hash function must
/// therefore be consistent with the equality operator: Identical objects must
/// have the same hash value.
///
/// @note The cache always accesses the collections with the default
/// geometry context, since the detector is always built in this context.
/// @note The cache identifies a collection by its address. The collections
/// must not move while the cache is in use (e.g. during a detector build).
/// @note The collections must only grow while the cache is in use. If a
/// collection gets smaller, the cache for this collection is reset. The
/// owner of the cache (e.g. the detector builder) has to clear it, before the
/// collections are destroyed or replaced.
class deduplication_cache {
 public:
  /// Find all objects in @p coll that are identical to @p obj
  ///
  /// Use this to find blocks of objects: A collection can contain identical
  /// objects, if they were added as part of different blocks.
  ///
  /// @tparam coll_t the type of the collection
  /// @tparam obj_t the type of the object
  /// @tparam hasher_t the type of the hash function for the objects
  ///
  /// @param coll the collection that is searched
  /// @param obj the object to be found
  /// @param hasher the hash function for the objects
  ///
  /// @returns the indices of the identical objects in ascending order
  template <typename coll_t, typename obj_t, typename hasher_t>
  auto find_all(const coll_t &coll, const obj_t &obj, const hasher_t &hasher)
      -> std::vector<dindex> {
    const auto &lookup = update(coll, hasher).m_lookup;

    std::vector<dindex> result{};
    const auto [first, last] = lookup.equal_range(hasher(obj));
    for (auto itr = first; itr != last; ++itr) {
      const dindex idx{itr->second};
      if (obj == coll.at(idx)) {
        result.push_back(idx);
      }
    }
    std::ranges::sort(result);

    return result;
  }

  /// Find an object in @p coll that is identical to @p obj
  ///
  /// Use this for collections in which every object is added through the
  /// cache one by one: Such a collection holds at most one object that is
  /// identical to @p obj.
  ///
  /// @tparam coll_t the type of the collection
  /// @tparam obj_t the type of the object
  /// @tparam hasher_t the type of the hash function for the objects
  ///
  /// @param coll the collection that is searched
  /// @param obj the object to be found
  /// @param hasher the hash function for the objects
  ///
  /// @returns the index of an identical object or @c dindex_invalid
  template <typename coll_t, typename obj_t, typename hasher_t>
  auto find(const coll_t &coll, const obj_t &obj, const hasher_t &hasher)
      -> dindex {
    const auto &lookup = update(coll, hasher).m_lookup;

    // The hash only selects the candidates (hash collisions are possible)
    const auto [first, last] = lookup.equal_range(hasher(obj));
    for (auto itr = first; itr != last; ++itr) {
      if (obj == coll.at(itr->second)) {
        return itr->second;
      }
    }

    return dindex_invalid;
  }

  /// @returns the number of collections that the cache holds a map for
  std::size_t size() const { return m_entries.size(); }

  /// Remove all hash maps
  void clear() { m_entries.clear(); }

 private:
  /// Identifies a hash map: Collection address and hash function type
  struct key {
    const void *m_coll{nullptr};
    std::type_index m_hasher;

    bool operator==(const key &) const = default;
  };

  /// Hash function for the hash map keys
  struct key_hash {
    std::size_t operator()(const key &k) const {
      return hash_combine(std::hash<const void *>{}(k.m_coll), k.m_hasher);
    }
  };

  /// Hash map of one collection
  struct entry {
    /// Number of collection objects that are in the hash map
    std::size_t m_n_hashed{0u};
    /// Map from the object hash to the object index in the collection
    std::unordered_multimap<std::size_t, dindex> m_lookup{};
  };

  /// Add the objects of @p coll that are not yet in the hash map
  ///
  /// @returns the up-to-date hash map entry of the collection
  template <typename coll_t, typename hasher_t>
  auto update(const coll_t &coll, const hasher_t &hasher) -> entry & {
    entry &e = m_entries[key{&coll, std::type_index(typeid(hasher_t))}];

    const auto n_objs{static_cast<std::size_t>(coll.size())};
    // The collection was reset: Start again
    if (n_objs < e.m_n_hashed) {
      e = entry{};
    }
    for (std::size_t i = e.m_n_hashed; i < n_objs; ++i) {
      const auto idx{static_cast<dindex>(i)};
      e.m_lookup.emplace(hasher(coll.at(idx)), idx);
    }
    e.m_n_hashed = n_objs;

    return e;
  }

  /// The hash maps of all collections
  std::unordered_map<key, entry, key_hash> m_entries{};
};

}  // namespace detray::detail
