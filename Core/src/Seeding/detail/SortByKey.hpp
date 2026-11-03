// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <bitset>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <numeric>
#include <span>
#include <utility>

namespace Acts::detail {

/// Sorts at most @p kMaxSize items stably by key with a rank sort: an item's
/// new position is the number of smaller keys, plus the number of equal keys
/// before it.
template <std::size_t kMaxSize, typename item_t, typename key_of_t>
void rankSortByKey(std::span<item_t> items, const key_of_t& keyOf) {
  // A smaller rank sort while the items fit in half, down to two items
  if constexpr (kMaxSize > 2) {
    if (items.size() <= kMaxSize / 2) {
      rankSortByKey<kMaxSize / 2>(items, keyOf);
      return;
    }
  }

  // The type keyOf returns for an item
  using key_t = decltype(keyOf(items.front()));
  const std::size_t n = items.size();

  // The keys, padded to kMaxSize with the largest value: the comparison loop
  // below then has a fixed length, which the compiler can vectorise.
  std::array<key_t, kMaxSize> keys{};
  std::ranges::fill(keys, std::numeric_limits<key_t>::max());
  std::ranges::transform(items, keys.begin(), keyOf);

  // Each item's place is the number of smaller keys. Every place taken is
  // marked as one bit of the mask taken.
  std::array<std::uint32_t, kMaxSize> places{};
  std::bitset<kMaxSize> taken;
  for (std::size_t i = 0; i < n; ++i) {
    const key_t key = keys[i];
    std::uint32_t place = 0;
    for (const key_t other : keys) {
      place += other < key;
    }
    places[i] = place;

    // Mark the place as taken
    taken[place] = true;
  }

  // Equal keys share a place, so fewer places are taken than there are
  // items: only then is each item moved past the equal keys before it.
  if (taken.count() != n) {
    for (std::size_t i = 0; i < n; ++i) {
      const key_t key = keys[i];
      std::uint32_t place = places[i];
      for (std::size_t j = 0; j < i; ++j) {
        place += keys[j] == key;
      }
      places[i] = place;
    }
  }

  // Each item moved to its place, then all of them back in order
  std::array<item_t, kMaxSize> sorted{};
  for (std::size_t i = 0; i < n; ++i) {
    sorted[places[i]] = std::move(items[i]);
  }

  std::ranges::move(std::span(sorted).first(n), items.begin());
}

/// Sorts @p items stably by key, given the smallest and the largest key: they
/// are spread over buckets across the key range, then each bucket is put in
/// order with an insertion pass. A bucket that gets more than
/// @p kQuadraticLinearCrossover items is sorted with std::ranges::stable_sort
/// instead, so that uneven keys cannot make it slow. The keys must lie well
/// within float's range.
template <std::size_t kQuadraticLinearCrossover, typename item_t,
          typename storage_t, typename key_of_t, typename key_t>
void bucketSortByKey(std::span<item_t> items, storage_t& storage,
                     const key_of_t& keyOf, const key_t minKey,
                     const key_t maxKey) {
  if (maxKey == minKey) {
    return;  // all keys equal
  }

  // Two buckets per item keep most of them at one item or none.
  const std::size_t n = items.size();
  const std::size_t nBuckets = 2 * std::bit_ceil(n);

  // The bucket arithmetic is done in float for any key type: the conversion
  // keeps the keys' order, which is all the buckets need.
  const auto lo = static_cast<float>(minKey);
  const auto hi = static_cast<float>(maxKey);

  // A range of about zero, too small to divide by, as when different keys
  // round to the same float, falls back to std::ranges::stable_sort
  const float range = hi - lo;
  if (range <
      static_cast<float>(nBuckets) / std::numeric_limits<float>::max()) {
    std::ranges::stable_sort(items, {}, keyOf);
    return;
  }

  const float scale = static_cast<float>(nBuckets) / range;
  const auto bucketOf = [&](const item_t& item) {
    const auto key = static_cast<float>(keyOf(item));
    return std::min(static_cast<std::size_t>((key - lo) * scale), nBuckets - 1);
  };

  // The storage holds each bucket's first position, then each item's bucket
  storage.buckets.resize(nBuckets + n);
  const std::span<std::uint32_t> starts =
      std::span(storage.buckets).first(nBuckets);
  const std::span<std::uint32_t> bucketOfItem =
      std::span(storage.buckets).subspan(nBuckets, n);

  // Count the items in each bucket, and the most in any one
  std::ranges::fill(starts, 0u);
  std::uint32_t largest = 0;
  for (std::size_t i = 0; i < n; ++i) {
    bucketOfItem[i] = static_cast<std::uint32_t>(bucketOf(items[i]));
    ++starts[bucketOfItem[i]];
    largest = std::max(largest, starts[bucketOfItem[i]]);
  }

  // Turn the counts into each bucket's first position
  std::exclusive_scan(starts.begin(), starts.end(), starts.begin(),
                      std::uint32_t{0});

  // Spread the items over the buckets, in order
  auto& spread = storage.scratch;
  spread.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    spread[starts[bucketOfItem[i]]] = std::move(items[i]);
    ++starts[bucketOfItem[i]];
  }

  // A crowded bucket, as when one key is far from the others, is sorted on its
  // own: the insertion below is quadratic in a bucket's size, and the worst
  // case stays O(n log n).
  if (largest > kQuadraticLinearCrossover) {
    std::uint32_t begin = 0;
    for (const std::uint32_t end : starts) {
      if (end - begin > kQuadraticLinearCrossover) {
        std::ranges::stable_sort(std::span(spread).subspan(begin, end - begin),
                                 {}, keyOf);
      }
      begin = end;
    }
  }

  // Back in order: each item moves only past the larger keys of its bucket.
  for (std::size_t i = 0; i < n; ++i) {
    const key_t key = keyOf(spread[i]);
    std::size_t j = i;
    for (; j > 0 && key < keyOf(items[j - 1]); --j) {
      items[j] = std::move(items[j - 1]);
    }

    items[j] = std::move(spread[i]);
  }
}

/// Sorts @p items stably by the number @p keyOf gives each, of any arithmetic
/// type: the rank sort up to @p kQuadraticLinearCrossover items, the bucket
/// sort above, and std::ranges::stable_sort where a key is not finite.
///
/// @tparam kQuadraticLinearCrossover Up to how many items a quadratic sort is
///         used: the rank sort for the whole list, the insertion pass for one
///         bucket
/// @param items The items to sort
/// @param storage Kept between calls: `scratch`, a vector of items, and
///        `buckets`, a vector of std::uint32_t, both of any content
/// @param keyOf Returns the key of an item
template <std::size_t kQuadraticLinearCrossover, typename item_t,
          typename storage_t, typename key_of_t>
void sortByKey(std::span<item_t> items, storage_t& storage,
               const key_of_t& keyOf) {
  // The type keyOf returns for an item
  using key_t = decltype(keyOf(items.front()));

  // Find the key range for the bucket sort. A list with a NaN or infinite key
  // is sorted by std::ranges::stable_sort instead; isfinite is checked before
  // the key is compared here, as comparing a NaN traps where floating point
  // exceptions are enabled.
  bool finite = true;
  key_t minKey = std::numeric_limits<key_t>::max();
  key_t maxKey = std::numeric_limits<key_t>::lowest();

  for (const item_t& item : items) {
    const key_t key = keyOf(item);
    if (std::isfinite(key)) {
      minKey = std::min(minKey, key);
      maxKey = std::max(maxKey, key);
    } else {
      finite = false;
    }
  }

  if (!finite) {
    std::ranges::stable_sort(items, {}, keyOf);
    return;
  }

  if (items.size() <= kQuadraticLinearCrossover) {
    rankSortByKey<kQuadraticLinearCrossover>(items, keyOf);
  } else {
    bucketSortByKey<kQuadraticLinearCrossover>(items, storage, keyOf, minKey,
                                               maxKey);
  }
}

/// Sorts @p items stably by key, and the items with equal keys by
/// @p orderOfEqualKeys. With float keys equal keys are rare, so this costs
/// little more than the sort by key.
///
/// @tparam kQuadraticLinearCrossover Up to how many items a quadratic sort is
///         used: the rank sort for the whole list, the insertion pass for one
///         bucket
/// @param items The items to sort
/// @param storage Kept between calls: `scratch`, a vector of items, and
///        `buckets`, a vector of std::uint32_t, both of any content
/// @param keyOf Returns the key of an item
/// @param orderOfEqualKeys Whether the first of two items with equal keys goes
///        before the second, a strict weak order; it is called only on items
///        with equal keys
template <std::size_t kQuadraticLinearCrossover, typename item_t,
          typename storage_t, typename key_of_t, typename order_t>
void sortByKey(std::span<item_t> items, storage_t& storage,
               const key_of_t& keyOf, const order_t& orderOfEqualKeys) {
  sortByKey<kQuadraticLinearCrossover>(items, storage, keyOf);

  for (auto first = items.begin(); first != items.end();) {
    const auto key = keyOf(*first);
    const auto last =
        std::find_if(std::next(first), items.end(),
                     [&](const item_t& item) { return keyOf(item) != key; });
    if (std::distance(first, last) > 1) {
      std::stable_sort(first, last, orderOfEqualKeys);
    }
    first = last;
  }
}

}  // namespace Acts::detail
