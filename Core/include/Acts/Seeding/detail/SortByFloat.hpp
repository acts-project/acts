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
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace Acts::detail {

/// Sorts @p items stably by the float key @p keyOf gives each: each item is
/// moved back past the larger keys before it. For a few items, where it is
/// faster than the other sorts here.
///
/// @param items The items to sort
/// @param keyOf Returns the float key of an item
template <typename item_t, typename key_of_t>
void insertionSortByFloat(std::vector<item_t>& items, key_of_t keyOf) {
  for (std::size_t i = 1; i < items.size(); ++i) {
    item_t item = std::move(items[i]);
    const float key = keyOf(item);
    std::size_t j = i;
    for (; j > 0 && keyOf(items[j - 1]) > key; --j) {
      items[j] = std::move(items[j - 1]);
    }
    items[j] = std::move(item);
  }
}

/// Sorts at most @p MaxSize @p items stably by the float key @p keyOf gives
/// each: an item's place is the number of smaller keys, plus, only if two
/// items share a place, the number of equal keys before it. Quadratic, but
/// without a branch, so faster than a comparison sort for a few tens of
/// items. The keys are compared over all @p MaxSize slots, the unused ones
/// holding +inf: a loop length known at compile time lets compilers
/// vectorise it at -O2 too.
///
/// Every key must be finite, and @p item_t default-constructible.
///
/// @param items The items to sort, at most @p MaxSize of them
/// @param keyOf Returns the float key of an item
template <std::size_t MaxSize, typename item_t, typename key_of_t>
void rankSortByFloat(std::vector<item_t>& items, key_of_t keyOf) {
  static_assert(MaxSize <= 64, "the places must fit in a 64-bit mask");
  const std::size_t n = items.size();
  assert(n <= MaxSize && "too many items for a rank sort");
  // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
  std::array<float, MaxSize> keys;
  std::array<item_t, MaxSize> sorted;
  std::array<std::uint32_t, MaxSize> places;
  // NOLINTEND(cppcoreguidelines-pro-type-member-init)
  for (std::size_t i = 0; i < n; ++i) {
    keys[i] = keyOf(items[i]);
  }
  for (std::size_t i = n; i < MaxSize; ++i) {
    keys[i] = std::numeric_limits<float>::infinity();
  }
  std::uint64_t taken = 0;
  for (std::size_t i = 0; i < n; ++i) {
    std::uint32_t place = 0;
    for (std::size_t j = 0; j < MaxSize; ++j) {
      place += keys[j] < keys[i] ? 1 : 0;
    }
    places[i] = place;
    taken |= std::uint64_t{1} << place;
  }
  // Equal keys share a place: only then count the equal keys before each
  // item, which keeps their order.
  if (static_cast<std::size_t>(std::popcount(taken)) != n) {
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < i; ++j) {
        places[i] += keys[j] == keys[i] ? 1 : 0;
      }
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    sorted[places[i]] = std::move(items[i]);
  }
  std::move(sorted.begin(), sorted.begin() + n, items.begin());
}

/// Sorts @p items stably by the float key @p keyOf gives each, given the
/// smallest and largest key: the items are spread over 2 * bit_ceil(n)
/// buckets across the key range, then an insertion pass orders them within
/// each bucket. Lists longer than 1024 items, key ranges that cannot be
/// scaled, and lists where one bucket would hold more than 32 items (keys
/// far from even, e.g. one outlier) use std::ranges::stable_sort instead, so
/// that the insertion pass stays short.
///
/// Every key must be finite.
///
/// @param items The items to sort
/// @param scratch Storage of the same type, of any content, reused between calls
/// @param keyOf Returns the float key of an item
/// @param minKey The smallest key
/// @param maxKey The largest key
template <typename item_t, typename key_of_t>
void bucketSortByFloat(std::vector<item_t>& items, std::vector<item_t>& scratch,
                       key_of_t keyOf, float minKey, float maxKey) {
  if (!(maxKey > minKey)) {
    return;  // fewer than two items, or all keys equal
  }

  // Tuned on AMD EPYC 7702 and 7713 (Zen 2 and 3) on ITk ttbar events: more
  // items, or fewer buckets per item, sorted no faster.
  constexpr std::size_t kMaxItems = 1024;
  constexpr std::size_t kBucketsPerItem = 2;
  constexpr std::uint32_t kMaxItemsPerBucket = 32;
  // Keep maxKey - minKey and the scale finite: an overflow traps where
  // floating point exceptions are enabled.
  const bool scalable =
      maxKey < 0x1p125f && minKey > -0x1p125f && maxKey - minKey >= 0x1p-100f;
  if (items.size() > kMaxItems || !scalable) {
    std::ranges::stable_sort(items, {}, keyOf);
    return;
  }

  const std::size_t nBuckets = kBucketsPerItem * std::bit_ceil(items.size());
  const float scale = static_cast<float>(nBuckets) / (maxKey - minKey);
  // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
  std::array<std::uint32_t, kBucketsPerItem * kMaxItems> bucketStart;
  std::array<std::uint16_t, kMaxItems> bucketOf;
  // NOLINTEND(cppcoreguidelines-pro-type-member-init)
  std::fill_n(bucketStart.begin(), nBuckets, 0u);
  for (std::size_t i = 0; i < items.size(); ++i) {
    const auto bucket =
        static_cast<std::size_t>((keyOf(items[i]) - minKey) * scale);
    bucketOf[i] = static_cast<std::uint16_t>(std::min(bucket, nBuckets - 1));
    ++bucketStart[bucketOf[i]];
  }
  // Turn the counts into each bucket's first position.
  std::uint32_t start = 0;
  std::uint32_t largest = 0;
  for (std::size_t b = 0; b < nBuckets; ++b) {
    largest = std::max(largest, bucketStart[b]);
    start += std::exchange(bucketStart[b], start);
  }
  if (largest > kMaxItemsPerBucket) {
    std::ranges::stable_sort(items, {}, keyOf);
    return;
  }
  scratch.resize(items.size());
  for (std::size_t i = 0; i < items.size(); ++i) {
    scratch[bucketStart[bucketOf[i]]++] = std::move(items[i]);
  }
  items.swap(scratch);

  // Insertion pass: moves an item only past larger keys, so it stays stable.
  insertionSortByFloat(items, keyOf);
}

/// Fills @p items with the @p n items @p itemAt gives, in order, and sorts
/// them stably by the float key @p keyOf gives each. The key range the
/// bucket sort needs is found while the list is filled, so every key is read
/// once. By size, measured on AMD EPYC 7702 with GCC 13 and 15 at -O2 and
/// -O3: up to 4 items an insertion sort is fastest, up to 40 the rank sort,
/// above that the bucket sort. A list with a key that is not finite is
/// sorted by std::ranges::stable_sort.
///
/// @param items Filled with the sorted items
/// @param scratch Storage of the same type, of any content, reused between calls
/// @param n The number of items
/// @param itemAt Returns item k, for k from 0 to n - 1
/// @param keyOf Returns the float key of an item
template <typename item_t, typename item_at_t, typename key_of_t>
void fillAndSortByFloat(std::vector<item_t>& items,
                        std::vector<item_t>& scratch, std::size_t n,
                        item_at_t itemAt, key_of_t keyOf) {
  constexpr std::size_t kInsertionSortMaxSize = 4;
  constexpr std::size_t kRankSortMaxSize = 40;
  items.clear();
  items.reserve(n);
  float minKey = std::numeric_limits<float>::max();
  float maxKey = std::numeric_limits<float>::lowest();
  bool finite = true;
  for (std::size_t k = 0; k < n; ++k) {
    items.push_back(itemAt(k));
    const float key = keyOf(items.back());
    // Checked first: an ordered comparison with a NaN traps where floating
    // point exceptions are enabled.
    if (std::isfinite(key)) {
      minKey = std::min(minKey, key);
      maxKey = std::max(maxKey, key);
    } else {
      finite = false;
    }
  }
  if (!finite) {
    std::ranges::stable_sort(items, {}, keyOf);
  } else if (n <= kInsertionSortMaxSize) {
    insertionSortByFloat(items, keyOf);
  } else if (n <= kRankSortMaxSize) {
    rankSortByFloat<kRankSortMaxSize>(items, keyOf);
  } else {
    bucketSortByFloat(items, scratch, keyOf, minKey, maxKey);
  }
}

}  // namespace Acts::detail
