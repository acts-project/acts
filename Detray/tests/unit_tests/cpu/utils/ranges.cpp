// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// detray core
#include "detray/utils/ranges.hpp"

#include "detray/definitions/containers.hpp"

// Vecmem include(s)
#include "vecmem/containers/device_vector.hpp"
#include "vecmem/containers/jagged_device_vector.hpp"
#include "vecmem/containers/jagged_vector.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <forward_list>
#include <list>
#include <type_traits>

using namespace detray;

// Test basic ranges definitions
GTEST_TEST(detray_utils, ranges) {
  //
  // detray containers
  //
  // std::vector
  static_assert(detray::ranges::range<dvector<float>>);
  static_assert(!detray::ranges::view<dvector<float>>);
  static_assert(detray::ranges::random_access_range<dvector<float>>);
  static_assert(detray::ranges::common_range<dvector<float>>);

  // std::array
  static_assert(detray::ranges::range<darray<float, 1>>);
  static_assert(!detray::ranges::view<darray<float, 1>>);
  static_assert(detray::ranges::random_access_range<darray<float, 1>>);
  static_assert(detray::ranges::common_range<darray<float, 1>>);

  // std::map
  static_assert(detray::ranges::range<dmap<int, float>>);
  static_assert(!detray::ranges::view<dmap<int, float>>);
  static_assert(detray::ranges::bidirectional_range<dmap<int, float>>);
  static_assert(!detray::ranges::random_access_range<dmap<int, float>>);
  static_assert(detray::ranges::common_range<dmap<int, float>>);

  // std::tuple
  static_assert(!detray::ranges::range<dtuple<int, float>>);
  static_assert(!detray::ranges::view<dtuple<int, float>>);

  //
  // vecmem containers
  //

  // vecmem::device_vector
  static_assert(detray::ranges::range<vecmem::device_vector<float>>);
  static_assert(!detray::ranges::view<vecmem::device_vector<float>>);
  static_assert(
      detray::ranges::random_access_range<vecmem::device_vector<float>>);
  static_assert(detray::ranges::common_range<vecmem::device_vector<float>>);

  // vecmem::jagged_vector
  static_assert(detray::ranges::range<djagged_vector<float>>);
  static_assert(!detray::ranges::view<djagged_vector<float>>);
  static_assert(detray::ranges::bidirectional_range<djagged_vector<float>>);
  static_assert(detray::ranges::common_range<djagged_vector<float>>);

  // vecmem::jagged_device_vector
  static_assert(detray::ranges::range<vecmem::jagged_device_vector<float>>);
  static_assert(!detray::ranges::view<vecmem::jagged_device_vector<float>>);
  static_assert(
      detray::ranges::bidirectional_range<vecmem::jagged_device_vector<float>>);
  static_assert(
      detray::ranges::common_range<vecmem::jagged_device_vector<float>>);

  //
  // Some additional STL containers
  //

  // std::forward_list
  static_assert(detray::ranges::range<std::forward_list<float>>);
  static_assert(!detray::ranges::view<std::forward_list<float>>);
  static_assert(detray::ranges::forward_range<std::forward_list<float>>);
  static_assert(!detray::ranges::bidirectional_range<std::forward_list<float>>);
  static_assert(detray::ranges::common_range<std::forward_list<float>>);

  // std::list
  static_assert(detray::ranges::range<std::list<float>>);
  static_assert(!detray::ranges::view<std::list<float>>);
  static_assert(detray::ranges::bidirectional_range<std::list<float>>);
  static_assert(!detray::ranges::random_access_range<std::list<float>>);
  static_assert(detray::ranges::common_range<std::list<float>>);
}

// Unittest for an empty view
GTEST_TEST(detray_utils, ranges_empty) {
  auto ev = detray::views::empty<float>();

  // general tests
  static_assert(std::is_copy_assignable_v<decltype(ev)>);
  static_assert(detray::ranges::range<decltype(ev)>);
  static_assert(detray::ranges::view<decltype(ev)>);
  static_assert(detray::ranges::random_access_range<decltype(ev)>);

  // Test prerequisites for LagacyIterator
  static_assert(
      std::is_copy_constructible_v<typename decltype(ev)::iterator_t>);
  static_assert(std::is_copy_assignable_v<typename decltype(ev)::iterator_t>);
  static_assert(std::is_destructible_v<typename decltype(ev)::iterator_t>);

  // Test inherited member functions
  ASSERT_EQ(ev.size(), 0u);

  for (const auto i : ev) {
    ASSERT_TRUE(i != i);
  }
}

// Unittest for the generation of a single element sequence
GTEST_TEST(detray_utils, ranges_single) {
  const dindex value{251u};

  auto sngl = detray::views::single(value);

  // general tests
  static_assert(std::is_copy_assignable_v<decltype(sngl)>);
  static_assert(detray::ranges::range<decltype(sngl)>);
  static_assert(detray::ranges::view<decltype(sngl)>);
  static_assert(detray::ranges::random_access_range<decltype(sngl)>);

  // Test prerequisites for LagacyIterator
  static_assert(
      std::is_copy_constructible_v<typename decltype(sngl)::iterator_t>);
  static_assert(std::is_copy_assignable_v<typename decltype(sngl)::iterator_t>);
  static_assert(std::is_destructible_v<typename decltype(sngl)::iterator_t>);

  // Test inherited member functions
  ASSERT_EQ(sngl[0], value);
  ASSERT_EQ(sngl.size(), 1u);
  ASSERT_EQ(sngl.front(), 251u);
  ASSERT_EQ(sngl.back(), 251u);

  for (auto i : sngl) {
    ASSERT_EQ(251u, i);
  }
}

// Unittest for the generation of a pointer to an element
GTEST_TEST(detray_utils, ranges_pointer) {
  const dindex value{251u};

  auto ptr = detray::views::pointer(value);

  // general tests
  static_assert(std::is_copy_assignable_v<decltype(ptr)>);
  static_assert(detray::ranges::range<decltype(ptr)>);
  static_assert(detray::ranges::view<decltype(ptr)>);
  static_assert(detray::ranges::random_access_range<decltype(ptr)>);

  // Test prerequisites for LagacyIterator
  static_assert(
      std::is_copy_constructible_v<typename decltype(ptr)::iterator_t>);
  static_assert(std::is_copy_assignable_v<typename decltype(ptr)::iterator_t>);
  static_assert(std::is_destructible_v<typename decltype(ptr)::iterator_t>);

  // Test inherited member functions
  ASSERT_EQ(ptr[0], value);
  ASSERT_EQ(ptr.size(), 1u);
  ASSERT_EQ(ptr.front(), 251u);
  ASSERT_EQ(ptr.back(), 251u);

  for (auto i : ptr) {
    ASSERT_EQ(251u, i);
  }
}

// Unittest for the generation of a single element sequence
GTEST_TEST(detray_utils, ranges_iota_single) {
  dindex check{0u};
  dindex single{7u};

  auto seq = detray::views::iota(single);

  // general tests
  static_assert(std::is_copy_assignable_v<decltype(seq)>);
  static_assert(detray::ranges::range<decltype(seq)>);
  static_assert(detray::ranges::view<decltype(seq)>);
  static_assert(detray::ranges::forward_range<decltype(seq)>);
  static_assert(detray::ranges::bidirectional_range<decltype(seq)>);
  static_assert(detray::ranges::random_access_range<decltype(seq)>);

  // Test prerequisites for LagacyIterator
  static_assert(
      std::is_copy_constructible_v<typename decltype(seq)::iterator_t>);
  static_assert(std::is_copy_assignable_v<typename decltype(seq)::iterator_t>);
  static_assert(std::is_destructible_v<typename decltype(seq)::iterator_t>);

  // Test inherited member functions
  ASSERT_EQ(seq.size(), 1u);

  for (auto i : seq) {
    check += i;
  }
  ASSERT_EQ(check, single);
}

// Unittest for the generation of a sequence in an interval
GTEST_TEST(detray_utils, ranges_iota_interval) {
  darray<dindex, 2> interval = {2u, 7u};

  auto seq = detray::views::iota(interval);

  // general tests
  static_assert(detray::ranges::range<decltype(seq)>);
  static_assert(detray::ranges::view<decltype(seq)>);
  static_assert(std::is_copy_assignable_v<decltype(seq)>);
  static_assert(detray::ranges::forward_range<decltype(seq)>);
  static_assert(detray::ranges::bidirectional_range<decltype(seq)>);
  static_assert(detray::ranges::random_access_range<decltype(seq)>);

  // Test prerequisites for LagacyIterator
  static_assert(
      std::is_copy_constructible_v<typename decltype(seq)::iterator_t>);
  static_assert(std::is_copy_assignable_v<typename decltype(seq)::iterator_t>);
  static_assert(std::is_destructible_v<typename decltype(seq)::iterator_t>);

  // Test inherited member functions
  ASSERT_EQ(seq.size(), 5u);

  std::vector<dindex> reference = {2u, 3u, 4u, 5u, 6u};
  std::vector<dindex> check = {};
  for (auto i : seq) {
    check.push_back(i);
  }
  ASSERT_EQ(check.size(), reference.size());
  ASSERT_EQ(check, reference);
}

// Unittest for the generation of a cartesian product (trivial case)
GTEST_TEST(detray_utils, ranges_cartesian_product_trivial) {
  auto seq1 = detray::views::iota(dindex_range{1u, 2u});
  auto seq2 = detray::views::iota(dindex_range{2u, 3u});
  auto seq3 = detray::views::iota(dindex_range{3u, 4u});

  detray::views::cartesian_product cp{std::move(seq1), std::move(seq2),
                                      std::move(seq3)};

  // General tests
  static_assert(detray::ranges::range<decltype(cp)>);
  static_assert(detray::ranges::view<decltype(cp)>);
  static_assert(std::is_copy_assignable_v<decltype(cp)>);
  static_assert(detray::ranges::input_range<decltype(cp)>);
  static_assert(detray::ranges::forward_range<decltype(cp)>);
  static_assert(detray::ranges::bidirectional_range<decltype(cp)>);
  static_assert(!detray::ranges::random_access_range<decltype(cp)>);

  // Test prerequisites for LagacyIterator
  static_assert(
      std::is_copy_constructible_v<typename decltype(cp)::iterator_t>);
  static_assert(std::is_copy_assignable_v<typename decltype(cp)::iterator_t>);
  static_assert(std::is_destructible_v<typename decltype(cp)::iterator_t>);

  // Test size
  ASSERT_EQ(cp.size(), 1u);

  std::size_t r{0u};
  for (const auto [i, j, k] : cp) {
    ++r;
    ASSERT_EQ(i, 1u);
    ASSERT_EQ(j, 2u);
    ASSERT_EQ(k, 3u);
  }

  ASSERT_EQ(r, 1u);
}
// Unittest for the generation of a cartesian product from a range of intervals
GTEST_TEST(detray_utils, ranges_cartesian_product) {
  const dindex_range range1{2u, 7u};
  const dindex_range range2{1u, 10u};
  const dindex_range range3{3u, 4u};

  auto seq1 = detray::views::iota(range1);
  auto seq2 = detray::views::iota(range2);
  auto seq3 = detray::views::iota(range3);

  const std::size_t size{seq1.size() * seq2.size() * seq3.size()};

  detray::views::cartesian_product cp{std::move(seq1), std::move(seq2),
                                      std::move(seq3)};

  // General tests
  static_assert(detray::ranges::range<decltype(cp)>);
  static_assert(detray::ranges::view<decltype(cp)>);
  static_assert(std::is_copy_assignable_v<decltype(cp)>);
  static_assert(detray::ranges::input_range<decltype(cp)>);
  static_assert(detray::ranges::forward_range<decltype(cp)>);
  static_assert(detray::ranges::bidirectional_range<decltype(cp)>);
  static_assert(!detray::ranges::random_access_range<decltype(cp)>);

  // Test prerequisites for LagacyIterator
  static_assert(
      std::is_copy_constructible_v<typename decltype(cp)::iterator_t>);
  static_assert(std::is_copy_assignable_v<typename decltype(cp)::iterator_t>);
  static_assert(std::is_destructible_v<typename decltype(cp)::iterator_t>);

  // Test size
  ASSERT_EQ(cp.size(), size);

  // Generate truth
  std::vector<std::tuple<dindex, dindex, dindex>> result;
  for (auto i : detray::views::iota(range1)) {
    for (auto j : detray::views::iota(range2)) {
      for (auto k : detray::views::iota(range3)) {
        result.emplace_back(i, j, k);
      }
    }
  }

  std::size_t r{0u};
  for (const auto [i, j, k] : cp) {
    const auto [l, m, n] = result[r];
    ++r;
    ASSERT_EQ(i, l);
    ASSERT_EQ(j, m);
    ASSERT_EQ(k, n);
  }
}

// Unittest for the convenience enumeration of a range
GTEST_TEST(detray_utils, ranges_enumerate) {
  struct uint_holder {
    unsigned int ui{0u};
  };

  dvector<uint_holder> seq = {{0u}, {1u}, {2u}, {3u}, {4u}, {5u}};

  auto enumerator = detray::views::enumerate(seq);

  // general tests
  static_assert(detray::ranges::range<decltype(enumerator)>);
  static_assert(detray::ranges::view<decltype(enumerator)>);
  static_assert(std::is_copy_assignable_v<decltype(enumerator)>);
  static_assert(detray::ranges::random_access_range<decltype(enumerator)>);

  // Test prerequisites for LagacyIterator
  static_assert(
      std::is_copy_constructible_v<typename decltype(enumerator)::iterator_t>);
  static_assert(
      std::is_copy_assignable_v<typename decltype(enumerator)::iterator_t>);
  static_assert(
      std::is_destructible_v<typename decltype(enumerator)::iterator_t>);

  // Test inherited member functions
  const auto [i, v] = enumerator[2];
  ASSERT_EQ(i, 2u);
  ASSERT_EQ(v.ui, 2u);
  ASSERT_EQ(enumerator.size(), 6u);
  const auto [i_front, v_front] = enumerator.front();
  ASSERT_EQ(i_front, 0u);
  ASSERT_EQ(v_front.ui, 0u);
  const auto [i_back, v_back] = enumerator.back();
  ASSERT_EQ(i_back, 5u);
  ASSERT_EQ(v_back.ui, 5u);

  for (auto [j, w] : enumerator) {
    ASSERT_TRUE(j == w.ui);
    ASSERT_EQ(j, w.ui);
  }
}

// Integration test for the picking of indexed elements from another range
GTEST_TEST(detray_utils, ranges_pick) {
  // The indices of the iota elements to be picked
  std::vector<dindex> indices = {2u, 3u, 7u, 8u};
  std::vector<dindex> check = {};

  struct uint_holder {
    unsigned int ui{0u};
  };

  dvector<uint_holder> seq = {{0u}, {1u}, {2u}, {3u}, {4u},
                              {5u}, {6u}, {7u}, {8u}};

  auto selected = detray::views::pick(seq, indices);

  // general tests
  static_assert(detray::ranges::range<decltype(selected)>);
  static_assert(detray::ranges::view<decltype(selected)>);
  static_assert(std::is_copy_assignable_v<decltype(selected)>);
  static_assert(detray::ranges::random_access_range<decltype(selected)>);

  // Test prerequisites for LagacyIterator
  static_assert(
      std::is_copy_constructible_v<typename decltype(selected)::iterator_t>);
  static_assert(
      std::is_copy_assignable_v<typename decltype(selected)::iterator_t>);
  static_assert(
      std::is_destructible_v<typename decltype(selected)::iterator_t>);

  // Test inherited member functions
  const auto [i, v] = selected[2];
  ASSERT_EQ(i, 7u);
  ASSERT_EQ(v.ui, 7u);
  ASSERT_EQ(selected.size(), 4u);
  const auto [i_front, v_front] = selected.front();
  ASSERT_EQ(i_front, 2u);
  ASSERT_EQ(v_front.ui, 2u);
  const auto [i_back, v_back] = selected.back();
  ASSERT_EQ(i_back, 8u);
  ASSERT_EQ(v_back.ui, 8u);

  for (auto [j, w] : selected) {
    ASSERT_EQ(j, w.ui);
    check.push_back(w.ui);
  }
  ASSERT_EQ(check.size(), indices.size());
  ASSERT_EQ(check, indices);
}

// Unittest for the joining of multiple ranges
GTEST_TEST(detray_utils, ranges_join) {
  dvector<dindex> interval_0 = {};
  dvector<dindex> interval_1 = {2u, 3u, 4u};
  dvector<dindex> interval_2 = {7u, 8u, 9u};
  dvector<dindex> interval_3 = {10u, 11u, 12u, 13u};
  dvector<dvector<dindex>> intervals{interval_0, interval_0, interval_1,
                                     interval_2, interval_0, interval_3,
                                     interval_0};

  std::vector<dindex> reference = {2u, 3u, 4u, 7u, 8u, 9u, 10u, 11u, 12u, 13u};
  std::vector<dindex> check = {};

  auto joined = detray::views::join(intervals);

  // general tests
  static_assert(detray::ranges::range<decltype(joined)>);
  static_assert(detray::ranges::view<decltype(joined)>);
  static_assert(std::is_copy_assignable_v<decltype(joined)>);
  static_assert(detray::ranges::random_access_range<decltype(joined)>);

  // Test prerequisites for LagacyIterator
  static_assert(
      std::is_copy_constructible_v<typename decltype(joined)::iterator_t>);
  static_assert(
      std::is_copy_assignable_v<typename decltype(joined)::iterator_t>);
  static_assert(std::is_destructible_v<typename decltype(joined)::iterator_t>);

  // Test inherited member functions
  ASSERT_EQ(joined[1], 3u);
  ASSERT_EQ(joined[4], 8u);
  ASSERT_EQ(joined[7], 11u);
  ASSERT_EQ(joined.size(), 10u);
  ASSERT_EQ(joined.front(), 2u);
  ASSERT_EQ(joined.back(), 13u);

  for (auto j : joined) {
    static_assert(!std::is_const_v<decltype(j)>, "Non-const element in join");
    check.push_back(j);
  }
  ASSERT_EQ(check.size(), reference.size());
  ASSERT_EQ(check, reference);

  auto joined_comp = intervals | detray::views::join();

  check.clear();
  for (auto j : joined_comp) {
    static_assert(!std::is_const_v<decltype(j)>, "Non-const element in join");
    check.push_back(j);
  }
  ASSERT_EQ(check.size(), reference.size());
  ASSERT_EQ(check, reference);
}

// Unittest for the joining of multiple ranges
GTEST_TEST(detray_utils, ranges_static_join) {
  dvector<dindex> interval_1 = {2u, 3u, 4u};
  dvector<dindex> interval_2 = {7u, 8u, 9u};

  std::vector<dindex> reference = {2u, 3u, 4u, 7u, 8u, 9u};
  std::vector<dindex> check = {};

  auto joined = detray::views::static_join(interval_1, interval_2);

  // general tests
  static_assert(detray::ranges::range<decltype(joined)>);
  static_assert(detray::ranges::view<decltype(joined)>);
  static_assert(std::is_copy_assignable_v<decltype(joined)>);
  static_assert(detray::ranges::random_access_range<decltype(joined)>);

  // Test prerequisites for LagacyIterator
  static_assert(
      std::is_copy_constructible_v<typename decltype(joined)::iterator_t>);
  static_assert(
      std::is_copy_assignable_v<typename decltype(joined)::iterator_t>);
  static_assert(std::is_destructible_v<typename decltype(joined)::iterator_t>);

  // Test inherited member functions
  ASSERT_EQ(joined[1], 3u);
  ASSERT_EQ(joined[4], 8u);
  ASSERT_EQ(joined.size(), 6u);
  ASSERT_EQ(joined.front(), 2u);
  ASSERT_EQ(joined.back(), 9u);

  for (const auto j : joined) {
    check.push_back(j);
  }
  ASSERT_EQ(check.size(), reference.size());
  ASSERT_EQ(check, reference);
}

// Unittest for the subrange implementation
GTEST_TEST(detray_utils, ranges_subrange) {
  std::size_t begin{1u};
  std::size_t end{4u};
  std::array<std::size_t, 2> interval{begin, end};

  dvector<int> seq = {0, 1, 2, 3, 4, 5};

  auto sr = detray::ranges::subrange(seq, interval);

  // general tests
  static_assert(detray::ranges::range<decltype(sr)>);
  static_assert(detray::ranges::view<decltype(sr)>);
  static_assert(std::is_copy_assignable_v<decltype(sr)>);
  static_assert(detray::ranges::random_access_range<decltype(sr)>);

  // Test prerequisites for LagacyIterator
  static_assert(
      std::is_copy_constructible_v<typename decltype(sr)::iterator_t>);
  static_assert(std::is_copy_assignable_v<typename decltype(sr)::iterator_t>);
  static_assert(std::is_destructible_v<typename decltype(sr)::iterator_t>);

  ASSERT_EQ(sr[1], seq[begin + 1u]);
  ASSERT_EQ(sr.size(), 3u);
  ASSERT_EQ(sr.front(), 1u);
  ASSERT_EQ(sr.back(), 3u);

  // non-const iteration
  std::size_t i{1u};
  for (const auto& v : sr) {
    ASSERT_NE(v, 0u);
    ASSERT_NE(v, 4u);
    ASSERT_EQ(v, seq[i++]);
  }

  // const iteration
  const dvector<int> seq_c(seq);
  i = 1u;
  for (const auto& v : detray::ranges::subrange(seq_c, interval)) {
    ASSERT_EQ(v, seq[i++]);
  }

  // Check that an empty subrange does not get iterated
  interval[0] = 0u;
  interval[1] = 0u;
  for (const auto& v : detray::ranges::subrange(seq_c, interval)) {
    (void)v;
    ASSERT_TRUE(false);
  }
}
