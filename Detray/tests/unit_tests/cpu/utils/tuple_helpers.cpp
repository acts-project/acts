// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/utils/tuple_helpers.hpp"

// System include(s)
#include <cassert>
#include <string>
#include <tuple>
#include <type_traits>

// GoogleTest include(s).
#include <gtest/gtest.h>

GTEST_TEST(detray_utils, tuple_helpers) {
  using namespace detray;

  // std::tuple test
  auto s_tuple =
      detail::make_tuple<std::tuple>(2.0f, -3L, std::string("std::tuple"), 4UL);
  static_assert(
      std::is_same_v<std::tuple<float, long, std::string, unsigned long>,
                     decltype(s_tuple)>,
      "detail::make_tuple failed for std::tuple");

  static_assert(std::is_same_v<detail::tuple_element_t<2, decltype(s_tuple)>,
                               std::string>,
                "detail::tuple_element retrieval failed for std::tuple");

  const auto s_tuple_size = detail::tuple_size_v<decltype(s_tuple)>;
  EXPECT_EQ(s_tuple_size, 4UL);

  EXPECT_FLOAT_EQ(detail::get<0>(s_tuple), 2.0f);
  EXPECT_EQ(detail::get<1>(s_tuple), -3L);
  EXPECT_EQ(detail::get<2>(s_tuple), std::string("std::tuple"));
  EXPECT_EQ(detail::get<3>(s_tuple), 4UL);
  EXPECT_FLOAT_EQ(detail::get<float>(s_tuple), 2.0f);
  EXPECT_EQ(detail::get<long>(s_tuple), -3L);
  EXPECT_EQ(detail::get<std::string>(s_tuple), std::string("std::tuple"));
  EXPECT_EQ(detail::get<unsigned long>(s_tuple), 4UL);

  // dtuple test
  auto d_tuple = detail::make_tuple<dtuple>(1.0f, 2UL, std::string("dtuple"));
  static_assert(std::is_same_v<dtuple<float, unsigned long, std::string>,
                               decltype(d_tuple)>,
                "detail::make_tuple failed for dtuple");

  static_assert(std::is_same_v<detail::tuple_element_t<1, decltype(d_tuple)>,
                               unsigned long>,
                "detail::tuple_element retrieval failed for dtuple");

  const auto d_tuple_size = detail::tuple_size_v<decltype(d_tuple)>;
  EXPECT_EQ(d_tuple_size, 3UL);

  EXPECT_FLOAT_EQ(detail::get<0>(d_tuple), 1.0f);
  EXPECT_EQ(detail::get<1>(d_tuple), 2UL);
  EXPECT_EQ(detail::get<2>(d_tuple), std::string("dtuple"));
  EXPECT_FLOAT_EQ(detail::get<float>(d_tuple), 1.0f);
  EXPECT_EQ(detail::get<unsigned long>(d_tuple), 2UL);
  EXPECT_EQ(detail::get<std::string>(d_tuple), std::string("dtuple"));

  // Check type concatenation
  static_assert(
      std::same_as<
          detail::tuple_cat_t<std::tuple<int, double>, std::tuple<float>,
                              std::tuple<>, std::tuple<char, bool, double>>,
          std::tuple<int, double, float, char, bool, double>>);

  static_assert(
      std::same_as<detail::tuple_cat_t<dtuple<int, double>, dtuple<float>,
                                       dtuple<>, dtuple<char, bool, double>>,
                   dtuple<int, double, float, char, bool, double>>);

  // Permutation check
  static_assert(detail::is_permutation_v<dtuple<>, dtuple<>>);
  static_assert(detail::is_permutation_v<dtuple<int>, dtuple<int>>);
  static_assert(!detail::is_permutation_v<dtuple<int>, dtuple<int, int>>);
  static_assert(detail::is_permutation_v<dtuple<int, float, double>,
                                         dtuple<double, int, float>>);
  static_assert(!detail::is_permutation_v<dtuple<int, float, double>,
                                          dtuple<char, int, float>>);
  static_assert(!detail::is_permutation_v<dtuple<int, float, double>,
                                          dtuple<int, int, double, float>>);

  // Check unique element tuple type
  static_assert(std::same_as<detail::unique_t<std::tuple<>>, std::tuple<>>);
  static_assert(
      std::same_as<detail::unique_t<std::tuple<int>>, std::tuple<int>>);
  static_assert(std::same_as<detail::unique_t<std::tuple<int, float, double>>,
                             std::tuple<int, float, double>>);
  static_assert(
      std::same_as<detail::unique_t<std::tuple<int, float, int, int, double>>,
                   std::tuple<int, float, double>>);

  static_assert(std::same_as<detail::unique_t<dtuple<>>, dtuple<>>);
  static_assert(std::same_as<detail::unique_t<dtuple<int>>, dtuple<int>>);
  static_assert(std::same_as<detail::unique_t<dtuple<int, float, double>>,
                             dtuple<int, float, double>>);
  static_assert(std::same_as<detail::unique_t<dtuple<int, int>>, dtuple<int>>);
  static_assert(
      std::same_as<detail::unique_t<dtuple<int, int, int>>, dtuple<int>>);
  static_assert(std::same_as<detail::unique_t<dtuple<int, float, int, int>>,
                             dtuple<int, float>>);
  static_assert(
      std::same_as<detail::unique_t<dtuple<int, float, int, int, double>>,
                   dtuple<int, float, double>>);
}
