// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s).
#include "detray/utils/invalid_values.hpp"

// Google Test include(s).
#include <gtest/gtest.h>

// System include(s)
#include <cstdint>
#include <limits>

using namespace detray;

// Test sort functions
GTEST_TEST(detray_utils, invalid_values) {
  ASSERT_EQ(std::numeric_limits<int>::max(), detail::invalid_value<int>());
  ASSERT_EQ(std::numeric_limits<std::int_least8_t>::max(),
            detail::invalid_value<std::int_least8_t>());
  ASSERT_EQ(std::numeric_limits<std::int_least16_t>::max(),
            detail::invalid_value<std::int_least16_t>());
  ASSERT_EQ(std::numeric_limits<std::int_least32_t>::max(),
            detail::invalid_value<std::int_least32_t>());
  ASSERT_EQ(std::numeric_limits<std::int_least64_t>::max(),
            detail::invalid_value<std::int_least64_t>());
  ASSERT_EQ(std::numeric_limits<std::intmax_t>::max(),
            detail::invalid_value<std::intmax_t>());

  ASSERT_EQ(std::numeric_limits<unsigned int>::max(),
            detail::invalid_value<unsigned int>());
  ASSERT_EQ(std::numeric_limits<std::size_t>::max(),
            detail::invalid_value<std::size_t>());
  ASSERT_EQ(std::numeric_limits<std::uint_least8_t>::max(),
            detail::invalid_value<std::uint_least8_t>());
  ASSERT_EQ(std::numeric_limits<std::uint_least16_t>::max(),
            detail::invalid_value<std::uint_least16_t>());
  ASSERT_EQ(std::numeric_limits<std::uint_least32_t>::max(),
            detail::invalid_value<std::uint_least32_t>());
  ASSERT_EQ(std::numeric_limits<std::uint_least64_t>::max(),
            detail::invalid_value<std::uint_least64_t>());
  ASSERT_EQ(std::numeric_limits<std::uintmax_t>::max(),
            detail::invalid_value<std::uintmax_t>());

  ASSERT_EQ(std::numeric_limits<float>::max(), detail::invalid_value<float>());
  ASSERT_EQ(std::numeric_limits<double>::max(),
            detail::invalid_value<double>());

  ASSERT_TRUE(detail::is_invalid_value(std::numeric_limits<int>::max()));
  ASSERT_TRUE(
      detail::is_invalid_value(std::numeric_limits<std::int_least8_t>::max()));
  ASSERT_TRUE(
      detail::is_invalid_value(std::numeric_limits<std::int_least16_t>::max()));
  ASSERT_TRUE(
      detail::is_invalid_value(std::numeric_limits<std::int_least32_t>::max()));
  ASSERT_TRUE(
      detail::is_invalid_value(std::numeric_limits<std::int_least64_t>::max()));
  ASSERT_TRUE(
      detail::is_invalid_value(std::numeric_limits<unsigned int>::max()));
  ASSERT_TRUE(
      detail::is_invalid_value(std::numeric_limits<std::size_t>::max()));
  ASSERT_TRUE(
      detail::is_invalid_value(std::numeric_limits<std::uint_least8_t>::max()));
  ASSERT_TRUE(detail::is_invalid_value(
      std::numeric_limits<std::uint_least16_t>::max()));
  ASSERT_TRUE(detail::is_invalid_value(
      std::numeric_limits<std::uint_least32_t>::max()));
  ASSERT_TRUE(detail::is_invalid_value(
      std::numeric_limits<std::uint_least64_t>::max()));
  ASSERT_TRUE(
      detail::is_invalid_value(std::numeric_limits<std::uintmax_t>::max()));
  ASSERT_TRUE(detail::is_invalid_value(std::numeric_limits<float>::max()));
  ASSERT_TRUE(detail::is_invalid_value(std::numeric_limits<double>::max()));

  ASSERT_FALSE(detail::is_invalid_value(1));
  ASSERT_FALSE(detail::is_invalid_value(0));
  ASSERT_FALSE(detail::is_invalid_value(-1));
  ASSERT_FALSE(detail::is_invalid_value(1u));
  ASSERT_FALSE(detail::is_invalid_value(1ul));
  ASSERT_FALSE(detail::is_invalid_value(1.));
  ASSERT_FALSE(detail::is_invalid_value(1.f));
}
