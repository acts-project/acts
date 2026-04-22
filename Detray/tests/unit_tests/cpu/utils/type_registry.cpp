// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/utils/type_registry.hpp"

// Google Test include(s).
#include <gtest/gtest.h>

// System include(s)
#include <cstdint>
#include <type_traits>

// Test IDs
enum class test_ids : std::uint_least8_t {
  e_int = 0u,
  e_double = 1u,
  e_float1 = 2u,
  e_float2 = 3u,
};

// Test type registry implementation
GTEST_TEST(detray_utils, type_registry) {
  using namespace detray;
  using enum test_ids;

  constexpr auto int_idx{static_cast<std::size_t>(e_int)};
  constexpr auto double_idx{static_cast<std::size_t>(e_double)};
  constexpr auto float1_idx{static_cast<std::size_t>(e_float1)};
  constexpr auto float2_idx{static_cast<std::size_t>(e_float2)};

  using type_registry_t = types::registry<test_ids, int, double, float, float>;

  static_assert(type_registry_t::n_types == 4u,
                "Incorrect number of types in registry");
  static_assert(types::size<type_registry_t> == 4u,
                "Incorrect number of types in registry");

  // contains
  static_assert(types::contains<type_registry_t, int>,
                "'contains' failed for 'int' type");
  static_assert(types::contains<type_registry_t, double>,
                "'contains' failed for 'double' type");
  static_assert(types::contains<type_registry_t, float>,
                "'contains' failed for 'float' type");

  static_assert(!types::contains<type_registry_t, char>,
                "'contains' failed for 'char' type");
  static_assert(!types::contains<type_registry_t, void>,
                "'contains' failed for 'void' type");

  // Is valid
  static_assert(type_registry_t::is_valid(int_idx), "Index for 'int' invalid");
  static_assert(type_registry_t::is_valid(double_idx),
                "Index for 'double' invalid");
  static_assert(type_registry_t::is_valid(float1_idx),
                "Index for 'float 1' invalid");
  static_assert(type_registry_t::is_valid(float2_idx),
                "Index for 'float 2' invalid");
  static_assert(!type_registry_t::is_valid(5u), "ID '5' not invalid");

  static_assert(type_registry_t::is_valid(e_int), "ID for 'int' invalid");
  static_assert(type_registry_t::is_valid(e_double), "ID for 'double' invalid");
  static_assert(type_registry_t::is_valid(e_float1),
                "ID for 'float 1' invalid");
  static_assert(type_registry_t::is_valid(e_float2),
                "ID for 'float 2' invalid");
  static_assert(!type_registry_t::is_valid(5u), "ID '5' not invalid");

  // Get index (position)
  static_assert(types::position<type_registry_t, int> == int_idx,
                "Position for type 'int' incorrect");
  static_assert(types::position<type_registry_t, double> == double_idx,
                "Position for type 'double' incorrect");
  static_assert(types::position<type_registry_t, float> == float1_idx,
                "Position for type 'float' incorrect");

  // Get ID
  static_assert(types::id<type_registry_t, int> == e_int,
                "ID for type 'int' incorrect");
  static_assert(types::id<type_registry_t, double> == e_double,
                "ID for type 'double' incorrect");
  static_assert(types::id<type_registry_t, float> == e_float1,
                "ID for type 'float' incorrect");

  // Get type at position
  static_assert(std::same_as<int, types::at<type_registry_t, int_idx>>,
                "Got incorrect type at position 0");
  static_assert(std::same_as<double, types::at<type_registry_t, double_idx>>,
                "Got incorrect type at position 1");
  static_assert(std::same_as<float, types::at<type_registry_t, float1_idx>>,
                "Got incorrect type at position 2");
  static_assert(std::same_as<float, types::at<type_registry_t, float2_idx>>,
                "Got incorrect type at position 3");

  // Get type from ID
  static_assert(std::same_as<int, types::get<type_registry_t, e_int>>,
                "Got incorrect type for 'e_int'");
  static_assert(std::same_as<double, types::get<type_registry_t, e_double>>,
                "Got incorrect type for 'e_double'");
  static_assert(std::same_as<float, types::get<type_registry_t, e_float1>>,
                "Got incorrect type for 'e_float1'");
  static_assert(std::same_as<float, types::get<type_registry_t, e_float2>>,
                "Got incorrect type for 'e_float2'");
}

namespace detray::test {

struct visitor {
  double operator()(const int &, int arg1, double arg2) const {
    return static_cast<double>(arg1) + arg2;
  }

  double operator()(const double &, int arg1, double arg2) const {
    return static_cast<double>(arg1) * arg2;
  }

  double operator()(const float &, int arg1, double arg2) const {
    return static_cast<double>(arg1) / arg2;
  }
};

}  // namespace detray::test

// Test the type registry visitor
GTEST_TEST(detray_utils, visit_type_registry) {
  using namespace detray;

  using type_registry_t = types::registry<test_ids, int, double, float, float>;

  double result = types::visit<type_registry_t, test::visitor>(0u, 2, 3.f);
  ASSERT_EQ(result, 5.);

  result = types::visit<type_registry_t, test::visitor>(1u, 2, 3.f);
  ASSERT_EQ(result, 6.);

  result = types::visit<type_registry_t, test::visitor>(2u, 2, 3.f);
  ASSERT_EQ(result, 2. / 3.);

  result = types::visit<type_registry_t, test::visitor>(3u, 2, 3.f);
  ASSERT_EQ(result, 2. / 3.);
}
