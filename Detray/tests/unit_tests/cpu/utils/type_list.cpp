// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s).
#include "detray/utils/type_list.hpp"

// Google Test include(s).
#include <gtest/gtest.h>

// System include(s)
#include <iostream>

// Test type list implementation
GTEST_TEST(detray_utils, type_list) {
  using namespace detray;

  using list = types::list<float, int, double>;
  static_assert(types::contains<list, float>,
                "List should contain type 'float'");
  static_assert(types::position<list, float> == 0u,
                "Wrong position for type 'float'");

  static_assert(types::contains<list, int>, "List should contain type 'int'");
  static_assert(types::position<list, int> == 1u,
                "Wrong position for type 'int'");

  static_assert(types::contains<list, double>,
                "List should contain type 'double'");
  static_assert(types::position<list, double> == 2u,
                "Wrong position for type 'double'");

  static_assert(!types::contains<list, char>,
                "List should not contain type 'char'");
  static_assert(
      types::position<list, char> == std::numeric_limits<std::size_t>::max(),
      "Wrong position for type 'char'");

  static_assert(!types::contains<list, bool>,
                "List should not contain type 'bool'");
  static_assert(
      types::position<list, bool> == std::numeric_limits<std::size_t>::max(),
      "Wrong position for type 'bool'");

  static_assert(std::is_same_v<types::front<list>, float>,
                "Could not access type list front");

  static_assert(std::is_same_v<types::back<list>, double>,
                "Could not access type list back");

  static_assert(std::is_same_v<types::push_back<list, char>,
                               types::list<float, int, double, char>>,
                "Failed to push back new type");

  static_assert(std::is_same_v<types::push_back_unique<list, char>,
                               types::list<float, int, double, char>>,
                "Failed to push back unique new type");

  static_assert(std::is_same_v<types::push_back_unique<list, int>, list>,
                "Failed to preserve list for duplicate pushed back type");

  static_assert(std::is_same_v<types::push_front<list, char>,
                               types::list<char, float, int, double>>,
                "Failed to push front new type");

  static_assert(types::size<list> == 3ul, "Incorrect size");

  static_assert(std::is_same_v<types::at<list, 1>, int>, "Failed access type");

  types::print<list>();
  types::print<list>(false);

  // Print with template params
  std::clog << types::get_name<list>(true) << std::endl;
  // Print without template params
  std::clog << types::get_name<list>() << std::endl;
}
