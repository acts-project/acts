// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/FunctionComposition.hpp"

#include <memory>
#include <string>

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_AUTO_TEST_CASE(ComposeTwoFunctions) {
  auto addOne = [](int x) { return x + 1; };
  auto doubleValue = [](int x) { return 2 * x; };

  auto composed = Acts::compose(doubleValue, addOne);

  BOOST_CHECK_EQUAL(composed(3), 8);
}

BOOST_AUTO_TEST_CASE(ComposeMultipleFunctions) {
  auto decorate = [](const std::string& s) { return "[" + s + "]"; };
  auto appendB = [](const std::string& s) { return s + "b"; };
  auto appendA = [](const std::string& s) { return s + "a"; };

  auto composed = Acts::compose(decorate, appendB, appendA);

  BOOST_CHECK_EQUAL(composed(std::string{"x"}), "[xab]");
}

BOOST_AUTO_TEST_CASE(ComposeForwardsMoveOnlyInput) {
  auto timesTwo = [](int x) { return x * 2; };
  auto takeOwnership = [](std::unique_ptr<int> ptr) { return *ptr + 1; };

  auto composed = Acts::compose(timesTwo, takeOwnership);

  BOOST_CHECK_EQUAL(composed(std::make_unique<int>(4)), 10);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
