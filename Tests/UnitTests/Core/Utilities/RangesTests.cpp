// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Ranges.hpp"

#include <array>
#include <ranges>
#include <set>
#include <vector>

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(RangesTests)

BOOST_AUTO_TEST_CASE(ToVectorFromView) {
  std::array<int, 6> values = {1, 2, 3, 4, 5, 6};

  auto evenTimesTen =
      values | std::views::filter([](int value) { return value % 2 == 0; }) |
      std::views::transform([](int value) { return value * 10; });

  auto result = evenTimesTen | Acts::Ranges::to<std::vector>;

  std::vector<int> expected = {20, 40, 60};
  BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                result.begin(), result.end());
}

BOOST_AUTO_TEST_CASE(ToSetRemovesDuplicatesAndOrders) {
  std::vector<int> values = {4, 1, 4, 3, 2, 2};

  auto result = values | Acts::Ranges::to<std::set>;

  std::set<int> expected = {1, 2, 3, 4};
  BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                result.begin(), result.end());
}

BOOST_AUTO_TEST_CASE(ToAdaptorCallableSyntax) {
  std::vector<int> values = {5, 6, 7};

  auto adaptor = Acts::Ranges::to_adaptor<std::vector>{};
  auto result = adaptor(values);

  std::vector<int> expected = {5, 6, 7};
  BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                result.begin(), result.end());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
