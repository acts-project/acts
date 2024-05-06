// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Propagator/detail/Auctioneer.hpp"

#include <array>

namespace Acts::Test {

BOOST_AUTO_TEST_CASE(AuctioneerTest_VoidAuctioneer) {
  // Build arbitrary vector
  std::array<int, 4> vecArb = {0, 2, -5, 4};
  std::array<bool, 4> vecRes = {false, true, false, true};
  // Let it run through auction
  detail::VoidAuctioneer va;
  std::array<bool, 4> resultVa = va(vecArb);
  // Test that vector did not change
  BOOST_CHECK_EQUAL_COLLECTIONS(vecRes.begin(), vecRes.end(), resultVa.begin(),
                                resultVa.end());
}

BOOST_AUTO_TEST_CASE(AuctioneerTest_FirstValidAuctioneer) {
  // Build arbitrary vector
  std::array<int, 4> vecArb = {0, 1, -2, 4};
  // Let it run through auction
  detail::FirstValidAuctioneer fva;
  std::array<bool, 4> resultFva = fva(vecArb);
  std::array<bool, 4> expected = {false, true, false, false};
  // Test that vector did not change
  BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                resultFva.begin(), resultFva.end());
}

BOOST_AUTO_TEST_CASE(AuctioneerTest_HighestValidAuctioneer) {
  // Build arbitrary vector
  std::array<int, 4> vecArb = {0, 1, -2, 4};
  // Let it run through auction
  detail::HighestValidAuctioneer fva;
  std::array<bool, 4> resultFva = fva(vecArb);
  std::array<bool, 4> expected = {false, false, false, true};
  // Test that vector did not change
  BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                resultFva.begin(), resultFva.end());
}
}  // namespace Acts::Test
