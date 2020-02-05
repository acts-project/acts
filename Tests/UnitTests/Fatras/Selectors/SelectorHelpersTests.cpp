// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsFatras/Selectors/KinematicCasts.hpp"
#include "ActsFatras/Selectors/SelectorHelpers.hpp"
#include "Dataset.hpp"

namespace {
const auto& backward = Dataset::backwardPion;
const auto& central = Dataset::centralPion;
const auto& forward = Dataset::forwardPion;
}  // namespace

BOOST_AUTO_TEST_SUITE(FatrasSelectorHelpers)

BOOST_AUTO_TEST_CASE(Min) {
  // require a minimum eta value of 0.5
  ActsFatras::Min<ActsFatras::Casts::Eta> minEta{0.5};
  BOOST_TEST(not minEta(backward));
  BOOST_TEST(not minEta(central));
  BOOST_TEST(minEta(forward));

  // require a mininum absolute eta value of 0.5
  ActsFatras::Min<ActsFatras::Casts::AbsEta> minAbsEta{0.5};
  BOOST_TEST(minAbsEta(backward));
  BOOST_TEST(not minAbsEta(central));
  BOOST_TEST(minAbsEta(forward));
}

BOOST_AUTO_TEST_CASE(Max) {
  // require a maximum eta value of 0.5
  ActsFatras::Max<ActsFatras::Casts::Eta> maxEta{0.5};
  BOOST_TEST(maxEta(backward));
  BOOST_TEST(maxEta(central));
  BOOST_TEST(not maxEta(forward));

  // require a maximum absolute eta value of 0.5
  ActsFatras::Max<ActsFatras::Casts::AbsEta> maxAbsEta{0.5};
  BOOST_TEST(not maxAbsEta(backward));
  BOOST_TEST(maxAbsEta(central));
  BOOST_TEST(not maxAbsEta(forward));
}

BOOST_AUTO_TEST_CASE(Range) {
  ActsFatras::Range<ActsFatras::Casts::Eta> rangeEta{-6.0, -0.5};
  BOOST_TEST(rangeEta(backward));
  BOOST_TEST(not rangeEta(central));
  BOOST_TEST(not rangeEta(forward));

  ActsFatras::Range<ActsFatras::Casts::AbsEta> rangeAbsEta{0.5, 6.0};
  BOOST_TEST(rangeAbsEta(backward));
  BOOST_TEST(not rangeAbsEta(central));
  BOOST_TEST(rangeAbsEta(forward));
}

BOOST_AUTO_TEST_SUITE_END()
