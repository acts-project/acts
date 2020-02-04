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
const auto& detector = Dataset::thinSlab;
const auto& backward = Dataset::backwardPion;
const auto& central = Dataset::centralPion;
const auto& forward = Dataset::forwardPion;
}  // namespace

BOOST_AUTO_TEST_SUITE(FatrasSelectorHelpers)

BOOST_AUTO_TEST_CASE(Min) {
  // require a minimum eta value of 0.5
  ActsFatras::Min<ActsFatras::Casts::Eta> minEta;
  minEta.valMin = 0.5;

  BOOST_TEST(not minEta(detector, backward));
  BOOST_TEST(not minEta(detector, central));
  BOOST_TEST(minEta(detector, forward));

  // require a mininum absolute eta value of 0.5
  ActsFatras::Min<ActsFatras::Casts::AbsEta> minAbsEta;
  minAbsEta.valMin = 0.5;

  BOOST_TEST(minAbsEta(detector, backward));
  BOOST_TEST(not minAbsEta(detector, central));
  BOOST_TEST(minAbsEta(detector, forward));
}

BOOST_AUTO_TEST_CASE(Max) {
  // require a maximum eta value of 0.5
  ActsFatras::Max<ActsFatras::Casts::Eta> maxEta;
  maxEta.valMax = 0.5;

  BOOST_TEST(maxEta(detector, backward));
  BOOST_TEST(maxEta(detector, central));
  BOOST_TEST(not maxEta(detector, forward));

  // require a maximum absolute eta value of 0.5
  ActsFatras::Max<ActsFatras::Casts::AbsEta> maxAbsEta;
  maxAbsEta.valMax = 0.5;

  BOOST_TEST(not maxAbsEta(detector, backward));
  BOOST_TEST(maxAbsEta(detector, central));
  BOOST_TEST(not maxAbsEta(detector, forward));
}

BOOST_AUTO_TEST_CASE(Range) {
  ActsFatras::Range<ActsFatras::Casts::Eta> rangeEta;
  rangeEta.valMin = -6.;
  rangeEta.valMax = -0.5;

  BOOST_TEST(rangeEta(detector, backward));
  BOOST_TEST(not rangeEta(detector, central));
  BOOST_TEST(not rangeEta(detector, forward));

  ActsFatras::Range<ActsFatras::Casts::AbsEta> rangeAbsEta;
  rangeAbsEta.valMin = 0.5;
  rangeAbsEta.valMax = 6.0;

  BOOST_TEST(rangeAbsEta(detector, backward));
  BOOST_TEST(not rangeAbsEta(detector, central));
  BOOST_TEST(rangeAbsEta(detector, forward));
}

BOOST_AUTO_TEST_SUITE_END()
