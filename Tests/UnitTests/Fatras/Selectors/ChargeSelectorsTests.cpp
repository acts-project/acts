// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsFatras/Selectors/ChargeSelectors.hpp"
#include "Dataset.hpp"

using namespace ActsFatras;

BOOST_AUTO_TEST_SUITE(FatrasChargeSelectors)

BOOST_AUTO_TEST_CASE(NegativeParticle) {
  const auto& particle = Dataset::centralElectron;
  BOOST_TEST(not NeutralSelector()(particle));
  BOOST_TEST(ChargedSelector()(particle));
  BOOST_TEST(not PositiveSelector()(particle));
  BOOST_TEST(NegativeSelector()(particle));
}

BOOST_AUTO_TEST_CASE(NeutralParticle) {
  const auto& particle = Dataset::centralNeutron;
  BOOST_TEST(NeutralSelector()(particle));
  BOOST_TEST(not ChargedSelector()(particle));
  BOOST_TEST(not PositiveSelector()(particle));
  BOOST_TEST(not NegativeSelector()(particle));
}

BOOST_AUTO_TEST_CASE(PositiveParticle) {
  const auto& particle = Dataset::centralPositron;
  BOOST_TEST(not NeutralSelector()(particle));
  BOOST_TEST(ChargedSelector()(particle));
  BOOST_TEST(PositiveSelector()(particle));
  BOOST_TEST(not NegativeSelector()(particle));
}

BOOST_AUTO_TEST_SUITE_END()
