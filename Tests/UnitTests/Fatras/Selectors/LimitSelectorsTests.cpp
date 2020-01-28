// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/MaterialProperties.hpp"
#include "ActsFatras/Selectors/LimitSelectors.hpp"
#include "Dataset.hpp"

using namespace Acts::UnitLiterals;

BOOST_AUTO_TEST_SUITE(LimitSelectors)

BOOST_AUTO_TEST_CASE(X0L0Limits) {
  ActsFatras::X0Limit selectX0;
  ActsFatras::L0Limit selectL0;

  // create particle and move it close to the X0/L0 limit
  auto particle = Dataset::centralPion;
  particle.setLimits(0.15, 0.45);
  particle.update(particle.position(), particle.momentum(), 0.10, 0.34);

  // particle is still within limits for thin block
  BOOST_TEST(not selectX0(Dataset::thinSlab, particle));
  BOOST_TEST(not selectL0(Dataset::thinSlab, particle));
  // particle would pass limits for thick block
  BOOST_TEST(selectX0(Dataset::thickSlab, particle));
  BOOST_TEST(selectL0(Dataset::thickSlab, particle));
}

BOOST_AUTO_TEST_SUITE_END()
