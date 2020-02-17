// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "ActsFatras/Selectors/LimitSelectors.hpp"
#include "Dataset.hpp"

namespace {
// Construct a particle that is close to its X0/L0 path limit.
//
// Passing a thin slab should still be Ok, but the thick slab should not.
ActsFatras::Particle makeParticleCloseToLimit() {
  // create particle and move it close to the X0/L0 limit
  auto particle = Dataset::centralPion;
  particle.setMaterialPassed(0.125, 0.0125);
  // limit is a bit above 1% in radiation/interaction length
  particle.setMaterialLimits(0.125 + 0.0125, 0.0125 + 0.0125);
  return particle;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(LimitSelectors)

BOOST_AUTO_TEST_CASE(PathLimitX0) {
  ActsFatras::PathLimitX0 select;
  auto particle = makeParticleCloseToLimit();
  // particle is still within limits for thin block
  BOOST_TEST(not select(particle, Acts::Test::makePercentSlab()));
  // particle would pass limits for thick block
  BOOST_TEST(select(particle, Acts::Test::makeUnitSlab()));
}

BOOST_AUTO_TEST_CASE(PathLimitL0) {
  ActsFatras::PathLimitL0 select;
  auto particle = makeParticleCloseToLimit();
  // particle is still within limits for thin block
  BOOST_TEST(not select(particle, Acts::Test::makePercentSlab()));
  // particle would pass limits for thick block
  BOOST_TEST(select(particle, Acts::Test::makeUnitSlab()));
}

BOOST_AUTO_TEST_SUITE_END()
