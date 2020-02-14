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

// Construct a particle that is close to its X0/L0 path limit.
//
// Passing a thin slab should still be Ok, but the thick slab should not.
static ActsFatras::Particle makeParticleCloseToLimit() {
  // create particle and move it close to the X0/L0 limit
  auto particle = Dataset::centralPion;
  particle.addPassedMaterial(0.125, 0.0125);
  particle.setMaterialLimits(
      0.125 + 1.125 * Dataset::thinSlab.thicknessInX0(),
      0.0125 + 1.125 * Dataset::thinSlab.thicknessInL0());
  return particle;
}

BOOST_AUTO_TEST_CASE(PathLimitX0) {
  ActsFatras::PathLimitX0 select;
  auto particle = makeParticleCloseToLimit();
  // particle is still within limits for thin block
  BOOST_TEST(not select(Dataset::thinSlab, particle));
  // particle would pass limits for thick block
  BOOST_TEST(select(Dataset::thickSlab, particle));
}

BOOST_AUTO_TEST_CASE(PathLimitL0) {
  ActsFatras::PathLimitL0 select;
  auto particle = makeParticleCloseToLimit();
  // particle is still within limits for thin block
  BOOST_TEST(not select(Dataset::thinSlab, particle));
  // particle would pass limits for thick block
  BOOST_TEST(select(Dataset::thickSlab, particle));
}

BOOST_AUTO_TEST_SUITE_END()
