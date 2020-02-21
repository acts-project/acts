// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <random>

#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Physics/EnergyLoss/BetheBloch.hpp"
#include "ActsFatras/Physics/EnergyLoss/BetheHeitler.hpp"
#include "Dataset.hpp"

using Generator = std::ranlux48;

BOOST_AUTO_TEST_SUITE(FatrasEnergyLoss)

BOOST_DATA_TEST_CASE(BetheBloch, Dataset::parameters, pdg, phi, lambda, p,
                     seed) {
  Generator gen(seed);
  ActsFatras::Particle before = Dataset::makeParticle(pdg, phi, lambda, p);
  ActsFatras::Particle after = before;

  ActsFatras::BetheBloch process;
  const auto outgoing = process(gen, Acts::Test::makeUnitSlab(), after);
  // energy loss changes momentum and energy
  BOOST_TEST(after.absMomentum() < before.absMomentum());
  BOOST_TEST(after.energy() < before.energy());
  // energy loss creates no new particles
  BOOST_TEST(outgoing.empty());
}

BOOST_DATA_TEST_CASE(BetheHeitler, Dataset::parameters, pdg, phi, lambda, p,
                     seed) {
  Generator gen(seed);
  ActsFatras::Particle before = Dataset::makeParticle(pdg, phi, lambda, p);
  ActsFatras::Particle after = before;

  ActsFatras::BetheHeitler process;
  const auto outgoing = process(gen, Acts::Test::makeUnitSlab(), after);
  // energy loss changes momentum and energy
  BOOST_TEST(after.absMomentum() < before.absMomentum());
  BOOST_TEST(after.energy() < before.energy());
  // energy loss creates no new particles
  BOOST_TEST(outgoing.empty());
}

BOOST_AUTO_TEST_SUITE_END()
