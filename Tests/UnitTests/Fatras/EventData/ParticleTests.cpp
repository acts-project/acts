// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <limits>

using Acts::PdgParticle;
using ActsFatras::Barcode;
using ActsFatras::Particle;

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace ActsFatras;

namespace {
constexpr auto eps = std::numeric_limits<double>::epsilon();
}

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EventDataSuite)

BOOST_AUTO_TEST_CASE(Construct) {
  const auto pid = Barcode().withVertexPrimary(1).withParticle(42);
  const auto particle = Particle(pid, PdgParticle::eProton, 1_e, 1_GeV);

  BOOST_CHECK_EQUAL(particle.particleId(), pid);
  BOOST_CHECK_EQUAL(particle.pdg(), PdgParticle::eProton);
  // particle is at rest at the origin
  BOOST_CHECK_EQUAL(particle.fourPosition(), Vector4::Zero());
  BOOST_CHECK_EQUAL(particle.position(), Vector3::Zero());
  BOOST_CHECK_EQUAL(particle.time(), 0.);
  BOOST_CHECK_EQUAL(particle.fourPosition().x(), particle.position().x());
  BOOST_CHECK_EQUAL(particle.fourPosition().y(), particle.position().y());
  BOOST_CHECK_EQUAL(particle.fourPosition().z(), particle.position().z());
  BOOST_CHECK_EQUAL(particle.fourPosition().w(), particle.time());
  // particle direction is undefined, but must be normalized
  CHECK_CLOSE_REL(particle.direction().norm(), 1, eps);
  BOOST_CHECK_EQUAL(particle.transverseMomentum(), 0.);
  BOOST_CHECK_EQUAL(particle.absoluteMomentum(), 0.);
}

BOOST_AUTO_TEST_CASE(CorrectEnergy) {
  const auto pid = Barcode().withVertexPrimary(1).withParticle(42);
  auto particle = Particle(pid, PdgParticle::eProton, 1_e, 1_GeV)
                      .setDirection(Vector3::UnitX())
                      .setAbsoluteMomentum(2_GeV);

  BOOST_CHECK_EQUAL(particle.mass(), 1_GeV);
  // check that the particle has some input energy
  BOOST_CHECK_EQUAL(particle.fourMomentum().x(), 2_GeV);
  BOOST_CHECK_EQUAL(particle.fourMomentum().y(), 0_GeV);
  BOOST_CHECK_EQUAL(particle.fourMomentum().z(), 0_GeV);
  BOOST_CHECK_EQUAL(particle.fourMomentum().w(), std::hypot(1_GeV, 2_GeV));
  BOOST_CHECK_EQUAL(particle.transverseMomentum(), 2_GeV);
  BOOST_CHECK_EQUAL(particle.absoluteMomentum(), 2_GeV);
  BOOST_CHECK_EQUAL(particle.energy(), std::hypot(1_GeV, 2_GeV));
  // particle direction must be normalized
  CHECK_CLOSE_REL(particle.direction().norm(), 1, eps);

  // lose some energy
  particle.loseEnergy(100_MeV);
  BOOST_CHECK_LT(particle.transverseMomentum(), 2_GeV);
  BOOST_CHECK_LT(particle.absoluteMomentum(), 2_GeV);
  CHECK_CLOSE_REL(particle.energy(), std::hypot(1_GeV, 2_GeV) - 100_MeV, eps);
  CHECK_CLOSE_REL(particle.direction().norm(), 1, eps);

  // lose some more energy
  particle.loseEnergy(200_MeV);
  BOOST_CHECK_LT(particle.transverseMomentum(), 2_GeV);
  BOOST_CHECK_LT(particle.absoluteMomentum(), 2_GeV);
  CHECK_CLOSE_REL(particle.energy(), std::hypot(1_GeV, 2_GeV) - 300_MeV, eps);
  CHECK_CLOSE_REL(particle.direction().norm(), 1, eps);

  // lose a lot of energy
  particle.loseEnergy(3_GeV);
  BOOST_CHECK_EQUAL(particle.transverseMomentum(), 0.);
  BOOST_CHECK_EQUAL(particle.absoluteMomentum(), 0.);
  BOOST_CHECK_EQUAL(particle.energy(), particle.mass());
  CHECK_CLOSE_REL(particle.direction().norm(), 1, eps);

  // losing even more energy does nothing
  BOOST_CHECK_THROW(particle.loseEnergy(10_GeV), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
