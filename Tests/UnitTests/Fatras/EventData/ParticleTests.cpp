// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <random>

#include "Acts/Utilities/Units.hpp"
#include "ActsFatras/EventData/Particle.hpp"

using Acts::PdgParticle;
using ActsFatras::Barcode;
using ActsFatras::Particle;
using namespace Acts::UnitLiterals;

BOOST_AUTO_TEST_SUITE(FatrasParticle)

BOOST_AUTO_TEST_CASE(Construct) {
  const auto id = Barcode().setVertexPrimary(1).setParticle(42);
  const auto particle = Particle(id, PdgParticle::eProton, 1_GeV, 1_e);

  BOOST_TEST(particle.id() == id);
  BOOST_TEST(particle.pdg() == PdgParticle::eProton);
  // particle is at rest at the origin
  BOOST_TEST(particle.position4() == Particle::Vector4::Zero());
  BOOST_TEST(particle.position() == Particle::Vector3::Zero());
  BOOST_TEST(particle.time() == Particle::Scalar(0));
  BOOST_TEST(particle.position4().x() == particle.position().x());
  BOOST_TEST(particle.position4().y() == particle.position().y());
  BOOST_TEST(particle.position4().z() == particle.position().z());
  BOOST_TEST(particle.position4().w() == particle.time());
  // particle direction is undefined
  BOOST_TEST(particle.momentum() == Particle::Scalar(0));
  // particle is created at rest and thus not alive
  BOOST_TEST(not particle);
}

BOOST_AUTO_TEST_CASE(CorrectEnergy) {
  const auto id = Barcode().setVertexPrimary(1).setParticle(42);
  auto particle = Particle(id, PdgParticle::eProton, 1_GeV, 1_e)
                      .setDirection(Particle::Vector3::UnitX())
                      .setMomentum(2_GeV);

  BOOST_TEST(particle.mass() == 1_GeV);
  // check that the particle has some input energy
  BOOST_TEST(particle.momentum4().x() == 2_GeV);
  BOOST_TEST(particle.momentum4().y() == 0_GeV);
  BOOST_TEST(particle.momentum4().z() == 0_GeV);
  BOOST_TEST(particle.momentum4().w() == std::hypot(1_GeV, 2_GeV));
  BOOST_TEST(particle.momentum() == 2_GeV);
  BOOST_TEST(particle.energy() == std::hypot(1_GeV, 2_GeV));
  // loose some energy
  particle.correctEnergy(-100_MeV);
  BOOST_TEST(particle.momentum() < 2_GeV);
  BOOST_TEST(particle.energy() ==
             Particle::Scalar(std::hypot(1_GeV, 2_GeV) - 100_MeV));
  // particle is still alive
  BOOST_TEST(particle);
  // loose a lot of energy
  particle.correctEnergy(-3_GeV);
  BOOST_TEST(particle.momentum() == Particle::Scalar(0));
  BOOST_TEST(particle.energy() == particle.mass());
  // lossing even more energy does nothing
  particle.correctEnergy(-10_GeV);
  BOOST_TEST(particle.momentum() == Particle::Scalar(0));
  BOOST_TEST(particle.energy() == particle.mass());
  // particle is not alive anymore
  BOOST_TEST(not particle);
}

BOOST_AUTO_TEST_SUITE_END()
