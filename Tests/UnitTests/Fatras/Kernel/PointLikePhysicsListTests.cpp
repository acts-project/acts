// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Kernel/PointLikePhysicsList.hpp"

#include <limits>
#include <random>
#include <utility>
#include <vector>

using namespace Acts::UnitLiterals;
using Particle = ActsFatras::Particle;
using Scalar = ActsFatras::Particle::Scalar;
using Rng = std::ranlux48;

namespace {

/// EM-like physics process that triggers on X0 and leaves the particle alive.
///
/// Each run call creates one descendant particle.
struct X0Process {
  std::pair<Scalar, Scalar> generatePathLimits(Rng &, const Particle &) const {
    return {0.5, std::numeric_limits<Scalar>::infinity()};
  }

  bool run(Rng &, Particle &particle, std::vector<Particle> &generated) const {
    auto pid0 = particle.particleId().makeDescendant(0);
    generated.emplace_back(particle.withParticleId(pid0));
    return false;
  }
};

/// Nuclear-like physics process that triggers on L0 and destroys the particle.
///
/// Each run call creates two descendant particles.
struct L0Process {
  std::pair<Scalar, Scalar> generatePathLimits(Rng &, const Particle &) const {
    return {std::numeric_limits<Scalar>::infinity(), 1.5};
  }

  bool run(Rng &, Particle &particle, std::vector<Particle> &generated) const {
    auto pid0 = particle.particleId().makeDescendant(0);
    auto pid1 = particle.particleId().makeDescendant(1);
    generated.emplace_back(particle.withParticleId(pid0));
    generated.emplace_back(particle.withParticleId(pid1));
    return true;
  }
};

struct Fixture {
  Rng rng{23};
  Particle incoming;
  std::vector<Particle> outgoing;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(FatrasPointLikePhysicsList)

BOOST_AUTO_TEST_CASE(Empty) {
  Fixture f;
  ActsFatras::PointLikePhysicsList<> pl;

  // w/o processes the list should never abort
  auto sel = pl.arm(f.rng, f.incoming);
  BOOST_CHECK_EQUAL(sel.x0Limit, std::numeric_limits<Scalar>::infinity());
  BOOST_CHECK_EQUAL(sel.l0Limit, std::numeric_limits<Scalar>::infinity());
  BOOST_CHECK_EQUAL(sel.x0Process, SIZE_MAX);
  BOOST_CHECK_EQUAL(sel.l0Process, SIZE_MAX);

  // running with an invalid process index should do nothing
  // physics lists is empty and 0 should already be invalid
  BOOST_CHECK(not pl.run(f.rng, 0u, f.incoming, f.outgoing));
  BOOST_CHECK_EQUAL(f.outgoing.size(), 0u);
  // SIZE_MAX should always be an invalid index
  BOOST_CHECK(not pl.run(f.rng, SIZE_MAX, f.incoming, f.outgoing));
  BOOST_CHECK_EQUAL(f.outgoing.size(), 0u);
}

BOOST_AUTO_TEST_CASE(SingleX0) {
  Fixture f;
  ActsFatras::PointLikePhysicsList<X0Process> pl;

  // w/o processes the list should never abort
  auto sel = pl.arm(f.rng, f.incoming);
  BOOST_CHECK_EQUAL(sel.x0Limit, 0.5);
  BOOST_CHECK_EQUAL(sel.l0Limit, std::numeric_limits<Scalar>::infinity());
  BOOST_CHECK_EQUAL(sel.x0Process, 0u);
  BOOST_CHECK_EQUAL(sel.l0Process, SIZE_MAX);

  // valid index, X0Process leaves the particle alive
  BOOST_CHECK(not pl.run(f.rng, 0u, f.incoming, f.outgoing));
  BOOST_CHECK_EQUAL(f.outgoing.size(), 1u);
  // invalid index, should do nothing
  BOOST_CHECK(not pl.run(f.rng, SIZE_MAX, f.incoming, f.outgoing));
  BOOST_CHECK_EQUAL(f.outgoing.size(), 1u);
}

BOOST_AUTO_TEST_CASE(SingleL0) {
  Fixture f;
  ActsFatras::PointLikePhysicsList<L0Process> pl;

  // w/o processes the list should never abort
  auto sel = pl.arm(f.rng, f.incoming);
  BOOST_CHECK_EQUAL(sel.x0Limit, std::numeric_limits<Scalar>::infinity());
  BOOST_CHECK_EQUAL(sel.l0Limit, 1.5);
  BOOST_CHECK_EQUAL(sel.x0Process, SIZE_MAX);
  BOOST_CHECK_EQUAL(sel.l0Process, 0u);

  // valid index, L0Process kills the particles and creates 2 descendants
  BOOST_CHECK(pl.run(f.rng, 0u, f.incoming, f.outgoing));
  BOOST_CHECK_EQUAL(f.outgoing.size(), 2u);
  // invalid index, should do nothing
  BOOST_CHECK(not pl.run(f.rng, SIZE_MAX, f.incoming, f.outgoing));
  BOOST_CHECK_EQUAL(f.outgoing.size(), 2u);
}

BOOST_AUTO_TEST_CASE(BothX0L0) {
  Fixture f;
  ActsFatras::PointLikePhysicsList<X0Process, L0Process> pl;

  // w/o processes the list should never abort
  auto sel = pl.arm(f.rng, f.incoming);
  BOOST_CHECK_EQUAL(sel.x0Limit, 0.5);
  BOOST_CHECK_EQUAL(sel.l0Limit, 1.5);
  BOOST_CHECK_EQUAL(sel.x0Process, 0u);
  BOOST_CHECK_EQUAL(sel.l0Process, 1u);

  // valid index, X0Process leaves the particle alive
  BOOST_CHECK(not pl.run(f.rng, 0u, f.incoming, f.outgoing));
  BOOST_CHECK_EQUAL(f.outgoing.size(), 1u);
  // valid index, L0Process kills the particles and creates 2 descendants
  BOOST_CHECK(pl.run(f.rng, 1u, f.incoming, f.outgoing));
  BOOST_CHECK_EQUAL(f.outgoing.size(), 3u);
  // invalid index, should do nothing
  BOOST_CHECK(not pl.run(f.rng, SIZE_MAX, f.incoming, f.outgoing));
  BOOST_CHECK_EQUAL(f.outgoing.size(), 3u);
}

BOOST_AUTO_TEST_CASE(Disable) {
  Fixture f;
  ActsFatras::PointLikePhysicsList<X0Process, L0Process> pl;

  // disabled X0Process should participate in arming procedure
  pl.disable<X0Process>();
  {
    auto sel = pl.arm(f.rng, f.incoming);
    BOOST_CHECK_EQUAL(sel.x0Limit, std::numeric_limits<Scalar>::infinity());
    BOOST_CHECK_EQUAL(sel.l0Limit, 1.5);
    BOOST_CHECK_EQUAL(sel.x0Process, SIZE_MAX);
    BOOST_CHECK_EQUAL(sel.l0Process, 1u);
  }

  // disabling L0Process is equivalent to an empty list
  pl.disable<L0Process>();
  {
    auto sel = pl.arm(f.rng, f.incoming);
    BOOST_CHECK_EQUAL(sel.x0Limit, std::numeric_limits<Scalar>::infinity());
    BOOST_CHECK_EQUAL(sel.l0Limit, std::numeric_limits<Scalar>::infinity());
    BOOST_CHECK_EQUAL(sel.x0Process, SIZE_MAX);
    BOOST_CHECK_EQUAL(sel.l0Process, SIZE_MAX);
  }

  // valid index for X0Process, should do nothing since its disabled
  BOOST_CHECK(not pl.run(f.rng, 0u, f.incoming, f.outgoing));
  BOOST_CHECK_EQUAL(f.outgoing.size(), 0u);
  // valid index for L0Process, should do nothing since its disables
  BOOST_CHECK(not pl.run(f.rng, 1u, f.incoming, f.outgoing));
  BOOST_CHECK_EQUAL(f.outgoing.size(), 0u);
  // invalid index, should do nothing
  BOOST_CHECK(not pl.run(f.rng, SIZE_MAX, f.incoming, f.outgoing));
  BOOST_CHECK_EQUAL(f.outgoing.size(), 0u);
}

BOOST_AUTO_TEST_SUITE_END()
