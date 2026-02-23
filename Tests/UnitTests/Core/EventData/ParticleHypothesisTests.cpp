// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <limits>

using namespace Acts;
using namespace Acts::UnitLiterals;

static auto eps = std::numeric_limits<double>::epsilon();

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EventDataSuite)

ACTS_PUSH_IGNORE_DEPRECATED()

BOOST_AUTO_TEST_CASE(Neutral) {
  auto p = NeutralParticleHypothesis::pion0();

  BOOST_CHECK_EQUAL(p.extractCharge(1.23), 0_e);
  BOOST_CHECK_EQUAL(p.extractCharge(2.54), 0_e);
  BOOST_CHECK_EQUAL(p.extractCharge(-1.98), 0_e);
  BOOST_CHECK_EQUAL(p.extractCharge(-2.23), 0_e);
  CHECK_CLOSE_REL(p.extractMomentum(1 / 64_GeV), 64_GeV, eps);
  CHECK_CLOSE_REL(p.extractMomentum(1 / 128_MeV), 128_MeV, eps);
}

BOOST_AUTO_TEST_CASE(SinglyCharged) {
  auto p = SinglyChargedParticleHypothesis::pion();

  BOOST_CHECK_EQUAL(p.extractCharge(1.23), 1_e);
  BOOST_CHECK_EQUAL(p.extractCharge(2.54), 1_e);
  BOOST_CHECK_EQUAL(p.extractCharge(-1.98), -1_e);
  BOOST_CHECK_EQUAL(p.extractCharge(-2.23), -1_e);
  CHECK_CLOSE_REL(p.extractMomentum(1_e / 64_GeV), 64_GeV, eps);
  CHECK_CLOSE_REL(p.extractMomentum(1_e / 128_MeV), 128_MeV, eps);
  CHECK_CLOSE_REL(p.extractMomentum(-1_e / 128_MeV), 128_MeV, eps);
}

BOOST_AUTO_TEST_CASE(NonNeutralChargeSingle) {
  auto p = NonNeutralChargedParticleHypothesis::pion();

  BOOST_CHECK_EQUAL(p.extractCharge(1.23), 1_e);
  BOOST_CHECK_EQUAL(p.extractCharge(2.54), 1_e);
  BOOST_CHECK_EQUAL(p.extractCharge(-1.98), -1_e);
  BOOST_CHECK_EQUAL(p.extractCharge(-2.23), -1_e);
  CHECK_CLOSE_REL(p.extractMomentum(1_e / 64_GeV), 64_GeV, eps);
  CHECK_CLOSE_REL(p.extractMomentum(1_e / 128_MeV), 128_MeV, eps);
  CHECK_CLOSE_REL(p.extractMomentum(-1_e / 128_MeV), 128_MeV, eps);
}

BOOST_AUTO_TEST_CASE(NonNeutralChargeMultiple) {
  auto p = NonNeutralChargedParticleHypothesis::pionLike(3_e);

  BOOST_CHECK_EQUAL(p.extractCharge(1.23), 3_e);
  BOOST_CHECK_EQUAL(p.extractCharge(2.54), 3_e);
  BOOST_CHECK_EQUAL(p.extractCharge(-1.98), -3_e);
  BOOST_CHECK_EQUAL(p.extractCharge(-2.23), -3_e);
  CHECK_CLOSE_REL(p.extractMomentum(3_e / 64_GeV), 64_GeV, eps);
  CHECK_CLOSE_REL(p.extractMomentum(3_e / 128_MeV), 128_MeV, eps);
  CHECK_CLOSE_REL(p.extractMomentum(-3_e / 128_MeV), 128_MeV, eps);
}

ACTS_POP_IGNORE_DEPRECATED()

BOOST_AUTO_TEST_CASE(AnyChargeNeutral) {
  auto p = ParticleHypothesis::pion0();

  BOOST_CHECK_EQUAL(p.extractCharge(1.23), 0_e);
  BOOST_CHECK_EQUAL(p.extractCharge(2.54), 0_e);
  BOOST_CHECK_EQUAL(p.extractCharge(-1.98), 0_e);
  BOOST_CHECK_EQUAL(p.extractCharge(-2.23), 0_e);
  CHECK_CLOSE_REL(p.extractMomentum(1 / 64_GeV), 64_GeV, eps);
  CHECK_CLOSE_REL(p.extractMomentum(1 / 128_MeV), 128_MeV, eps);

  // negative inputs should not occur for neutral particles
  // the result is not defined, but we check it anyway
  CHECK_CLOSE_REL(p.extractMomentum(-1 / 128_MeV), -128_MeV, eps);
}

BOOST_AUTO_TEST_CASE(AnyChargeSingle) {
  auto p = ParticleHypothesis::pion();

  BOOST_CHECK_EQUAL(p.extractCharge(1.23), 1_e);
  BOOST_CHECK_EQUAL(p.extractCharge(2.54), 1_e);
  BOOST_CHECK_EQUAL(p.extractCharge(-1.98), -1_e);
  BOOST_CHECK_EQUAL(p.extractCharge(-2.23), -1_e);
  CHECK_CLOSE_REL(p.extractMomentum(1_e / 64_GeV), 64_GeV, eps);
  CHECK_CLOSE_REL(p.extractMomentum(1_e / 128_MeV), 128_MeV, eps);
  CHECK_CLOSE_REL(p.extractMomentum(-1_e / 128_MeV), 128_MeV, eps);
}

BOOST_AUTO_TEST_CASE(AnyChargeMultiple) {
  auto p = ParticleHypothesis::pionLike(3_e);

  BOOST_CHECK_EQUAL(p.extractCharge(1.23), 3_e);
  BOOST_CHECK_EQUAL(p.extractCharge(2.54), 3_e);
  BOOST_CHECK_EQUAL(p.extractCharge(-1.98), -3_e);
  BOOST_CHECK_EQUAL(p.extractCharge(-2.23), -3_e);
  CHECK_CLOSE_REL(p.extractMomentum(3_e / 64_GeV), 64_GeV, eps);
  CHECK_CLOSE_REL(p.extractMomentum(3_e / 128_MeV), 128_MeV, eps);
  CHECK_CLOSE_REL(p.extractMomentum(-3_e / 128_MeV), 128_MeV, eps);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
