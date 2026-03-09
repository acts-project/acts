// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/Charge.hpp"
#include "Acts/Utilities/Diagnostics.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <limits>
#include <type_traits>

using namespace Acts;
using namespace Acts::UnitLiterals;

static auto eps = std::numeric_limits<double>::epsilon();

ACTS_PUSH_IGNORE_DEPRECATED()
BOOST_TEST_DONT_PRINT_LOG_VALUE(Neutral)
BOOST_TEST_DONT_PRINT_LOG_VALUE(SinglyCharged)
BOOST_TEST_DONT_PRINT_LOG_VALUE(NonNeutralCharge)
ACTS_POP_IGNORE_DEPRECATED()
BOOST_TEST_DONT_PRINT_LOG_VALUE(AnyCharge)

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EventDataSuite)

ACTS_PUSH_IGNORE_DEPRECATED()

BOOST_AUTO_TEST_CASE(Constructibility) {
  BOOST_CHECK(std::is_trivially_default_constructible_v<Neutral>);
  BOOST_CHECK(std::is_trivially_default_constructible_v<SinglyCharged>);
  BOOST_CHECK(std::is_nothrow_default_constructible_v<Neutral>);
  BOOST_CHECK(std::is_nothrow_default_constructible_v<SinglyCharged>);
  BOOST_CHECK(std::is_trivially_constructible_v<Neutral>);
  BOOST_CHECK(std::is_trivially_constructible_v<SinglyCharged>);
  // BOOST_CHECK(std::is_trivially_constructible_v<NonNeutralCharge>);
  // BOOST_CHECK(std::is_trivially_constructible_v<AnyCharge>);
  BOOST_CHECK(std::is_nothrow_constructible_v<Neutral>);
  BOOST_CHECK(std::is_nothrow_constructible_v<SinglyCharged>);
  // BOOST_CHECK(std::is_nothrow_constructible_v<NonNeutralCharge>);
  // BOOST_CHECK(std::is_nothrow_constructible_v<AnyCharge>);
}

BOOST_AUTO_TEST_CASE(NeutralTest) {
  Neutral q;

  BOOST_CHECK_EQUAL(q.extractCharge(1.23), 0_e);
  BOOST_CHECK_EQUAL(q.extractCharge(2.54), 0_e);
  BOOST_CHECK_EQUAL(q.extractCharge(-1.98), 0_e);
  BOOST_CHECK_EQUAL(q.extractCharge(-2.23), 0_e);
  CHECK_CLOSE_REL(q.extractMomentum(1 / 64_GeV), 64_GeV, eps);
  CHECK_CLOSE_REL(q.extractMomentum(1 / 128_MeV), 128_MeV, eps);

  // negative inputs should not occur for neutral particles
  // the result is not defined, but we check it anyway
  // update: this is asserted now
  // CHECK_CLOSE_REL(q.extractMomentum(-1 / 128_MeV), -128_MeV, eps);

  BOOST_CHECK_EQUAL(q, Neutral());
  BOOST_CHECK_EQUAL(Neutral(), q);
  BOOST_CHECK_EQUAL(q, Neutral(0_e));
  BOOST_CHECK_EQUAL(Neutral(0_e), q);
}

BOOST_AUTO_TEST_CASE(SinglyChargedTest) {
  SinglyCharged q;

  BOOST_CHECK_EQUAL(q.extractCharge(1.23), 1_e);
  BOOST_CHECK_EQUAL(q.extractCharge(2.54), 1_e);
  BOOST_CHECK_EQUAL(q.extractCharge(-1.98), -1_e);
  BOOST_CHECK_EQUAL(q.extractCharge(-2.23), -1_e);
  CHECK_CLOSE_REL(q.extractMomentum(1_e / 64_GeV), 64_GeV, eps);
  CHECK_CLOSE_REL(q.extractMomentum(1_e / 128_MeV), 128_MeV, eps);
  CHECK_CLOSE_REL(q.extractMomentum(-1_e / 128_MeV), 128_MeV, eps);

  BOOST_CHECK_EQUAL(q, SinglyCharged());
  BOOST_CHECK_EQUAL(SinglyCharged(), q);
  BOOST_CHECK_EQUAL(q, SinglyCharged(1_e));
  BOOST_CHECK_EQUAL(SinglyCharged(1_e), q);
}

BOOST_AUTO_TEST_CASE(NonNeutralChargeSingle) {
  NonNeutralCharge q(1_e);

  BOOST_CHECK_EQUAL(q.extractCharge(1.23), 1_e);
  BOOST_CHECK_EQUAL(q.extractCharge(2.54), 1_e);
  BOOST_CHECK_EQUAL(q.extractCharge(-1.98), -1_e);
  BOOST_CHECK_EQUAL(q.extractCharge(-2.23), -1_e);
  CHECK_CLOSE_REL(q.extractMomentum(1_e / 64_GeV), 64_GeV, eps);
  CHECK_CLOSE_REL(q.extractMomentum(1_e / 128_MeV), 128_MeV, eps);
  CHECK_CLOSE_REL(q.extractMomentum(-1_e / 128_MeV), 128_MeV, eps);

  BOOST_CHECK_EQUAL(q, NonNeutralCharge(1_e));
  BOOST_CHECK_EQUAL(NonNeutralCharge(1_e), q);
}

BOOST_AUTO_TEST_CASE(NonNeutralChargeMultiple) {
  NonNeutralCharge q(3_e);

  BOOST_CHECK_EQUAL(q.extractCharge(1.23), 3_e);
  BOOST_CHECK_EQUAL(q.extractCharge(2.54), 3_e);
  BOOST_CHECK_EQUAL(q.extractCharge(-1.98), -3_e);
  BOOST_CHECK_EQUAL(q.extractCharge(-2.23), -3_e);
  CHECK_CLOSE_REL(q.extractMomentum(3_e / 64_GeV), 64_GeV, eps);
  CHECK_CLOSE_REL(q.extractMomentum(3_e / 128_MeV), 128_MeV, eps);
  CHECK_CLOSE_REL(q.extractMomentum(-3_e / 128_MeV), 128_MeV, eps);

  BOOST_CHECK(!(q == NonNeutralCharge(1_e)));
  BOOST_CHECK(!(NonNeutralCharge(1_e) == q));
  BOOST_CHECK(!(q == NonNeutralCharge(2_e)));
  BOOST_CHECK(!(NonNeutralCharge(2_e) == q));
  BOOST_CHECK(q == NonNeutralCharge(3_e));
  BOOST_CHECK(NonNeutralCharge(3_e) == q);
}

ACTS_POP_IGNORE_DEPRECATED()

BOOST_AUTO_TEST_CASE(AnyChargeNeutral) {
  AnyCharge q(0_e);

  BOOST_CHECK_EQUAL(q.extractCharge(1.23), 0_e);
  BOOST_CHECK_EQUAL(q.extractCharge(2.54), 0_e);
  BOOST_CHECK_EQUAL(q.extractCharge(-1.98), 0_e);
  BOOST_CHECK_EQUAL(q.extractCharge(-2.23), 0_e);
  CHECK_CLOSE_REL(q.extractMomentum(1 / 64_GeV), 64_GeV, eps);
  CHECK_CLOSE_REL(q.extractMomentum(1 / 128_MeV), 128_MeV, eps);

  // negative inputs should not occur for neutral particles
  // the result is not defined, but we check it anyway
  CHECK_CLOSE_REL(q.extractMomentum(-1 / 128_MeV), -128_MeV, eps);

  BOOST_CHECK(q == AnyCharge(0_e));
  BOOST_CHECK(AnyCharge(0_e) == q);
  BOOST_CHECK(!(q == AnyCharge(1_e)));
  BOOST_CHECK(!(AnyCharge(1_e) == q));
  BOOST_CHECK(!(q == AnyCharge(2_e)));
  BOOST_CHECK(!(AnyCharge(2_e) == q));
}

BOOST_AUTO_TEST_CASE(AnyChargeSingle) {
  AnyCharge q(1_e);

  BOOST_CHECK_EQUAL(q.extractCharge(1.23), 1_e);
  BOOST_CHECK_EQUAL(q.extractCharge(2.54), 1_e);
  BOOST_CHECK_EQUAL(q.extractCharge(-1.98), -1_e);
  BOOST_CHECK_EQUAL(q.extractCharge(-2.23), -1_e);
  CHECK_CLOSE_REL(q.extractMomentum(1_e / 64_GeV), 64_GeV, eps);
  CHECK_CLOSE_REL(q.extractMomentum(1_e / 128_MeV), 128_MeV, eps);
  CHECK_CLOSE_REL(q.extractMomentum(-1_e / 128_MeV), 128_MeV, eps);

  BOOST_CHECK(!(q == AnyCharge(0_e)));
  BOOST_CHECK(!(AnyCharge(0_e) == q));
  BOOST_CHECK(q == AnyCharge(1_e));
  BOOST_CHECK(AnyCharge(1_e) == q);
  BOOST_CHECK(!(q == AnyCharge(2_e)));
  BOOST_CHECK(!(AnyCharge(2_e) == q));
}

BOOST_AUTO_TEST_CASE(AnyChargeMultiple) {
  AnyCharge q(3_e);

  BOOST_CHECK_EQUAL(q.extractCharge(1.23), 3_e);
  BOOST_CHECK_EQUAL(q.extractCharge(2.54), 3_e);
  BOOST_CHECK_EQUAL(q.extractCharge(-1.98), -3_e);
  BOOST_CHECK_EQUAL(q.extractCharge(-2.23), -3_e);
  CHECK_CLOSE_REL(q.extractMomentum(3_e / 64_GeV), 64_GeV, eps);
  CHECK_CLOSE_REL(q.extractMomentum(3_e / 128_MeV), 128_MeV, eps);
  CHECK_CLOSE_REL(q.extractMomentum(-3_e / 128_MeV), 128_MeV, eps);

  BOOST_CHECK(!(q == AnyCharge(0_e)));
  BOOST_CHECK(!(AnyCharge(0_e) == q));
  BOOST_CHECK(!(q == AnyCharge(1_e)));
  BOOST_CHECK(!(AnyCharge(1_e) == q));
  BOOST_CHECK(!(q == AnyCharge(2_e)));
  BOOST_CHECK(!(AnyCharge(2_e) == q));
  BOOST_CHECK(q == AnyCharge(3_e));
  BOOST_CHECK(AnyCharge(3_e) == q);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
