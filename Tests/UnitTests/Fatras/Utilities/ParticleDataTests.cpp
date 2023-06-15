// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/PdgParticle.hpp"
#include "ActsFatras/Utilities/ParticleData.hpp"

#include <cmath>
#include <string_view>

using Acts::PdgParticle;
using namespace Acts::UnitLiterals;
using namespace ActsFatras;

namespace {
// NOTE: the used mass comparison values are not as exact as the data values
static constexpr float eps = 0.001f;
}  // namespace

BOOST_AUTO_TEST_SUITE(FatrasParticleData)

BOOST_AUTO_TEST_CASE(InvalidInput) {
  BOOST_CHECK(std::isnan(findCharge(PdgParticle::eInvalid)));
  BOOST_CHECK_EQUAL(findMass(PdgParticle::eInvalid), 0.0f);
  BOOST_CHECK(findName(PdgParticle::eInvalid).empty());
}

BOOST_AUTO_TEST_CASE(Electron) {
  BOOST_CHECK_EQUAL(findCharge(PdgParticle::eAntiElectron), 1_e);
  CHECK_CLOSE_REL(findMass(PdgParticle::eAntiElectron), 511_keV, eps);
  BOOST_CHECK_EQUAL(findName(PdgParticle::eAntiElectron), "e+");
  BOOST_CHECK_EQUAL(findCharge(PdgParticle::eElectron), -1_e);
  CHECK_CLOSE_REL(findMass(PdgParticle::eElectron), 511_keV, eps);
  BOOST_CHECK_EQUAL(findName(PdgParticle::eElectron), "e-");
  BOOST_CHECK_EQUAL(findCharge(PdgParticle::ePositron), 1_e);
  CHECK_CLOSE_REL(findMass(PdgParticle::ePositron), 511_keV, eps);
  BOOST_CHECK_EQUAL(findName(PdgParticle::ePositron), "e+");
}

BOOST_AUTO_TEST_CASE(Gamma) {
  BOOST_CHECK_EQUAL(findCharge(PdgParticle::eGamma), 0);
  BOOST_CHECK_EQUAL(findMass(PdgParticle::eGamma), 0);
  BOOST_CHECK_EQUAL(findName(PdgParticle::eGamma), "gamma");
}

BOOST_AUTO_TEST_CASE(Pion) {
  BOOST_CHECK_EQUAL(findCharge(PdgParticle::ePionMinus), -1_e);
  CHECK_CLOSE_REL(findMass(PdgParticle::ePionMinus), 139.57_MeV, eps);
  BOOST_CHECK_EQUAL(findName(PdgParticle::ePionMinus), "pi-");
  BOOST_CHECK_EQUAL(findCharge(PdgParticle::ePionPlus), 1_e);
  CHECK_CLOSE_REL(findMass(PdgParticle::ePionPlus), 139.57_MeV, eps);
  BOOST_CHECK_EQUAL(findName(PdgParticle::ePionPlus), "pi+");
  BOOST_CHECK_EQUAL(findCharge(PdgParticle::ePionZero), 0);
  CHECK_CLOSE_REL(findMass(PdgParticle::ePionZero), 134.98_MeV, eps);
  BOOST_CHECK_EQUAL(findName(PdgParticle::ePionZero), "pi0");
}

BOOST_AUTO_TEST_SUITE_END()
