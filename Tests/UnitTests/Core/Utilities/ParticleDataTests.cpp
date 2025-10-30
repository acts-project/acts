// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/ParticleData.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace Acts::ParticleIdHelper;

namespace {
// NOTE: the used mass comparison values are not as exact as the data values
static constexpr float eps = 0.001f;
}  // namespace

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_AUTO_TEST_CASE(InvalidInput) {
  BOOST_CHECK(!findCharge(PdgParticle::eInvalid));
  BOOST_CHECK(!findMass(PdgParticle::eInvalid));
  BOOST_CHECK(!findName(PdgParticle::eInvalid));
}

BOOST_AUTO_TEST_CASE(Electron) {
  BOOST_CHECK_EQUAL(*findCharge(PdgParticle::eAntiElectron), 1_e);
  CHECK_CLOSE_REL(*findMass(PdgParticle::eAntiElectron), 511_keV, eps);
  BOOST_CHECK_EQUAL(*findName(PdgParticle::eAntiElectron), "e+");
  BOOST_CHECK_EQUAL(*findCharge(PdgParticle::eElectron), -1_e);
  CHECK_CLOSE_REL(*findMass(PdgParticle::eElectron), 511_keV, eps);
  BOOST_CHECK_EQUAL(*findName(PdgParticle::eElectron), "e-");
  BOOST_CHECK_EQUAL(*findCharge(PdgParticle::ePositron), 1_e);
  CHECK_CLOSE_REL(*findMass(PdgParticle::ePositron), 511_keV, eps);
  BOOST_CHECK_EQUAL(*findName(PdgParticle::ePositron), "e+");
}

BOOST_AUTO_TEST_CASE(Gamma) {
  BOOST_CHECK_EQUAL(*findCharge(PdgParticle::eGamma), 0);
  BOOST_CHECK_EQUAL(*findMass(PdgParticle::eGamma), 0);
  BOOST_CHECK_EQUAL(*findName(PdgParticle::eGamma), "gamma");
}

BOOST_AUTO_TEST_CASE(Pion) {
  BOOST_CHECK_EQUAL(*findCharge(PdgParticle::ePionMinus), -1_e);
  CHECK_CLOSE_REL(*findMass(PdgParticle::ePionMinus), 139.57_MeV, eps);
  BOOST_CHECK_EQUAL(*findName(PdgParticle::ePionMinus), "pi-");
  BOOST_CHECK_EQUAL(*findCharge(PdgParticle::ePionPlus), 1_e);
  CHECK_CLOSE_REL(*findMass(PdgParticle::ePionPlus), 139.57_MeV, eps);
  BOOST_CHECK_EQUAL(*findName(PdgParticle::ePionPlus), "pi+");
  BOOST_CHECK_EQUAL(*findCharge(PdgParticle::ePionZero), 0);
  CHECK_CLOSE_REL(*findMass(PdgParticle::ePionZero), 134.98_MeV, eps);
  BOOST_CHECK_EQUAL(*findName(PdgParticle::ePionZero), "pi0");
}

BOOST_AUTO_TEST_CASE(ParticleIdHelperentification) {
  // Test leptons
  BOOST_CHECK(ParticleIdHelper::isLepton(PdgParticle::eElectron));
  BOOST_CHECK(ParticleIdHelper::isLepton(PdgParticle::ePositron));
  BOOST_CHECK(ParticleIdHelper::isLepton(PdgParticle::eMuon));
  BOOST_CHECK(ParticleIdHelper::isLepton(PdgParticle::eTau));
  BOOST_CHECK(!ParticleIdHelper::isLepton(PdgParticle::eProton));
  BOOST_CHECK(!ParticleIdHelper::isLepton(PdgParticle::ePionPlus));

  // Test hadrons
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eProton));
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eNeutron));
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::ePionPlus));
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::ePionZero));
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eKaonPlus));
  BOOST_CHECK(!ParticleIdHelper::isHadron(PdgParticle::eElectron));
  BOOST_CHECK(!ParticleIdHelper::isHadron(PdgParticle::eGamma));

  // Test interacting particles
  BOOST_CHECK(ParticleIdHelper::isInteracting(PdgParticle::eProton));
  BOOST_CHECK(ParticleIdHelper::isInteracting(PdgParticle::ePionPlus));
  BOOST_CHECK(ParticleIdHelper::isInteracting(PdgParticle::eGamma));
  BOOST_CHECK(ParticleIdHelper::isInteracting(PdgParticle::eElectron));
  // Neutrinos are weakly interacting
  BOOST_CHECK(!ParticleIdHelper::isInteracting(PdgParticle::eNeutrinoE));
  BOOST_CHECK(!ParticleIdHelper::isInteracting(PdgParticle::eNeutrinoMu));
  BOOST_CHECK(!ParticleIdHelper::isInteracting(PdgParticle::eNeutrinoTau));
  // Anti-neutrinos should also be non-interacting
  BOOST_CHECK(!ParticleIdHelper::isInteracting(PdgParticle::eAntiNeutrinoE));
  BOOST_CHECK(!ParticleIdHelper::isInteracting(PdgParticle::eAntiNeutrinoMu));
  BOOST_CHECK(!ParticleIdHelper::isInteracting(PdgParticle::eAntiNeutrinoTau));

  // Test hadron types
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eProton),
                    HadronType::LightBaryon);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::ePionPlus),
                    HadronType::LightMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eKaonPlus),
                    HadronType::StrangeMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eElectron),
                    HadronType::Unknown);
}

BOOST_AUTO_TEST_CASE(MesonIdentification) {
  // BBbar mesons (bottomonium)
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eUpsilon1S));
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eUpsilon2S));
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eUpsilon3S));
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eChiB));
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eUpsilon1S),
                    HadronType::BBbarMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eUpsilon2S),
                    HadronType::BBbarMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eUpsilon3S),
                    HadronType::BBbarMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eChiB),
                    HadronType::BBbarMeson);

  // CCbar mesons (charmonium)
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eJPsi));
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::ePsi2S));
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eChiC));
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eEta_c));
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eJPsi),
                    HadronType::CCbarMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::ePsi2S),
                    HadronType::CCbarMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eChiC),
                    HadronType::CCbarMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eEta_c),
                    HadronType::CCbarMeson);

  // B mesons (bottom quark)
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eB0));
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eBPlus));
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eB0),
                    HadronType::BottomMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eBPlus),
                    HadronType::BottomMeson);

  // D mesons (charm quark)
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eD0));
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eDPlus));
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eD0),
                    HadronType::CharmedMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eDPlus),
                    HadronType::CharmedMeson);

  // Strange mesons
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eKaonPlus),
                    HadronType::StrangeMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eKaonMinus),
                    HadronType::StrangeMeson);

  // Light mesons (pions)
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::ePionPlus),
                    HadronType::LightMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::ePionMinus),
                    HadronType::LightMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::ePionZero),
                    HadronType::LightMeson);

  // Test anti-particles have same hadron type
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eB0),
                    ParticleIdHelper::hadronType(PdgParticle::eAntiB0));
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eD0),
                    ParticleIdHelper::hadronType(PdgParticle::eAntiD0));
}

BOOST_AUTO_TEST_CASE(HeavyBaryonIdentification) {
  // Bottom baryons
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eLambdaB));
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eSigmaB));
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eXiB));
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eLambdaB),
                    HadronType::BottomBaryon);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eSigmaB),
                    HadronType::BottomBaryon);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eXiB),
                    HadronType::BottomBaryon);

  // Charmed baryons
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eLambdaC));
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eSigmaC));
  BOOST_CHECK(ParticleIdHelper::isHadron(PdgParticle::eXiC));
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eLambdaC),
                    HadronType::CharmedBaryon);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eSigmaC),
                    HadronType::CharmedBaryon);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eXiC),
                    HadronType::CharmedBaryon);

  // Test anti-particles have same hadron type
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eLambdaB),
                    ParticleIdHelper::hadronType(PdgParticle::eAntiLambdaB));
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(PdgParticle::eLambdaC),
                    ParticleIdHelper::hadronType(PdgParticle::eAntiLambdaC));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
