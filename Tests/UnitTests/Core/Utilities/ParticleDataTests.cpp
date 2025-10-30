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
  constexpr int eNeutrinoE = 12;    // electron neutrino
  constexpr int eNeutrinoMu = 14;   // muon neutrino
  constexpr int eNeutrinoTau = 16;  // tau neutrino
  BOOST_CHECK(!ParticleIdHelper::isInteracting(eNeutrinoE));
  BOOST_CHECK(!ParticleIdHelper::isInteracting(eNeutrinoMu));
  BOOST_CHECK(!ParticleIdHelper::isInteracting(eNeutrinoTau));
  // Anti-neutrinos should also be non-interacting
  BOOST_CHECK(!ParticleIdHelper::isInteracting(-eNeutrinoE));
  BOOST_CHECK(!ParticleIdHelper::isInteracting(-eNeutrinoMu));
  BOOST_CHECK(!ParticleIdHelper::isInteracting(-eNeutrinoTau));

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
  constexpr int Upsilon1S = 553;     // ϒ(1S)
  constexpr int Upsilon2S = 100553;  // ϒ(2S)
  constexpr int Upsilon3S = 200553;  // ϒ(3S)
  constexpr int ChiB = 10551;        // χb0
  BOOST_CHECK(ParticleIdHelper::isHadron(Upsilon1S));
  BOOST_CHECK(ParticleIdHelper::isHadron(Upsilon2S));
  BOOST_CHECK(ParticleIdHelper::isHadron(Upsilon3S));
  BOOST_CHECK(ParticleIdHelper::isHadron(ChiB));
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(Upsilon1S), HadronType::BBbarMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(Upsilon2S), HadronType::BBbarMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(Upsilon3S), HadronType::BBbarMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(ChiB), HadronType::BBbarMeson);

  // CCbar mesons (charmonium)
  constexpr int JPsi = 443;      // J/ψ
  constexpr int Psi2S = 100443;  // ψ(2S)
  constexpr int ChiC = 10441;    // χc0
  constexpr int Eta_c = 441;     // ηc
  BOOST_CHECK(ParticleIdHelper::isHadron(JPsi));
  BOOST_CHECK(ParticleIdHelper::isHadron(Psi2S));
  BOOST_CHECK(ParticleIdHelper::isHadron(ChiC));
  BOOST_CHECK(ParticleIdHelper::isHadron(Eta_c));
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(JPsi), HadronType::CCbarMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(Psi2S), HadronType::CCbarMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(ChiC), HadronType::CCbarMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(Eta_c), HadronType::CCbarMeson);

  // B mesons (bottom quark)
  constexpr int B0 = 511;     // B0 meson (bd)
  constexpr int BPlus = 521;  // B+ meson (bu)
  BOOST_CHECK(ParticleIdHelper::isHadron(B0));
  BOOST_CHECK(ParticleIdHelper::isHadron(BPlus));
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(B0), HadronType::BottomMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(BPlus), HadronType::BottomMeson);

  // D mesons (charm quark)
  constexpr int D0 = 421;     // D0 meson (cu)
  constexpr int DPlus = 411;  // D+ meson (cd)
  BOOST_CHECK(ParticleIdHelper::isHadron(D0));
  BOOST_CHECK(ParticleIdHelper::isHadron(DPlus));
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(D0), HadronType::CharmedMeson);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(DPlus), HadronType::CharmedMeson);

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
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(B0), ParticleIdHelper::hadronType(-B0));
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(D0), ParticleIdHelper::hadronType(-D0));
}

BOOST_AUTO_TEST_CASE(HeavyBaryonIdentification) {
  // Bottom baryons
  constexpr int LambdaB = 5122;  // Λb0 (udb)
  constexpr int SigmaB = 5222;   // Σb+ (uub)
  constexpr int XiB = 5322;      // Ξb0 (usb)
  BOOST_CHECK(ParticleIdHelper::isHadron(LambdaB));
  BOOST_CHECK(ParticleIdHelper::isHadron(SigmaB));
  BOOST_CHECK(ParticleIdHelper::isHadron(XiB));
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(LambdaB), HadronType::BottomBaryon);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(SigmaB), HadronType::BottomBaryon);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(XiB), HadronType::BottomBaryon);

  // Charmed baryons
  constexpr int LambdaC = 4122;  // Λc+ (udc)
  constexpr int SigmaC = 4222;   // Σc++ (uuc)
  constexpr int XiC = 4322;      // Ξc0 (usc)
  BOOST_CHECK(ParticleIdHelper::isHadron(LambdaC));
  BOOST_CHECK(ParticleIdHelper::isHadron(SigmaC));
  BOOST_CHECK(ParticleIdHelper::isHadron(XiC));
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(LambdaC), HadronType::CharmedBaryon);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(SigmaC), HadronType::CharmedBaryon);
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(XiC), HadronType::CharmedBaryon);

  // Test anti-particles have same hadron type
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(LambdaB),
                    ParticleIdHelper::hadronType(-LambdaB));
  BOOST_CHECK_EQUAL(ParticleIdHelper::hadronType(LambdaC),
                    ParticleIdHelper::hadronType(-LambdaC));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
