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
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace Acts::ParticleId;  // Add this for easier access to particle ID
                                   // functions

namespace {
// NOTE: the used mass comparison values are not as exact as the data values
static constexpr float eps = 0.001f;
}  // namespace

BOOST_AUTO_TEST_SUITE(FatrasParticleData)

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

BOOST_AUTO_TEST_CASE(ParticleIdentification) {
  // Test leptons
  BOOST_CHECK(ParticleId::isLepton(PdgParticle::eElectron));
  BOOST_CHECK(ParticleId::isLepton(PdgParticle::ePositron));
  BOOST_CHECK(ParticleId::isLepton(PdgParticle::eMuon));
  BOOST_CHECK(ParticleId::isLepton(PdgParticle::eTau));
  BOOST_CHECK(!ParticleId::isLepton(PdgParticle::eProton));
  BOOST_CHECK(!ParticleId::isLepton(PdgParticle::ePionPlus));

  // Test hadrons
  BOOST_CHECK(ParticleId::isHadron(PdgParticle::eProton));
  BOOST_CHECK(ParticleId::isHadron(PdgParticle::eNeutron));
  BOOST_CHECK(ParticleId::isHadron(PdgParticle::ePionPlus));
  BOOST_CHECK(ParticleId::isHadron(PdgParticle::ePionZero));
  BOOST_CHECK(ParticleId::isHadron(PdgParticle::eKaonPlus));
  BOOST_CHECK(!ParticleId::isHadron(PdgParticle::eElectron));
  BOOST_CHECK(!ParticleId::isHadron(PdgParticle::eGamma));

  // Test interacting particles
  BOOST_CHECK(ParticleId::isInteracting(PdgParticle::eProton));
  BOOST_CHECK(ParticleId::isInteracting(PdgParticle::ePionPlus));
  BOOST_CHECK(ParticleId::isInteracting(PdgParticle::eGamma));
  BOOST_CHECK(ParticleId::isInteracting(PdgParticle::eElectron));
  // Neutrinos are weakly interacting
  constexpr int eNeutrinoE = 12;    // electron neutrino
  constexpr int eNeutrinoMu = 14;   // muon neutrino
  constexpr int eNeutrinoTau = 16;  // tau neutrino
  BOOST_CHECK(!ParticleId::isInteracting(eNeutrinoE));
  BOOST_CHECK(!ParticleId::isInteracting(eNeutrinoMu));
  BOOST_CHECK(!ParticleId::isInteracting(eNeutrinoTau));
  // Anti-neutrinos should also be non-interacting
  BOOST_CHECK(!ParticleId::isInteracting(-eNeutrinoE));
  BOOST_CHECK(!ParticleId::isInteracting(-eNeutrinoMu));
  BOOST_CHECK(!ParticleId::isInteracting(-eNeutrinoTau));

  // Test hadron types
  BOOST_CHECK_EQUAL(ParticleId::hadronType(PdgParticle::eProton),
                    HadronType::LightBaryon);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(PdgParticle::ePionPlus),
                    HadronType::LightMeson);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(PdgParticle::eKaonPlus),
                    HadronType::StrangeMeson);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(PdgParticle::eElectron),
                    HadronType::Unknown);
}

BOOST_AUTO_TEST_CASE(MesonIdentification) {
  // BBbar mesons (bottomonium)
  constexpr int Upsilon1S = 553;     // ϒ(1S)
  constexpr int Upsilon2S = 100553;  // ϒ(2S)
  constexpr int Upsilon3S = 200553;  // ϒ(3S)
  constexpr int ChiB = 10551;        // χb0
  BOOST_CHECK(ParticleId::isHadron(Upsilon1S));
  BOOST_CHECK(ParticleId::isHadron(Upsilon2S));
  BOOST_CHECK(ParticleId::isHadron(Upsilon3S));
  BOOST_CHECK(ParticleId::isHadron(ChiB));
  BOOST_CHECK_EQUAL(ParticleId::hadronType(Upsilon1S), HadronType::BBbarMeson);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(Upsilon2S), HadronType::BBbarMeson);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(Upsilon3S), HadronType::BBbarMeson);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(ChiB), HadronType::BBbarMeson);

  // CCbar mesons (charmonium)
  constexpr int JPsi = 443;      // J/ψ
  constexpr int Psi2S = 100443;  // ψ(2S)
  constexpr int ChiC = 10441;    // χc0
  constexpr int Eta_c = 441;     // ηc
  BOOST_CHECK(ParticleId::isHadron(JPsi));
  BOOST_CHECK(ParticleId::isHadron(Psi2S));
  BOOST_CHECK(ParticleId::isHadron(ChiC));
  BOOST_CHECK(ParticleId::isHadron(Eta_c));
  BOOST_CHECK_EQUAL(ParticleId::hadronType(JPsi), HadronType::CCbarMeson);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(Psi2S), HadronType::CCbarMeson);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(ChiC), HadronType::CCbarMeson);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(Eta_c), HadronType::CCbarMeson);

  // B mesons (bottom quark)
  constexpr int B0 = 511;     // B0 meson (bd)
  constexpr int BPlus = 521;  // B+ meson (bu)
  BOOST_CHECK(ParticleId::isHadron(B0));
  BOOST_CHECK(ParticleId::isHadron(BPlus));
  BOOST_CHECK_EQUAL(ParticleId::hadronType(B0), HadronType::BottomMeson);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(BPlus), HadronType::BottomMeson);

  // D mesons (charm quark)
  constexpr int D0 = 421;     // D0 meson (cu)
  constexpr int DPlus = 411;  // D+ meson (cd)
  BOOST_CHECK(ParticleId::isHadron(D0));
  BOOST_CHECK(ParticleId::isHadron(DPlus));
  BOOST_CHECK_EQUAL(ParticleId::hadronType(D0), HadronType::CharmedMeson);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(DPlus), HadronType::CharmedMeson);

  // Strange mesons
  BOOST_CHECK_EQUAL(ParticleId::hadronType(PdgParticle::eKaonPlus),
                    HadronType::StrangeMeson);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(PdgParticle::eKaonMinus),
                    HadronType::StrangeMeson);

  // Light mesons (pions)
  BOOST_CHECK_EQUAL(ParticleId::hadronType(PdgParticle::ePionPlus),
                    HadronType::LightMeson);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(PdgParticle::ePionMinus),
                    HadronType::LightMeson);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(PdgParticle::ePionZero),
                    HadronType::LightMeson);

  // Test anti-particles have same hadron type
  BOOST_CHECK_EQUAL(ParticleId::hadronType(B0), ParticleId::hadronType(-B0));
  BOOST_CHECK_EQUAL(ParticleId::hadronType(D0), ParticleId::hadronType(-D0));
}

BOOST_AUTO_TEST_CASE(HeavyBaryonIdentification) {
  // Bottom baryons
  constexpr int LambdaB = 5122;  // Λb0 (udb)
  constexpr int SigmaB = 5222;   // Σb+ (uub)
  constexpr int XiB = 5322;      // Ξb0 (usb)
  BOOST_CHECK(ParticleId::isHadron(LambdaB));
  BOOST_CHECK(ParticleId::isHadron(SigmaB));
  BOOST_CHECK(ParticleId::isHadron(XiB));
  BOOST_CHECK_EQUAL(ParticleId::hadronType(LambdaB), HadronType::BottomBaryon);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(SigmaB), HadronType::BottomBaryon);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(XiB), HadronType::BottomBaryon);

  // Charmed baryons
  constexpr int LambdaC = 4122;  // Λc+ (udc)
  constexpr int SigmaC = 4222;   // Σc++ (uuc)
  constexpr int XiC = 4322;      // Ξc0 (usc)
  BOOST_CHECK(ParticleId::isHadron(LambdaC));
  BOOST_CHECK(ParticleId::isHadron(SigmaC));
  BOOST_CHECK(ParticleId::isHadron(XiC));
  BOOST_CHECK_EQUAL(ParticleId::hadronType(LambdaC), HadronType::CharmedBaryon);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(SigmaC), HadronType::CharmedBaryon);
  BOOST_CHECK_EQUAL(ParticleId::hadronType(XiC), HadronType::CharmedBaryon);

  // Test anti-particles have same hadron type
  BOOST_CHECK_EQUAL(ParticleId::hadronType(LambdaB),
                    ParticleId::hadronType(-LambdaB));
  BOOST_CHECK_EQUAL(ParticleId::hadronType(LambdaC),
                    ParticleId::hadronType(-LambdaC));
}

BOOST_AUTO_TEST_SUITE_END()
