// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Geant4/PDGtoG4Converter.hpp"

#include <cstdlib>
#include <utility>

#include <G4AntiBMesonZero.hh>
#include <G4AntiBsMesonZero.hh>
#include <G4AntiDMesonZero.hh>
#include <G4AntiKaonZero.hh>
#include <G4AntiLambda.hh>
#include <G4AntiLambdacPlus.hh>
#include <G4AntiNeutrinoE.hh>
#include <G4AntiNeutrinoMu.hh>
#include <G4AntiNeutrinoTau.hh>
#include <G4AntiNeutron.hh>
#include <G4AntiOmegaMinus.hh>
#include <G4AntiOmegacZero.hh>
#include <G4AntiProton.hh>
#include <G4AntiSigmaMinus.hh>
#include <G4AntiSigmaPlus.hh>
#include <G4AntiSigmaZero.hh>
#include <G4AntiSigmacPlus.hh>
#include <G4AntiSigmacPlusPlus.hh>
#include <G4AntiSigmacZero.hh>
#include <G4AntiXiMinus.hh>
#include <G4AntiXiZero.hh>
#include <G4AntiXicPlus.hh>
#include <G4AntiXicZero.hh>
#include <G4BMesonMinus.hh>
#include <G4BMesonPlus.hh>
#include <G4BMesonZero.hh>
#include <G4BsMesonZero.hh>
#include <G4DMesonMinus.hh>
#include <G4DMesonPlus.hh>
#include <G4DMesonZero.hh>
#include <G4DsMesonMinus.hh>
#include <G4DsMesonPlus.hh>
#include <G4Electron.hh>
#include <G4Eta.hh>
#include <G4EtaPrime.hh>
#include <G4Gamma.hh>
#include <G4JPsi.hh>
#include <G4KaonMinus.hh>
#include <G4KaonPlus.hh>
#include <G4KaonZero.hh>
#include <G4KaonZeroLong.hh>
#include <G4KaonZeroShort.hh>
#include <G4Lambda.hh>
#include <G4LambdacPlus.hh>
#include <G4MuonMinus.hh>
#include <G4MuonPlus.hh>
#include <G4NeutrinoE.hh>
#include <G4NeutrinoMu.hh>
#include <G4NeutrinoTau.hh>
#include <G4Neutron.hh>
#include <G4OmegaMinus.hh>
#include <G4OmegacZero.hh>
#include <G4ParticleDefinition.hh>
#include <G4PionMinus.hh>
#include <G4PionPlus.hh>
#include <G4PionZero.hh>
#include <G4Positron.hh>
#include <G4Proton.hh>
#include <G4SigmaMinus.hh>
#include <G4SigmaPlus.hh>
#include <G4SigmaZero.hh>
#include <G4SigmacPlus.hh>
#include <G4SigmacPlusPlus.hh>
#include <G4SigmacZero.hh>
#include <G4TauMinus.hh>
#include <G4TauPlus.hh>
#include <G4XiMinus.hh>
#include <G4XiZero.hh>
#include <G4XicPlus.hh>
#include <G4XicZero.hh>

ActsFatras::PDGtoG4Converter::PDGtoG4Converter() : m_pdgG4ParticleMap() {
  /// Fill the storage
  fillPredefinedParticles();
}

G4ParticleDefinition* ActsFatras::PDGtoG4Converter::getParticleDefinition(
    Acts::PdgParticle pdgCode) const {
  std::unordered_map<Acts::PdgParticle, G4ParticleDefinition*>::const_iterator
      it = m_pdgG4ParticleMap.find(pdgCode);

  // Find the corresponding particle and return it
  if (it != m_pdgG4ParticleMap.end()) {
    return it->second;
  } else {
    // Rescue mechanism: If the particle is neutral and its anti-particle is
    // stored then return that one instead
    it = m_pdgG4ParticleMap.find(makeAbsolutePdgParticle(pdgCode));
    if (it != m_pdgG4ParticleMap.end() &&
        std::abs(it->second->GetPDGCharge()) < 0.1) {
      return it->second;
    }
  }
  // No particle found
  return nullptr;
}

void ActsFatras::PDGtoG4Converter::fillPredefinedParticles() {
  // Gauge and Higgs Bosons
  addParticle(G4Gamma::GammaDefinition());

  // Leptons
  addParticle(G4Electron::ElectronDefinition());
  addParticle(G4NeutrinoE::NeutrinoEDefinition());
  addParticle(G4MuonMinus::MuonMinusDefinition());
  addParticle(G4NeutrinoMu::NeutrinoMuDefinition());
  addParticle(G4TauMinus::TauMinusDefinition());
  addParticle(G4NeutrinoTau::NeutrinoTauDefinition());
  addParticle(G4Positron::PositronDefinition());
  addParticle(G4AntiNeutrinoE::AntiNeutrinoEDefinition());
  addParticle(G4MuonPlus::MuonPlusDefinition());
  addParticle(G4AntiNeutrinoMu::AntiNeutrinoMuDefinition());
  addParticle(G4TauPlus::TauPlusDefinition());
  addParticle(G4AntiNeutrinoTau::AntiNeutrinoTauDefinition());

  // Light I=1 Mesons
  addParticle(G4PionZero::PionZeroDefinition());
  addParticle(G4PionPlus::PionPlusDefinition());
  addParticle(G4PionMinus::PionMinusDefinition());

  // Light I=0 Mesons
  addParticle(G4Eta::EtaDefinition());
  addParticle(G4EtaPrime::EtaPrimeDefinition());

  // Strange Mesons
  addParticle(G4KaonZeroLong::KaonZeroLongDefinition());
  addParticle(G4KaonZeroShort::KaonZeroShortDefinition());
  addParticle(G4KaonZero::KaonZeroDefinition());
  addParticle(G4KaonPlus::KaonPlusDefinition());
  addParticle(G4AntiKaonZero::AntiKaonZeroDefinition());
  addParticle(G4KaonMinus::KaonMinusDefinition());

  // Charmed Mesons
  addParticle(G4DMesonPlus::DMesonPlusDefinition());
  addParticle(G4DMesonZero::DMesonZeroDefinition());
  addParticle(G4DsMesonPlus::DsMesonPlusDefinition());
  addParticle(G4DMesonMinus::DMesonMinusDefinition());
  addParticle(G4AntiDMesonZero::AntiDMesonZeroDefinition());
  addParticle(G4DsMesonMinus::DsMesonMinusDefinition());

  // Bottom Mesons
  addParticle(G4BMesonZero::BMesonZeroDefinition());
  addParticle(G4BMesonPlus::BMesonPlusDefinition());
  addParticle(G4BsMesonZero::BsMesonZeroDefinition());
  addParticle(G4AntiBMesonZero::AntiBMesonZeroDefinition());
  addParticle(G4BMesonMinus::BMesonMinusDefinition());
  addParticle(G4AntiBsMesonZero::AntiBsMesonZeroDefinition());

  // ccbar Mesons
  addParticle(G4JPsi::JPsiDefinition());

  // Light Baryons
  addParticle(G4Proton::ProtonDefinition());
  addParticle(G4Neutron::NeutronDefinition());
  addParticle(G4AntiProton::AntiProtonDefinition());
  addParticle(G4AntiNeutron::AntiNeutronDefinition());

  // Strange Baryons
  addParticle(G4Lambda::LambdaDefinition());
  addParticle(G4SigmaPlus::SigmaPlusDefinition());
  addParticle(G4SigmaZero::SigmaZeroDefinition());
  addParticle(G4SigmaMinus::SigmaMinusDefinition());
  addParticle(G4XiZero::XiZeroDefinition());
  addParticle(G4XiMinus::XiMinusDefinition());
  addParticle(G4OmegaMinus::OmegaMinusDefinition());
  addParticle(G4AntiLambda::AntiLambdaDefinition());
  addParticle(G4AntiSigmaPlus::AntiSigmaPlusDefinition());
  addParticle(G4AntiSigmaZero::AntiSigmaZeroDefinition());
  addParticle(G4AntiSigmaMinus::AntiSigmaMinusDefinition());
  addParticle(G4AntiXiZero::AntiXiZeroDefinition());
  addParticle(G4AntiXiMinus::AntiXiMinusDefinition());
  addParticle(G4AntiOmegaMinus::AntiOmegaMinusDefinition());

  // Charmed Baryons
  addParticle(G4LambdacPlus::LambdacPlusDefinition());
  addParticle(G4SigmacPlusPlus::SigmacPlusPlusDefinition());
  addParticle(G4SigmacPlus::SigmacPlusDefinition());
  addParticle(G4SigmacZero::SigmacZeroDefinition());
  addParticle(G4XicPlus::XicPlusDefinition());
  addParticle(G4XicZero::XicZeroDefinition());
  addParticle(G4OmegacZero::OmegacZeroDefinition());
  addParticle(G4AntiLambdacPlus::AntiLambdacPlusDefinition());
  addParticle(G4AntiSigmacPlusPlus::AntiSigmacPlusPlusDefinition());
  addParticle(G4AntiSigmacPlus::AntiSigmacPlusDefinition());
  addParticle(G4AntiSigmacZero::AntiSigmacZeroDefinition());
  addParticle(G4AntiXicPlus::AntiXicPlusDefinition());
  addParticle(G4AntiXicZero::AntiXicZeroDefinition());
  addParticle(G4AntiOmegacZero::AntiOmegacZeroDefinition());
}

void ActsFatras::PDGtoG4Converter::addParticle(G4ParticleDefinition* pDef) {
  if (pDef == nullptr) {
    return;
  }

  m_pdgG4ParticleMap[static_cast<Acts::PdgParticle>(pDef->GetPDGEncoding())] =
      pDef;
}
