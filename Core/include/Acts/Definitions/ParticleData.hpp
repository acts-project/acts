// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"

#include <iosfwd>
#include <optional>
#include <string_view>

namespace Acts {

struct ParticleData {
  float charge{};
  float mass{};
  std::string_view name;
};

/// Find the charge for a given PDG particle number.
/// @param pdg PDG particle number
/// @return Charge in native units.
std::optional<float> findCharge(PdgParticle pdg);

/// Find the charge for a given PDG particle number of a nucleus.
/// Try its ground state first, and ultimately get the proton number from PDG
/// @param pdg PDG particle number for the nucleus
/// @return Charge in native units.
float findChargeOfNucleus(PdgParticle pdg);

/// Find the mass for a given PDG particle number.
/// @param pdg PDG particle number
/// @return Mass in native units.
std::optional<float> findMass(PdgParticle pdg);

/// Find the mass for a given PDG particle number of a nucleus.
/// Try its ground state first, and ultimately get the mass from
/// Bethe-Weizsacker formula
/// @param pdg PDG particle number for the nucleus
/// @return Mass in native units
float findMassOfNucleus(PdgParticle pdg);

/// Calculate the mass of a nucleus using Bethe-Weizsacker formula
/// Parameters obtained from https://www.actaphys.uj.edu.pl/R/37/6/1833
///
/// @param pdg PDG particle number for the nucleus
/// @return Mass in native units
float calculateNucleusMass(PdgParticle pdg);

/// Find a descriptive particle name for a given PDG particle number.
/// @param pdg PDG particle number
/// @return Particle name.
std::optional<std::string_view> findName(PdgParticle pdg);

/// Find a descriptive particle name for a given PDG particle number of a
/// nucleus. Try to get the name from its ground state.
/// @param pdg PDG particle number for the nucleus
/// @return Particle name.
std::optional<std::string_view> findNameOfNucleus(PdgParticle pdg);

/// Find all known particle data for a given PDG particle number.
/// @param pdg PDG particle number
/// @return Particle data if found
std::optional<ParticleData> findParticleData(PdgParticle pdg);

/// Print PDG particle numbers with a descriptive name.
/// @param os Output stream
/// @param pdg PDG particle to output
/// @return Reference to output stream
std::ostream& operator<<(std::ostream& os, PdgParticle pdg);

/// Get short absolute string representation of PDG particle
/// @param pdg PDG particle number
/// @return Optional string view of particle name
std::optional<std::string_view> pdgToShortAbsString(PdgParticle pdg);

/// @brief Particle identification based on PDG number
namespace ParticleIdHelper {

/// @brief Check if the particle is a hadron
/// @param pdg PDG particle number
/// @return True if the particle is a hadron, false otherwise
bool isHadron(PdgParticle pdg);

/// @brief Check if the particle is a lepton
/// @param pdg PDG particle number
/// @return True if the particle is a lepton, false otherwise
bool isLepton(PdgParticle pdg);

/// @brief Check if the particle is a muon
/// @param pdg PDG particle number
/// @return True if the particle is a muon, false otherwise
bool isMuon(PdgParticle pdg);

/// @brief Check if the particle is an electron
/// @param pdg PDG particle number
/// @return True if the particle is an electron, false otherwise
bool isElectron(PdgParticle pdg);

/// @brief Check if the particle is a photon
/// @param pdg PDG particle number
/// @return True if the particle is a photon, false otherwise
bool isPhoton(PdgParticle pdg);

/// @brief Check if the particle is a tau
/// @param pdg PDG particle number
/// @return True if the particle is a tau, false otherwise
bool isTau(PdgParticle pdg);

/// @brief Check if the particle is a quark
/// @param pdg PDG particle number
/// @return True if the particle is a quark, false otherwise
bool isQuark(PdgParticle pdg);

/// @brief Check if the particle is interacting
/// @param pdg PDG particle number
/// @return True if the particle is interacting, false otherwise
bool isInteracting(PdgParticle pdg);

HadronType hadronType(PdgParticle pdg);

}  // namespace ParticleIdHelper

}  // namespace Acts
