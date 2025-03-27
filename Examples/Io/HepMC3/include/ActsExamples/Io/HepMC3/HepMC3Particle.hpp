// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"

#include <HepMC3/FourVector.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>

namespace ActsExamples::HepMC3Particle {

/// @brief Returns the barcode translated into Acts
/// @param particle HepMC3::GenParticle particle
/// @return corresponding Acts barcode
SimBarcode barcode(const HepMC3::GenParticle& particle);

/// @brief Returns the particle translated into Acts
/// @param particle HepMC3::GenParticle particle
/// @return corresponding Acts particle
SimParticle particle(const HepMC3::GenParticle& particle);

/// @brief Returns the id of the particle translated into Acts
/// @param particle HepMC3::GenParticle particle
/// @return id of the particle
int id(const HepMC3::GenParticle& particle);

/// @brief Returns the production vertex of the particle translated into Acts
/// @param particle HepMC3::GenParticle particle
/// @return production vertex of the particle
std::unique_ptr<SimVertex> productionVertex(
    const HepMC3::GenParticle& particle);

/// @brief Returns the end vertex of the particle translated into Acts
/// @param particle HepMC3::GenParticle particle
/// @return end vertex of the particle
std::unique_ptr<SimVertex> endVertex(const HepMC3::GenParticle& particle);

/// @brief Returns the PDG code of a particle translated into Acts
/// @param particle HepMC3::GenParticle particle
/// @return PDG code of the particle
int pdgID(const HepMC3::GenParticle& particle);

/// @brief Returns the momentum of a particle translated into Acts
/// @param particle HepMC3::GenParticle particle
/// @return momentum of the particle
Acts::Vector3 momentum(const HepMC3::GenParticle& particle);

/// @brief Returns the energy of a particle translated into Acts
/// @param particle HepMC3::GenParticle particle
/// @return energy of the particle
double energy(const HepMC3::GenParticle& particle);

/// @brief Returns the mass of a particle translated into Acts
/// @param particle HepMC3::GenParticle particle
/// @return mass of the particle
double mass(const HepMC3::GenParticle& particle);

/// @brief Returns the charge of a particle translated into Acts
/// @param particle HepMC3::GenParticle particle
/// @return charge of the particle
double charge(const HepMC3::GenParticle& particle);

/// @brief Sets the PDG code of a particle translated from Acts
/// @param particle HepMC3::GenParticle particle
/// @param pid PDG code that will be set
void pdgID(HepMC3::GenParticle& particle, const int pid);

/// @brief Sets the momentum of a particle translated from Acts
/// @param particle HepMC3::GenParticle particle
/// @param mom momentum that will be set
void momentum(HepMC3::GenParticle& particle, const Acts::Vector3& mom);

/// @brief Sets the energy of a particle translated from Acts
/// @param particle HepMC3::GenParticle particle
/// @param energy energy that will be set
void energy(HepMC3::GenParticle& particle, const double energy);

/// @brief Sets the mass of a particle translated from Acts
/// @param particle HepMC3::GenParticle particle
/// @param mass mass that will be set
void mass(HepMC3::GenParticle& particle, const double mass);

}  // namespace ActsExamples::HepMC3Particle
