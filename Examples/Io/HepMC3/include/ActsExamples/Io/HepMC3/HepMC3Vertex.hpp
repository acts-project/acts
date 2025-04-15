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

namespace ActsExamples::HepMC3Vertex {

/// @brief Returns a vertex translated into Acts
/// @param vertex vertex in HepMC data type
/// @return corresponding Acts vertex
std::unique_ptr<SimVertex> processVertex(const HepMC3::GenVertex& vertex);

/// @brief Returns a boolean expression if a vertex is in an event translated
/// into Acts
/// @param vertex vertex in HepMC data type
/// @return boolean expression if the vertex is in an event
bool inEvent(const HepMC3::GenVertex& vertex);

/// @brief Returns a vertex translated into Acts
/// @param vertex vertex in HepMC data type
/// @return id of the vertex
int id(const HepMC3::GenVertex& vertex);

/// @brief Returns the incoming particles of a vertex translated into Acts
/// @param vertex vertex in HepMC data type
/// @return incoming particles of the vertex
std::vector<SimParticle> particlesIn(const HepMC3::GenVertex& vertex);

/// @brief Returns the outgoing particles of a vertex translated into Acts
/// @param vertex vertex in HepMC data type
/// @return outgoing particles of the vertex
std::vector<SimParticle> particlesOut(const HepMC3::GenVertex& vertex);

/// @brief Returns the position of a vertex translated into Acts
/// @param vertex vertex in HepMC data type
/// @return position of the vertex
Acts::Vector3 position(const HepMC3::GenVertex& vertex);

/// @brief Returns the time of a vertex translated into Acts
/// @param vertex vertex in HepMC data type
/// @return time of the vertex
double time(const HepMC3::GenVertex& vertex);

/// @brief Adds an incoming particle to a vertex translated into Acts
/// @param vertex vertex in HepMC data type
/// @param particle incoming particle that will be added
void addParticleIn(HepMC3::GenVertex& vertex, const SimParticle& particle);

/// @brief Adds an outgoing particle to a vertex translated into Acts
/// @param vertex vertex in HepMC data type
/// @param particle outgoing particle that will be added
void addParticleOut(HepMC3::GenVertex& vertex, const SimParticle& particle);

/// @brief Removes an incoming particle from a vertex translated into Acts
/// @param vertex vertex in HepMC data type
/// @param particle incoming particle that will be removed
void removeParticleIn(HepMC3::GenVertex& vertex, const SimParticle& particle);

/// @brief Removes an outgoing particle from a vertex translated into Acts
/// @param vertex vertex in HepMC data type
/// @param particle outgoing particle that will be removed
void removeParticleOut(HepMC3::GenVertex& vertex, const SimParticle& particle);

/// @brief Sets the position of a vertex translated into Acts
/// @param vertex vertex in HepMC data type
/// @param pos new position of the vertex
void position(HepMC3::GenVertex& vertex, const Acts::Vector3& pos);

/// @brief Sets the time of a vertex translated into Acts
/// @param vertex vertex in HepMC data type
/// @param time new time of the vertex
void time(HepMC3::GenVertex& vertex, double time);

}  // namespace ActsExamples::HepMC3Vertex
