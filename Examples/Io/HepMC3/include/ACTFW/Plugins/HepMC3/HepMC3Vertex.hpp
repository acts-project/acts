// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <HepMC3/FourVector.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"

namespace FW {

/// Helper struct to convert HepMC3 vertex into the internal format.
struct HepMC3Vertex {
 public:
  /// @brief Returns a vertex translated into Acts
  /// @param vertex vertex in HepMC data type
  /// @return corresponding Acts vertex
  std::unique_ptr<SimVertex> processVertex(
      const std::shared_ptr<HepMC3::GenVertex> vertex);

  /// @brief Returns a boolean expression if a vertex is in an event translated
  /// into Acts
  /// @param vertex vertex in HepMC data type
  /// @return boolean expression if the vertex is in an event
  bool inEvent(const std::shared_ptr<HepMC3::GenVertex> vertex);

  /// @brief Returns a vertex translated into Acts
  /// @param vertex vertex in HepMC data type
  /// @return id of the vertex
  int id(const std::shared_ptr<HepMC3::GenVertex> vertex);

  /// @brief Returns the incoming particles of a vertex translated into Acts
  /// @param vertex vertex in HepMC data type
  /// @return incoming particles of the vertex
  std::vector<SimParticle> particlesIn(
      const std::shared_ptr<HepMC3::GenVertex> vertex);

  /// @brief Returns the outgoing particles of a vertex translated into Acts
  /// @param vertex vertex in HepMC data type
  /// @return outgoing particles of the vertex
  std::vector<SimParticle> particlesOut(
      const std::shared_ptr<HepMC3::GenVertex> vertex);

  /// @brief Returns the position of a vertex translated into Acts
  /// @param vertex vertex in HepMC data type
  /// @return position of the vertex
  Acts::Vector3D position(const std::shared_ptr<HepMC3::GenVertex> vertex);

  /// @brief Returns the time of a vertex translated into Acts
  /// @param vertex vertex in HepMC data type
  /// @return time of the vertex
  double time(const std::shared_ptr<HepMC3::GenVertex> vertex);

  /// @brief Adds an incoming particle to a vertex translated into Acts
  /// @param vertex vertex in HepMC data type
  /// @param particle incoming particle that will be added
  void addParticleIn(std::shared_ptr<HepMC3::GenVertex> vertex,
                     std::shared_ptr<SimParticle> particle);

  /// @brief Adds an outgoing particle to a vertex translated into Acts
  /// @param vertex vertex in HepMC data type
  /// @param particle outgoing particle that will be added
  void addParticleOut(std::shared_ptr<HepMC3::GenVertex> vertex,
                      std::shared_ptr<SimParticle> particle);

  /// @brief Removes an incoming particle from a vertex translated into Acts
  /// @param vertex vertex in HepMC data type
  /// @param particle incoming particle that will be removed
  void removeParticleIn(std::shared_ptr<HepMC3::GenVertex> vertex,
                        std::shared_ptr<SimParticle> particle);

  /// @brief Removes an outgoing particle from a vertex translated into Acts
  /// @param vertex vertex in HepMC data type
  /// @param particle outgoing particle that will be removed
  void removeParticleOut(std::shared_ptr<HepMC3::GenVertex> vertex,
                         std::shared_ptr<SimParticle> particle);

  /// @brief Sets the position of a vertex translated into Acts
  /// @param vertex vertex in HepMC data type
  /// @param pos new position of the vertex
  void position(const std::shared_ptr<HepMC3::GenVertex> vertex,
                Acts::Vector3D pos);

  /// @brief Sets the time of a vertex translated into Acts
  /// @param vertex vertex in HepMC data type
  /// @param time new time of the vertex
  void time(const std::shared_ptr<HepMC3::GenVertex> vertex, double time);

 private:
  /// @brief Converts HepMC3::GenParticle objects into Acts
  /// @param genParticles list of HepMC3::GenParticle objects
  /// @return converted list
  std::vector<SimParticle> genParticlesToActs(
      const std::vector<HepMC3::GenParticlePtr>& genParticles);

  /// @brief Converts an SimParticle into HepMC3::GenParticle
  /// @note The conversion ignores HepMC status codes
  /// @param actsParticle Acts particle that will be converted
  /// @return converted particle
  HepMC3::GenParticlePtr actsParticleToGen(
      std::shared_ptr<SimParticle> actsParticle);

  /// @brief Finds a HepMC3::GenParticle from a list that matches an
  /// SimParticle object
  /// @param genParticles list of HepMC particles
  /// @param actsParticle Acts particle
  /// @return HepMC particle that matched with the Acts particle or nullptr if
  /// no match was found
  HepMC3::GenParticlePtr matchParticles(
      const std::vector<HepMC3::GenParticlePtr>& genParticles,
      std::shared_ptr<SimParticle> actsParticle);
};
}  // namespace FW
