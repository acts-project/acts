// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <HepMC3/FourVector.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepPID/ParticleIDMethods.hh>

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "Acts/Utilities/Units.hpp"

namespace FW {

/// Helper struct to convert HepMC3 event to the internal format.
struct HepMC3Event {
 public:
  ///
  /// Setter
  ///

  /// @brief Sets new units for momentums
  /// @note The allowed units are MeV and Gev
  /// @param event event in HepMC data type
  /// @param momentumUnit new unit of momentum
  void momentumUnit(std::shared_ptr<HepMC3::GenEvent> event,
                    const double momentumUnit);

  /// @brief Sets new units for lengths
  /// @note The allowed units are mm and cm
  /// @param event event in HepMC data type
  /// @param lengthUnit new unit of length
  void lengthUnit(std::shared_ptr<HepMC3::GenEvent> event,
                  const double lengthUnit);

  /// @brief Shifts the positioning of an event in space and time
  /// @param event event in HepMC data type
  /// @param deltaPos relative spatial shift that will be applied
  /// @param deltaTime relative time shift that will be applied
  void shiftPositionBy(std::shared_ptr<HepMC3::GenEvent> event,
                       const Acts::Vector3D& deltaPos, const double deltaTime);

  /// @brief Shifts the positioning of an event to a paint in space and time
  /// @param event event in HepMC data type
  /// @param pos new position of the event
  /// @param time new time of the event
  void shiftPositionTo(std::shared_ptr<HepMC3::GenEvent> event,
                       const Acts::Vector3D& pos, const double time);

  /// @brief Shifts the positioning of an event to a paint in space
  /// @param event event in HepMC data type
  /// @param pos new position of the event
  void shiftPositionTo(std::shared_ptr<HepMC3::GenEvent> event,
                       const Acts::Vector3D& pos);

  /// @brief Shifts the positioning of an event to a paint in time
  /// @param event event in HepMC data type
  /// @param time new time of the event
  void shiftPositionTo(std::shared_ptr<HepMC3::GenEvent> event,
                       const double time);

  ///
  /// Adder
  ///

  /// @brief Adds a new particle
  /// @param event event in HepMC data type
  /// @param particle new particle that will be added
  void addParticle(std::shared_ptr<HepMC3::GenEvent> event,
                   std::shared_ptr<SimParticle> particle);

  /// @brief Adds a new vertex
  /// @param event event in HepMC data type
  /// @param vertex new vertex that will be added
  /// @note The statuses are not represented in Acts and therefore set to 0
  void addVertex(std::shared_ptr<HepMC3::GenEvent> event,
                 const std::shared_ptr<SimVertex> vertex);
  ///
  /// Remover
  ///

  /// @brief Removes a particle from the record
  /// @param event event in HepMC data type
  /// @param particle particle that will be removed
  void removeParticle(std::shared_ptr<HepMC3::GenEvent> event,
                      const std::shared_ptr<SimParticle>& particle);

  /// @brief Removes a vertex from the record
  /// @note The identification of the vertex is potentially unstable (c.f.
  /// HepMC3Event::compareVertices())
  /// @param event event in HepMC data type
  /// @param vertex vertex that will be removed
  void removeVertex(std::shared_ptr<HepMC3::GenEvent> event,
                    const std::shared_ptr<SimVertex>& vertex);

  ///
  /// Getter
  ///

  /// @brief Getter of the unit of momentum used
  /// @param event event in HepMC data type
  /// @return unit of momentum
  double momentumUnit(const std::shared_ptr<HepMC3::GenEvent> event);

  /// @brief Getter of the unit of length used
  /// @param event event in HepMC data type
  /// @return unit of length
  double lengthUnit(const std::shared_ptr<HepMC3::GenEvent> event);

  /// @brief Getter of the position of the event
  /// @param event event in HepMC data type
  /// @return vector to the location of the event
  Acts::Vector3D eventPos(const std::shared_ptr<HepMC3::GenEvent> event);

  /// @brief Getter of the time of the event
  /// @param event event in HepMC data type
  /// @return time of the event
  double eventTime(const std::shared_ptr<HepMC3::GenEvent> event);

  /// @brief Get list of particles
  /// @param event event in HepMC data type
  /// @return List of particles
  std::vector<std::unique_ptr<SimParticle>> particles(
      const std::shared_ptr<HepMC3::GenEvent> event);

  /// @brief Get list of vertices
  /// @param event event in HepMC data type
  /// @return List of vertices
  std::vector<std::unique_ptr<SimVertex>> vertices(
      const std::shared_ptr<HepMC3::GenEvent> event);

  /// @brief Get beam particles
  /// @param event event in HepMC data type
  /// @return List of beam particles
  std::vector<std::unique_ptr<SimParticle>> beams(
      const std::shared_ptr<HepMC3::GenEvent> event);

  /// @brief Get final state particles
  /// @param event event in HepMC data type
  /// @return List of final state particles
  std::vector<std::unique_ptr<SimParticle>> finalState(
      const std::shared_ptr<HepMC3::GenEvent> event);

 private:
  /// @brief Converts an SimParticle into HepMC3::GenParticle
  /// @note The conversion ignores HepMC status codes
  /// @param actsParticle Acts particle that will be converted
  /// @return converted particle
  HepMC3::GenParticlePtr actsParticleToGen(
      std::shared_ptr<SimParticle> actsParticle);

  /// @brief Converts an Acts vertex to a HepMC3::GenVertexPtr
  /// @note The conversion ignores HepMC status codes
  /// @param actsVertex Acts vertex that will be converted
  /// @return Converted Acts vertex to HepMC3::GenVertexPtr
  HepMC3::GenVertexPtr createGenVertex(
      const std::shared_ptr<SimVertex>& actsVertex);

  /// @brief Compares an Acts vertex with a HepMC3::GenVertex
  /// @note An Acts vertex does not store a barcode. Therefore the content of
  /// both vertices is compared. The position, time and number of incoming and
  /// outgoing particles will be compared. Since a second vertex could exist in
  /// the record with identical informations (although unlikely), this
  /// comparison could lead to false positive results. On the other hand, a
  /// numerical deviation of the parameters could lead to a false negative.
  /// @param actsVertex Acts vertex
  /// @param genVertex HepMC3::GenVertex
  /// @return boolean result if both vertices are identical
  bool compareVertices(const std::shared_ptr<SimVertex>& actsVertex,
                       const HepMC3::GenVertexPtr& genVertex);
};
}  // namespace FW
