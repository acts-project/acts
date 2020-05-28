// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <fstream>

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Plugins/HepMC3/HepMC3Event.hpp"
#include "ACTFW/Plugins/HepMC3/HepMC3Reader.hpp"
#include "HepMC3/ReaderAscii.h"
#include "HepPID/ParticleName.hh"

///
/// Straight forward example of reading a HepMC3 file.
///
int main(int /*argc*/, char** /*argv*/) {
  FW::HepMC3ReaderAscii simReader;

  std::cout << "Preparing reader " << std::flush;
  HepMC3::ReaderAscii reader("test.hepmc3");
  if (simReader.status(reader)) {
    std::cout << "succesful" << std::endl;
  } else {
    std::cout << "failed" << std::endl;
  }

  std::shared_ptr<HepMC3::GenEvent> genevt(new HepMC3::GenEvent());

  std::cout << "Reading event " << std::flush;
  if (simReader.readEvent(reader, genevt)) {
    std::cout << "succesful" << std::endl;
  } else {
    std::cout << "failed" << std::endl;
  }
  std::cout << std::endl;

  FW::HepMC3Event simEvent;

  std::cout << "Event data:" << std::endl;
  std::cout << "Units: ";
  if (simEvent.momentumUnit(genevt) == Acts::units::_GeV)
    std::cout << "[GEV], ";
  else if (simEvent.momentumUnit(genevt) == Acts::units::_MeV)
    std::cout << "[MeV], ";
  if (simEvent.lengthUnit(genevt) == Acts::units::_mm)
    std::cout << "[mm]" << std::endl;
  else if (simEvent.lengthUnit(genevt) == Acts::units::_cm)
    std::cout << "[cm]" << std::endl;
  Acts::Vector3D evtPos = simEvent.eventPos(genevt);
  std::cout << "Event position: " << evtPos(0) << ", " << evtPos(1) << ", "
            << evtPos(2) << std::endl;
  std::cout << "Event time: " << simEvent.eventTime(genevt) << std::endl;

  std::cout << "Beam particles: ";
  std::vector<std::unique_ptr<FW::SimParticle>> beam = simEvent.beams(genevt);
  if (beam.empty())
    std::cout << "none" << std::endl;
  else {
    for (auto& pbeam : beam)
      std::cout << HepPID::particleName(pbeam->pdg()) << " ";
    std::cout << std::endl;
  }

  std::cout << std::endl << "Vertices: ";
  std::vector<std::unique_ptr<FW::SimVertex>> vertices =
      simEvent.vertices(genevt);
  if (vertices.empty())
    std::cout << "none" << std::endl;
  else {
    std::cout << std::endl;
    for (auto& vertex : vertices) {
      for (auto& particle : vertex->incoming)
        std::cout << HepPID::particleName(particle.pdg()) << " ";
      std::cout << "-> ";
      for (auto& particle : vertex->outgoing)
        std::cout << HepPID::particleName(particle.pdg()) << " ";
      std::cout << "\t@(" << vertex->time() << ", " << vertex->position()(0)
                << ", " << vertex->position()(1) << ", "
                << vertex->position()(2) << ")" << std::endl;
    }
    std::cout << std::endl;
  }

  std::cout << "Total particle record:" << std::endl;
  std::vector<std::unique_ptr<FW::SimParticle>> particles =
      simEvent.particles(genevt);
  for (auto& particle : particles)
    std::cout << HepPID::particleName(particle->pdg())
              << "\tID:" << particle->particleId() << ", momentum: ("
              << particle->momentum4()(0) << ", " << particle->momentum4()(1)
              << ", " << particle->momentum4()(2)
              << "), mass:  " << particle->mass() << std::endl;

  std::cout << std::endl << "Initial to final state: ";
  std::vector<std::unique_ptr<FW::SimParticle>> fState =
      simEvent.finalState(genevt);
  for (auto& pbeam : beam)
    std::cout << HepPID::particleName(pbeam->pdg()) << " ";
  std::cout << "-> ";
  for (auto& fs : fState)
    std::cout << HepPID::particleName(fs->pdg()) << " ";
  std::cout << std::endl;
}
