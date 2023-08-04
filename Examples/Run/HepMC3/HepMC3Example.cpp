// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/ParticleData.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Event.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Reader.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/HepMC3Options.hpp"

#include <fstream>

#include <HepMC3/GenEvent.h>
#include <HepMC3/ReaderAscii.h>

///
/// Straight forward example of reading a HepMC3 file.
///
int main(int argc, char** argv) {
  // Declare the supported program options.
  // Setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  ActsExamples::Options::addHepMC3ReaderOptions(desc);
  auto vm = ActsExamples::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  auto logLevel = ActsExamples::Options::readLogLevel(vm);

  // Create the reader
  auto hepMC3ReaderConfig = ActsExamples::Options::readHepMC3ReaderOptions(vm);
  ActsExamples::HepMC3AsciiReader simReader(hepMC3ReaderConfig, logLevel);

  std::cout << "Preparing reader " << std::flush;
  HepMC3::ReaderAscii reader("test.hepmc3");
  if (simReader.status(reader)) {
    std::cout << "successful" << std::endl;
  } else {
    std::cout << "failed" << std::endl;
  }

  HepMC3::GenEvent genevt;

  std::cout << "Reading event " << std::flush;
  if (simReader.readEvent(reader, genevt)) {
    std::cout << "successful" << std::endl;
  } else {
    std::cout << "failed" << std::endl;
  }
  std::cout << std::endl;

  using namespace ActsExamples::HepMC3Event;
  std::cout << "Event data:" << std::endl;
  std::cout << "Units: ";
  if (momentumUnit(genevt) == Acts::UnitConstants::GeV) {
    std::cout << "[GEV], ";
  } else if (momentumUnit(genevt) == Acts::UnitConstants::MeV) {
    std::cout << "[MeV], ";
  }
  if (lengthUnit(genevt) == Acts::UnitConstants::mm) {
    std::cout << "[mm]" << std::endl;
  } else if (lengthUnit(genevt) == Acts::UnitConstants::cm) {
    std::cout << "[cm]" << std::endl;
  }
  Acts::Vector3 evtPos = eventPos(genevt);
  std::cout << "Event position: " << evtPos(0) << ", " << evtPos(1) << ", "
            << evtPos(2) << std::endl;
  std::cout << "Event time: " << eventTime(genevt) << std::endl;

  std::cout << "Beam particles: ";
  std::vector<ActsExamples::SimParticle> beam = beams(genevt);
  if (beam.empty()) {
    std::cout << "none" << std::endl;
  } else {
    for (auto& pbeam : beam) {
      std::cout << Acts::findName(static_cast<Acts::PdgParticle>(pbeam.pdg()))
                       .value_or("invalid")
                << " ";
    }
    std::cout << std::endl;
  }

  std::cout << std::endl << "Vertices: ";
  std::vector<std::unique_ptr<ActsExamples::SimVertex>> vertices =
      ActsExamples::HepMC3Event::vertices(genevt);
  if (vertices.empty()) {
    std::cout << "none" << std::endl;
  } else {
    std::cout << std::endl;
    for (auto& vertex : vertices) {
      for (auto& particle : vertex->incoming) {
        std::cout << Acts::findName(
                         static_cast<Acts::PdgParticle>(particle.pdg()))
                         .value_or("invalid")
                  << " ";
      }
      std::cout << "-> ";
      for (auto& particle : vertex->outgoing) {
        std::cout << Acts::findName(
                         static_cast<Acts::PdgParticle>(particle.pdg()))
                         .value_or("invalid")
                  << " ";
      }
      std::cout << "\t@(" << vertex->time() << ", " << vertex->position()(0)
                << ", " << vertex->position()(1) << ", "
                << vertex->position()(2) << ")" << std::endl;
    }
    std::cout << std::endl;
  }

  std::cout << "Total particle record:" << std::endl;
  std::vector<ActsExamples::SimParticle> particles =
      ActsExamples::HepMC3Event::particles(genevt);
  for (auto& particle : particles) {
    std::cout << Acts::findName(static_cast<Acts::PdgParticle>(particle.pdg()))
                     .value_or("invalid")
              << "\tID:" << particle.particleId() << ", momentum: ("
              << particle.fourMomentum()(0) << ", "
              << particle.fourMomentum()(1) << ", "
              << particle.fourMomentum()(2) << "), mass:  " << particle.mass()
              << std::endl;
  }

  std::cout << std::endl << "Initial to final state: ";
  std::vector<ActsExamples::SimParticle> fState = finalState(genevt);
  for (auto& pbeam : beam) {
    std::cout << Acts::findName(static_cast<Acts::PdgParticle>(pbeam.pdg()))
                     .value_or("invalid")
              << " ";
  }
  std::cout << "-> ";
  for (auto& fs : fState) {
    std::cout << Acts::findName(static_cast<Acts::PdgParticle>(fs.pdg()))
                     .value_or("invalid")
              << " ";
  }
  std::cout << std::endl;
}
