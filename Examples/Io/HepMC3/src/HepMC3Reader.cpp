// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Reader.hpp"

#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <HepMC3/Print.h>
#include <HepMC3/Units.h>

namespace ActsExamples {

bool HepMC3AsciiReader::readEvent(HepMC3::ReaderAscii& reader,
                                  HepMC3::GenEvent& event) {
  // Read event and store it
  return reader.read_event(event);
}

bool HepMC3AsciiReader::status(HepMC3::ReaderAscii& reader) {
  return !reader.failed();
}

HepMC3AsciiReader::HepMC3AsciiReader(const HepMC3AsciiReader::Config& cfg,
                                     Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_eventsRange(
          determineEventFilesRange(cfg.inputDir, cfg.inputStem + ".hepmc3")),
      m_logger(Acts::getDefaultLogger("HepMC3AsciiReader", lvl)) {
  if (m_cfg.inputStem.empty()) {
    throw std::invalid_argument("Missing input filename stem");
  }
  if (m_cfg.outputEvent.empty()) {
    throw std::invalid_argument("Missing output collection");
  }

  m_outputEvent.initialize(m_cfg.outputEvent);
}

std::string HepMC3AsciiReader::HepMC3AsciiReader::name() const {
  return "HepMC3AsciiReader";
}

std::pair<std::size_t, std::size_t> HepMC3AsciiReader::availableEvents() const {
  return m_eventsRange;
}

ProcessCode HepMC3AsciiReader::read(const AlgorithmContext& ctx) {
  auto event =
      std::make_shared<HepMC3::GenEvent>(HepMC3::Units::GEV, HepMC3::Units::MM);

  auto path = perEventFilepath(m_cfg.inputDir, m_cfg.inputStem + ".hepmc3",
                               ctx.eventNumber);

  ACTS_DEBUG("Attempting to read event from " << path);
  HepMC3::ReaderAscii reader(path);

  reader.read_event(*event);
  if (reader.failed()) {
    return ProcessCode::ABORT;
  }

  if (m_cfg.printListing) {
    ACTS_VERBOSE("Generated event:\n"
                 << [&]() {
                      std::stringstream ss;
                      HepMC3::Print::listing(ss, *event);
                      return ss.str();
                    }());
  }

  m_outputEvent(ctx, std::move(event));

  reader.close();
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
