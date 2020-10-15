// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Reader.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include <HepMC3/Units.h>

bool ActsExamples::HepMC3ReaderAscii::readEvent(
    HepMC3::ReaderAscii& reader, std::shared_ptr<HepMC3::GenEvent> event) {
  // Read event and store it
  return reader.read_event(*event);
}

bool ActsExamples::HepMC3ReaderAscii::status(HepMC3::ReaderAscii& reader) {
  return !reader.failed();
}

ActsExamples::HepMC3ReaderAscii::HepMC3ReaderAscii(
    const ActsExamples::HepMC3ReaderAscii::Config& cfg,
    Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_eventsRange(determineEventFilesRange(
          cfg.inputDir, cfg.inputStemFileName + ".hepmc3")),
      m_logger(Acts::getDefaultLogger("HepMC3ReaderAscii", lvl)) {
  if (m_cfg.inputStemFileName.empty()) {
    throw std::invalid_argument("Missing input filename stem");
  }
  if (m_cfg.outputCollection.empty()) {
    throw std::invalid_argument("Missing output collection");
  }
}

std::string ActsExamples::HepMC3ReaderAscii::HepMC3ReaderAscii::name() const {
  return "HepMC3ReaderAscii";
}

std::pair<size_t, size_t> ActsExamples::HepMC3ReaderAscii::availableEvents()
    const {
  return m_eventsRange;
}

ActsExamples::ProcessCode ActsExamples::HepMC3ReaderAscii::read(
    const ActsExamples::AlgorithmContext& ctx) {
  std::vector<std::shared_ptr<HepMC3::GenEvent>> events;
  HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);

  auto path = perEventFilepath(m_cfg.inputDir, m_cfg.inputStemFileName + ".csv",
                               ctx.eventNumber);
  HepMC3::ReaderAscii reader(path);

  while (reader.read_event(event)) {
    events.push_back(std::make_shared<HepMC3::GenEvent>(event));
    event.clear();
  }

  ctx.eventStore.add(m_cfg.outputCollection, std::move(events));

  return reader.failed() ? ActsExamples::ProcessCode::ABORT
                         : ActsExamples::ProcessCode::SUCCESS;
}