// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Reader.hpp"

#include "ActsExamples/Utilities/Paths.hpp"

#include <HepMC3/GenEvent.h>
#include <HepMC3/Print.h>
#include <HepMC3/Reader.h>
#include <HepMC3/ReaderFactory.h>
#include <HepMC3/Units.h>

namespace ActsExamples {

HepMC3Reader::HepMC3Reader(const HepMC3Reader::Config& cfg,
                           Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("HepMC3Reader", lvl)) {
  if (m_cfg.outputEvent.empty()) {
    throw std::invalid_argument("Missing output collection");
  }

  m_outputEvent.initialize(m_cfg.outputEvent);

  if (!m_cfg.perEvent) {
    // Create a single file reader
    m_reader = makeReader();
    if (m_cfg.numEvents.has_value()) {
      m_eventsRange = {0, m_cfg.numEvents.value()};
    } else {
      // Need to make a temporary reader to determine the number of events
      m_eventsRange = {0, determineNumEvents(*makeReader())};
    }
  } else {
    if (m_cfg.numEvents.has_value()) {
      throw std::invalid_argument(
          "perEvent and numEvents are mutually exclusive");
    }
    m_eventsRange = determineEventFilesRange(
        m_cfg.inputPath.parent_path(), m_cfg.inputPath.filename().string());
  }

  ACTS_DEBUG("HepMC3Reader: " << m_eventsRange.first << " - "
                              << m_eventsRange.second << " events");
}

std::shared_ptr<HepMC3::Reader> HepMC3Reader::makeReader() const {
  return HepMC3::deduce_reader(m_cfg.inputPath);
}

HepMC3Reader::~HepMC3Reader() = default;

std::string HepMC3Reader::name() const {
  return "HepMC3Reader";
}

std::pair<std::size_t, std::size_t> HepMC3Reader::availableEvents() const {
  return m_eventsRange;
}

std::size_t HepMC3Reader::determineNumEvents(HepMC3::Reader& reader) const {
  ACTS_INFO(
      "HepMC3Reader: Number of events not specified, will read the "
      "whole file to determine the number of events. This might take a while.");
  std::size_t numEvents = 0;
  auto event =
      std::make_shared<HepMC3::GenEvent>(HepMC3::Units::GEV, HepMC3::Units::MM);
  while (!reader.failed()) {
    reader.read_event(*event);
    if (!reader.failed()) {
      ++numEvents;
    }
  }
  return numEvents;
}

ProcessCode HepMC3Reader::read(const AlgorithmContext& ctx) {
  ACTS_VERBOSE("Reading event " << ctx.eventNumber << " from "
                                << m_cfg.inputPath);

  auto event =
      std::make_shared<HepMC3::GenEvent>(HepMC3::Units::GEV, HepMC3::Units::MM);

  auto read = [&](HepMC3::Reader& reader) {
    reader.read_event(*event);
    if (reader.failed()) {
      return ProcessCode::ABORT;
    }
    return ProcessCode::SUCCESS;
  };

  if (m_cfg.perEvent) {
    std::filesystem::path perEventFile =
        perEventFilepath(m_cfg.inputPath.parent_path(),
                         m_cfg.inputPath.filename().string(), ctx.eventNumber);
    ACTS_VERBOSE("Reading per-event file " << perEventFile);
    HepMC3::ReaderAscii reader(perEventFile);
    auto result = read(reader);
    reader.close();
    if (result == ProcessCode::ABORT) {
      return ProcessCode::ABORT;
    }
  } else {
    ACTS_VERBOSE("Reading from single file " << m_cfg.inputPath);
    // Take the lock until the end of the function
    std::scoped_lock lock(m_mutex);
    auto result = read(*m_reader);
    if (result == ProcessCode::ABORT) {
      ACTS_ERROR("HepMC3Reader: Error reading event " << ctx.eventNumber << " "
                                                      << m_cfg.inputPath);
      if (m_cfg.numEvents.has_value()) {
        ACTS_ERROR("HepMC3Reader: Expected " << m_cfg.numEvents.value()
                                             << " events, but got "
                                             << ctx.eventNumber << " events");
      }
      return ProcessCode::ABORT;
    }
  }

  if (m_cfg.printListing) {
    ACTS_VERBOSE("Read event:\n"
                 << [&]() {
                      std::stringstream ss;
                      HepMC3::Print::listing(ss, *event);
                      return ss.str();
                    }());
  }

  m_outputEvent(ctx, std::move(event));

  return ProcessCode::SUCCESS;
}

ProcessCode HepMC3Reader::finalize() {
  if (m_reader) {
    m_reader->close();
  }
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
