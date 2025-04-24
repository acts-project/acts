// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Reader.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <memory>

#include <HepMC3/GenEvent.h>
#include <HepMC3/Print.h>
#include <HepMC3/Reader.h>
#include <HepMC3/ReaderFactory.h>
#include <HepMC3/Units.h>
#include <boost/algorithm/string/join.hpp>

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
      Acts::ScopedTimer timer("Determining number of events by reading",
                              logger(), Acts::Logging::DEBUG);
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

std::shared_ptr<HepMC3::GenEvent> HepMC3Reader::makeEvent() {
  return std::make_shared<HepMC3::GenEvent>(HepMC3::Units::GEV,
                                            HepMC3::Units::MM);
}

ProcessCode HepMC3Reader::skip(std::size_t events) {
  using enum ProcessCode;
  if (events == 0) {
    return SUCCESS;
  }

  if (m_cfg.perEvent) {
    // nothing to do as we lookup the target file by filename
    return SUCCESS;
  }

  ACTS_DEBUG("Skipping " << events << " events");
  if (!m_reader->skip(static_cast<int>(events))) {
    ACTS_ERROR("Error skipping events " << events << " " << m_cfg.inputPath);
    return ABORT;
  }
  m_nextEvent = events;

  return SUCCESS;
}

ProcessCode HepMC3Reader::read(const AlgorithmContext& ctx) {
  std::shared_ptr<HepMC3::GenEvent> event;

  using enum ProcessCode;

  if (m_cfg.perEvent) {
    if (readPerEvent(ctx, event) != SUCCESS) {
      return ABORT;
    }
  } else {
    if (readSingleFile(ctx, event) != SUCCESS) {
      return ABORT;
    }
  }

  throw_assert(event != nullptr, "Event should not be null");

  if (m_cfg.printListing) {
    ACTS_VERBOSE("Read event:\n"
                 << [&]() {
                      std::stringstream ss;
                      HepMC3::Print::listing(ss, *event);
                      return ss.str();
                    }());
  }

  if (m_cfg.checkEventNumber &&
      static_cast<std::size_t>(event->event_number()) != ctx.eventNumber) {
    ACTS_ERROR("HepMC3Reader: Event number mismatch. Expected "
               << ctx.eventNumber << ", but got " << event->event_number()
               << ". You can turn this off in the configuration if your events "
                  "were not written in order.");
    return ABORT;
  }

  m_outputEvent(ctx, std::move(event));

  return SUCCESS;
}

ProcessCode HepMC3Reader::readPerEvent(
    const ActsExamples::AlgorithmContext& ctx,
    std::shared_ptr<HepMC3::GenEvent>& event) {
  using enum ProcessCode;
  std::filesystem::path perEventFile =
      perEventFilepath(m_cfg.inputPath.parent_path(),
                       m_cfg.inputPath.filename().string(), ctx.eventNumber);
  ACTS_VERBOSE("Reading per-event file " << perEventFile);
  HepMC3::ReaderAscii reader(perEventFile);
  event = makeEvent();
  reader.read_event(*event);
  if (reader.failed()) {
    return ABORT;
  }
  reader.close();

  return SUCCESS;
}

ProcessCode HepMC3Reader::readSingleFile(
    const ActsExamples::AlgorithmContext& ctx,
    std::shared_ptr<HepMC3::GenEvent>& event) {
  using enum ProcessCode;
  ACTS_VERBOSE("Reading from single file");
  std::scoped_lock lock(m_mutex);

  if (m_bufferError) {
    ACTS_ERROR("Buffer error (maybe in other thread), aborting");
    return ABORT;
  }

  // Check if we already read this event on another thread
  if (!m_events.empty() && ctx.eventNumber < m_nextEvent) {
    if (readCached(ctx, event) != SUCCESS) {
      ACTS_ERROR("Error reading event " << ctx.eventNumber);
      return ABORT;
    }
  } else {
    if (readBuffer(ctx, event) != SUCCESS) {
      ACTS_ERROR("Error reading event " << ctx.eventNumber);
      return ABORT;
    }
  }

  return SUCCESS;
}

ProcessCode HepMC3Reader::readCached(const ActsExamples::AlgorithmContext& ctx,
                                     std::shared_ptr<HepMC3::GenEvent>& event) {
  ACTS_VERBOSE("Already read event " << ctx.eventNumber);
  auto it = std::ranges::find_if(
      m_events, [&ctx](auto& elem) { return elem.first == ctx.eventNumber; });

  if (it == m_events.end()) {
    if (m_cfg.numEvents.has_value()) {
      ACTS_ERROR(
          "Event " << ctx.eventNumber
                   << " was not found in the queue. Most likely the manually "
                      "configured event count of the input file of "
                   << m_cfg.numEvents.value()
                   << " is larger than the number of events found in the file");
    } else {
      ACTS_ERROR(
          "Event " << ctx.eventNumber
                   << "  should be in the queue, but is not. This is a bug.");
    }
    m_bufferError = true;
    return ProcessCode::ABORT;
  }

  auto [num, genEvent] = std::move(*it);
  m_events.erase(it);
  event = std::move(genEvent);
  return ProcessCode::SUCCESS;
}

ProcessCode HepMC3Reader::readBuffer(const ActsExamples::AlgorithmContext& ctx,
                                     std::shared_ptr<HepMC3::GenEvent>& event) {
  using enum ProcessCode;

  ACTS_VERBOSE("Next event to read is: " << m_nextEvent);

  std::size_t eventsToRead = ctx.eventNumber - m_nextEvent + 1;
  ACTS_VERBOSE("event_number=" << ctx.eventNumber
                               << " next_event=" << m_nextEvent
                               << " events_to_read=" << eventsToRead);

  if (m_events.size() + eventsToRead > m_cfg.maxEventBufferSize) {
    ACTS_ERROR("Event buffer size of "
               << m_cfg.maxEventBufferSize
               << " would be exceeded. This can happen in case there are many "
                  "threads and processing happens strongly out-of-order");
    m_bufferError = true;
    return ABORT;
  }

  do {
    ACTS_VERBOSE("Reading next event as number " << m_nextEvent);
    std::size_t thisEvent = m_nextEvent;
    m_nextEvent += 1;

    auto genEvent = makeEvent();
    m_reader->read_event(*genEvent);
    if (m_reader->failed()) {
      ACTS_ERROR("HepMC3Reader: Error reading event " << thisEvent << " "
                                                      << m_cfg.inputPath);
      m_bufferError = true;
      return ABORT;
    }

    m_events.emplace_back(thisEvent, std::move(genEvent));
  } while (m_nextEvent <= ctx.eventNumber);

  ACTS_VERBOSE("Queue is now: [" << [&]() {
    std::vector<std::string> numbers;
    std::ranges::transform(
        m_events, std::back_inserter(numbers),
        [](const auto& v) { return std::to_string(v.first); });
    return boost::algorithm::join(numbers, ", ");
  }() << "]");

  throw_assert(!m_events.empty(), "Event should not be empty, this is a bug.");

  m_maxEventBufferSize = std::max(m_maxEventBufferSize, m_events.size());

  auto [num, genEvent] = std::move(m_events.back());
  m_events.pop_back();

  ACTS_VERBOSE("Popping event " << num << " from queue");

  event = std::move(genEvent);

  return SUCCESS;
}

ProcessCode HepMC3Reader::finalize() {
  if (m_reader) {
    m_reader->close();
  }

  ACTS_DEBUG("Maximum event buffer size was: "
             << m_maxEventBufferSize << ", limit=" << m_cfg.maxEventBufferSize);
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
