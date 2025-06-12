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
#include "ActsExamples/Io/HepMC3/HepMC3Util.hpp"
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

  if (!m_cfg.inputPaths.empty() && !m_cfg.inputPath.empty()) {
    throw std::invalid_argument(
        "inputPath and inputPaths are mutually exclusive");
  }

  if (!m_cfg.inputPath.empty()) {
    m_cfg.inputPaths.emplace_back(m_cfg.inputPath, 1);
  }

  if (m_cfg.inputPaths.empty()) {
    throw std::invalid_argument(
        "HepMC3 reader was not configured with any input files");
  }

  for (const auto& [path, numEvents] : m_cfg.inputPaths) {
    m_inputs.emplace_back(HepMC3::deduce_reader(path), numEvents, path);
  }

  if (m_cfg.numEvents.has_value()) {
    m_eventsRange = {0, m_cfg.numEvents.value()};
  } else {
    // Need to make a temporary reader to determine the number of events
    Acts::ScopedTimer timer("Determining number of events by reading", logger(),
                            Acts::Logging::DEBUG);
    // This uses the first reader that's configured, with the assumption that
    // this is the hard-scatter event
    auto reader = HepMC3::deduce_reader(m_inputs.front().path);
    m_eventsRange = {0, determineNumEvents(*reader)};
  }

  ACTS_DEBUG("HepMC3Reader: " << m_eventsRange.first << " - "
                              << m_eventsRange.second << " events");
}

HepMC3Reader::~HepMC3Reader() = default;

std::string HepMC3Reader::name() const {
  return "HepMC3Reader";
}

std::pair<std::size_t, std::size_t> HepMC3Reader::availableEvents() const {
  return m_eventsRange;
}

std::size_t HepMC3Reader::determineNumEvents(HepMC3::Reader& reader) const {
  ACTS_WARNING(
      "HepMC3Reader: Number of events not specified, will read the "
      "whole file to determine the number of events. This might take a while "
      "and is usually not what you want. Set the numEvents parameter to the "
      "number of events you want to read.");
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

  ACTS_DEBUG("Skipping " << events << " events");

  for (auto& [reader, numEvents, path] : m_inputs) {
    ACTS_VERBOSE("Skipping " << events << "*" << numEvents << "="
                             << events * numEvents << " events from " << path);
    if (!reader->skip(static_cast<int>(events * numEvents))) {
      ACTS_ERROR("Error skipping events " << events << " " << path);
      return ABORT;
    }
  }

  m_nextEvent = events;

  return SUCCESS;
}

ProcessCode HepMC3Reader::read(const AlgorithmContext& ctx) {
  std::shared_ptr<HepMC3::GenEvent> event;

  using enum ProcessCode;

  if (readSingleFile(ctx, event) != SUCCESS) {
    return ABORT;
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

ProcessCode HepMC3Reader::readSingleFile(
    const ActsExamples::AlgorithmContext& ctx,
    std::shared_ptr<HepMC3::GenEvent>& outputEvent) {
  using enum ProcessCode;
  ACTS_VERBOSE("Reading from single file");
  std::scoped_lock lock(m_mutex);

  if (m_bufferError) {
    ACTS_ERROR("Buffer error (maybe in other thread), aborting");
    return ABORT;
  }

  // Check if we already read this event on another thread
  if (!m_events.empty() && ctx.eventNumber < m_nextEvent) {
    if (readCached(ctx, outputEvent) != SUCCESS) {
      ACTS_ERROR("Error reading event " << ctx.eventNumber);
      return ABORT;
    }
  } else {
    if (readBuffer(ctx, outputEvent) != SUCCESS) {
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
          "Event "
          << ctx.eventNumber
          << " was not found in the queue. Most likely the manually "
             "configured event count of the input file(s) of "
          << m_cfg.numEvents.value()
          << " is larger than the number of events found in the file(s)");
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

ProcessCode HepMC3Reader::readBuffer(
    const ActsExamples::AlgorithmContext& ctx,
    std::shared_ptr<HepMC3::GenEvent>& outputEvent) {
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

    std::vector<std::shared_ptr<HepMC3::GenEvent>> events;

    if (readLogicalEvent(ctx, events) != SUCCESS) {
      ACTS_ERROR("Error reading event " << thisEvent);
      return ABORT;
    }

    if (events.empty()) {
      ACTS_ERROR("No events read from file, this is a bug");
      return ABORT;
    }

    auto genEvent = makeEvent();

    std::vector<const HepMC3::GenEvent*> eventPtrs;
    eventPtrs.reserve(events.size());
    std::ranges::transform(events, std::back_inserter(eventPtrs),
                           [](auto& event) { return event.get(); });
    HepMC3Util::mergeEvents(*genEvent, eventPtrs, logger());

    genEvent->set_event_number(events.front()->event_number());

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

  outputEvent = std::move(genEvent);

  return SUCCESS;
}

ProcessCode HepMC3Reader::readLogicalEvent(
    const ActsExamples::AlgorithmContext& ctx,
    std::vector<std::shared_ptr<HepMC3::GenEvent>>& events) {
  using enum ProcessCode;
  ACTS_VERBOSE("Reading logical event " << ctx.eventNumber);

  // @TODO: Add the index as an attribute to the event and it's content

  for (auto& [reader, numEvents, path] : m_inputs) {
    ACTS_VERBOSE("Reading " << numEvents << " events from " << path);
    for (std::size_t i = 0; i < numEvents; ++i) {
      auto event = makeEvent();

      reader->read_event(*event);
      if (reader->failed()) {
        ACTS_ERROR("Error reading event " << i << " from " << path);
        return ABORT;
      }
      events.push_back(std::move(event));
    }
  }

  ACTS_VERBOSE("Read " << events.size() << " events in total from all files");

  return SUCCESS;
}

ProcessCode HepMC3Reader::finalize() {
  ACTS_VERBOSE("Closing " << m_inputs.size() << " input files");
  for (auto& [reader, numEvents, path] : m_inputs) {
    reader->close();
  }

  ACTS_DEBUG("Maximum event buffer size was: "
             << m_maxEventBufferSize << ", limit=" << m_cfg.maxEventBufferSize);
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
