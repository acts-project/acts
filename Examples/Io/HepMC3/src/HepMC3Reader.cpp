// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Reader.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Metadata.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Util.hpp"
#include "ActsExamples/Utilities/MultiplicityGenerators.hpp"

#include <filesystem>
#include <memory>

#include <HepMC3/GenEvent.h>
#include <HepMC3/Print.h>
#include <HepMC3/Reader.h>
#include <HepMC3/Units.h>
#include <boost/algorithm/string/join.hpp>

using namespace Acts::UnitLiterals;

namespace ActsExamples {

HepMC3Reader::HepMC3Reader(const HepMC3Reader::Config& cfg,
                           Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("HepMC3Reader", lvl)) {
  if (m_cfg.outputEvent.empty()) {
    throw std::invalid_argument("Missing output collection");
  }

  m_outputEvent.initialize(m_cfg.outputEvent);

  // Validate: exactly one of inputPath or inputs must be set
  bool hasInputPath = m_cfg.inputPath.has_value();
  bool hasInputs = !m_cfg.inputs.empty();

  if (!hasInputPath && !hasInputs) {
    throw std::invalid_argument(
        "HepMC3 reader requires either 'inputPath' or 'inputs' to be set");
  }

  if (hasInputPath && hasInputs) {
    throw std::invalid_argument(
        "HepMC3 reader: 'inputPath' and 'inputs' are mutually exclusive. "
        "Use 'inputPath' for single file or 'inputs' for multiple files.");
  }

  // If inputPath is set, create a single input with default multiplicity
  // generator
  if (hasInputPath) {
    Input input;
    input.path = m_cfg.inputPath.value();
    input.multiplicityGenerator =
        std::make_shared<FixedMultiplicityGenerator>(1);

    auto reader = HepMC3Util::deduceReader(input.path);
    m_inputs.emplace_back(reader, input.path, input.multiplicityGenerator);
  } else {
    // Use the provided inputs
    for (const auto& input : m_cfg.inputs) {
      if (!input.multiplicityGenerator) {
        throw std::invalid_argument(
            "All Input objects must have a multiplicityGenerator set");
      }
      auto reader = HepMC3Util::deduceReader(input.path);
      m_inputs.emplace_back(reader, input.path, input.multiplicityGenerator);
    }
  }

  if (m_cfg.numEvents.has_value()) {
    m_eventsRange = {0, m_cfg.numEvents.value()};
  } else {
    // Need to make a temporary reader to determine the number of events
    Acts::ScopedTimer timer("Determining number of events by reading", logger(),
                            Acts::Logging::DEBUG);
    // This uses the first reader that's configured, with the assumption that
    // this is the hard-scatter event
    m_eventsRange = {0, determineNumEvents(m_inputs.front().path)};
  }

  // Check if any input uses a non-Fixed multiplicity generator
  m_hasNonFixedMultiplicity =
      std::ranges::any_of(m_inputs, [](const auto& input) {
        // Check if this is NOT a FixedMultiplicityGenerator
        return dynamic_cast<const FixedMultiplicityGenerator*>(
                   input.multiplicityGenerator.get()) == nullptr;
      });

  // Check if randomNumbers is required
  // RNG is needed if we have a vertex generator or any non-Fixed multiplicity
  // generator
  m_needsRng = (m_cfg.vertexGenerator != nullptr) || m_hasNonFixedMultiplicity;

  if (m_needsRng && m_cfg.randomNumbers == nullptr) {
    throw std::invalid_argument(
        "randomNumbers must be set if vertexGenerator or any non-Fixed "
        "multiplicityGenerator is used");
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

std::size_t HepMC3Reader::determineNumEvents(
    const std::filesystem::path& path) const {
  std::optional metadata = HepMC3Metadata::readSidecar(path);

  if (metadata.has_value()) {
    std::size_t numEvents = metadata.value().numEvents;
    ACTS_INFO("HepMC3Reader: Found sidecar metadata file for "
              << path << " with " << numEvents << " events");
    return numEvents;
  }

  auto reader = HepMC3Util::deduceReader(m_inputs.front().path);
  ACTS_INFO(
      "HepMC3Reader: Number of events not specified, will read the "
      "whole file to determine the number of events. This might take a while "
      "and is usually not what you want.");
  std::size_t numEvents = 0;
  auto event =
      std::make_shared<HepMC3::GenEvent>(HepMC3::Units::GEV, HepMC3::Units::MM);
  while (!reader->failed()) {
    reader->read_event(*event);
    if (!reader->failed()) {
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

  ACTS_DEBUG("Skipping " << events << " logical events");

  // Check if all multiplicity generators are Fixed (deterministic)
  // Use the precomputed flag
  if (m_hasNonFixedMultiplicity) {
    ACTS_ERROR(
        "Cannot skip events with non-Fixed multiplicityGenerator (e.g., "
        "PoissonMultiplicityGenerator). Skipping requires knowing the exact "
        "number of physical events to skip from each input file, which is only "
        "possible with deterministic (Fixed) multiplicity generators.");
    return ABORT;
  }

  // For each logical event to skip, evaluate the Fixed multiplicity generators
  // to determine how many physical events to skip from each input file
  for (std::size_t logicalEvent = 0; logicalEvent < events; ++logicalEvent) {
    for (const auto& input : m_inputs) {
      // Must be FixedMultiplicityGenerator (checked above), so downcast is safe
      const auto* fixedGen = static_cast<const FixedMultiplicityGenerator*>(
          input.multiplicityGenerator.get());
      std::size_t count = fixedGen->n;

      ACTS_VERBOSE("Skipping " << count << " events from " << input.path
                               << " for logical event "
                               << (m_nextEvent + logicalEvent));
      if (!input.reader->skip(static_cast<int>(count))) {
        ACTS_ERROR("Error skipping " << count << " events from " << input.path);
        return ABORT;
      }
    }
  }

  m_nextEvent += events;

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
    const AlgorithmContext& ctx,
    std::shared_ptr<HepMC3::GenEvent>& outputEvent) {
  using enum ProcessCode;
  ACTS_VERBOSE("Reading from single file");

  std::vector<std::shared_ptr<HepMC3::GenEvent>> events;
  {
    std::scoped_lock lock(m_queueMutex);

    if (m_bufferError) {
      ACTS_ERROR("Buffer error (maybe in other thread), aborting");
      return ABORT;
    }

    // Check if we already read this event on another thread
    if (!m_events.empty() && ctx.eventNumber < m_nextEvent) {
      if (readCached(ctx, events) != SUCCESS) {
        ACTS_ERROR("Error reading event " << ctx.eventNumber);
        return ABORT;
      }
    } else {
      if (readBuffer(ctx, events) != SUCCESS) {
        ACTS_ERROR("Error reading event " << ctx.eventNumber);
        return ABORT;
      }
    }
  }

  if (m_cfg.vertexGenerator != nullptr) {
    Acts::ScopedTimer timer("Shifting events to vertex", logger(),
                            Acts::Logging::DEBUG);
    auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);
    for (auto& event : events) {
      auto vertexPosition = (*m_cfg.vertexGenerator)(rng);

      ACTS_VERBOSE("Shifting event to " << vertexPosition.transpose());
      // Our internal time unit is ctau, so is HepMC3's, make sure we convert
      // to mm
      HepMC3::FourVector vtxPosHepMC(vertexPosition[Acts::eFreePos0] / 1_mm,
                                     vertexPosition[Acts::eFreePos1] / 1_mm,
                                     vertexPosition[Acts::eFreePos2] / 1_mm,
                                     vertexPosition[Acts::eFreeTime] / 1_mm);
      event->shift_position_to(vtxPosHepMC);
    }
  }

  outputEvent = makeEvent();

  std::vector<const HepMC3::GenEvent*> eventPtrs;
  eventPtrs.reserve(events.size());
  std::ranges::transform(events, std::back_inserter(eventPtrs),
                         [](auto& event) { return event.get(); });
  {
    Acts::ScopedTimer timer("Merging events", logger(), Acts::Logging::DEBUG);
    HepMC3Util::mergeEvents(*outputEvent, eventPtrs, logger());
  }

  outputEvent->set_event_number(events.front()->event_number());

  return SUCCESS;
}

ProcessCode HepMC3Reader::readCached(
    const AlgorithmContext& ctx,
    std::vector<std::shared_ptr<HepMC3::GenEvent>>& events) {
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

  auto [num, genEvents] = std::move(*it);
  m_events.erase(it);
  events = std::move(genEvents);
  return ProcessCode::SUCCESS;
}

ProcessCode HepMC3Reader::readBuffer(
    const AlgorithmContext& ctx,
    std::vector<std::shared_ptr<HepMC3::GenEvent>>& outputEvents) {
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

    m_events.emplace_back(thisEvent, std::move(events));

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

  auto [num, genEvents] = std::move(m_events.back());
  m_events.pop_back();

  ACTS_VERBOSE("Popping event " << num << " from queue");

  outputEvents = std::move(genEvents);

  return SUCCESS;
}

ProcessCode HepMC3Reader::readLogicalEvent(
    const AlgorithmContext& ctx,
    std::vector<std::shared_ptr<HepMC3::GenEvent>>& events) {
  using enum ProcessCode;
  ACTS_VERBOSE("Reading logical event " << ctx.eventNumber);
  Acts::ScopedTimer timer("Reading logical event", logger(),
                          Acts::Logging::DEBUG);

  // @TODO: Add the index as an attribute to the event and it's content

  // Spawn RNG for multiplicity generators if needed
  std::optional<RandomEngine> rng;
  if (m_needsRng) {
    rng = m_cfg.randomNumbers->spawnGenerator(ctx);
  }

  for (std::size_t inputIndex = 0; inputIndex < m_inputs.size(); ++inputIndex) {
    auto& reader = m_inputs[inputIndex].reader;
    auto& path = m_inputs[inputIndex].path;
    auto& multiplicityGenerator = m_inputs[inputIndex].multiplicityGenerator;

    // Use multiplicityGenerator to determine count
    std::size_t count = 0;
    if (rng.has_value()) {
      count = (*multiplicityGenerator)(*rng);
    } else {
      // Must be FixedMultiplicityGenerator if no RNG, so downcast is safe
      const auto& fixedGen = static_cast<const FixedMultiplicityGenerator&>(
          *multiplicityGenerator);
      count = fixedGen.n;
    }

    ACTS_VERBOSE("Reading " << count << " events from " << path);
    for (std::size_t i = 0; i < count; ++i) {
      auto event = makeEvent();

      reader->read_event(*event);
      if (reader->failed()) {
        ACTS_ERROR("Error reading event " << i << " (input index = "
                                          << inputIndex << ") from " << path);
        if (inputIndex > 0) {
          ACTS_ERROR("-> since this is input file index "
                     << inputIndex
                     << ", this probably means that the "
                        "input file has "
                        "fewer events than expected.");
        }
        return ABORT;
      }
      events.push_back(std::move(event));
      m_inputs[inputIndex].eventsRead++;
    }
  }

  ACTS_VERBOSE("Read " << events.size() << " events in total from all files");

  return SUCCESS;
}

ProcessCode HepMC3Reader::finalize() {
  ACTS_VERBOSE("Closing " << m_inputs.size() << " input files");

  std::size_t totalEventsRead = 0;

  // Print summary of events read from each input file
  if (m_inputs.size() > 1) {
    ACTS_INFO("Events read per input file:");
    for (std::size_t i = 0; i < m_inputs.size(); ++i) {
      const auto& input = m_inputs[i];
      ACTS_INFO("  [" << i << "] " << input.path.filename() << ": "
                      << input.eventsRead << " events");
      totalEventsRead += input.eventsRead;
      input.reader->close();
    }
    ACTS_INFO("Total physical events read: " << totalEventsRead);
  } else {
    // Single input file
    for (const auto& input : m_inputs) {
      ACTS_INFO("Events read from " << input.path.filename() << ": "
                                    << input.eventsRead);
      input.reader->close();
    }
  }

  ACTS_DEBUG("Maximum event buffer size was: "
             << m_maxEventBufferSize << ", limit=" << m_cfg.maxEventBufferSize);
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
