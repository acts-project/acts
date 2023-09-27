// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Framework/Sequencer.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/SequenceElement.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <algorithm>
#include <atomic>
#include <cfenv>
#include <chrono>
#include <exception>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <typeinfo>

#ifndef ACTS_EXAMPLES_NO_TBB
#include <TROOT.h>
#endif

#include <boost/algorithm/string.hpp>
#include <boost/core/demangle.hpp>
#include <dfe/dfe_io_dsv.hpp>
#include <dfe/dfe_namedtuple.hpp>

namespace ActsExamples {

namespace {

std::string_view getAlgorithmType(const SequenceElement& element) {
  if (dynamic_cast<const IWriter*>(&element) != nullptr) {
    return "Writer";
  }
  if (dynamic_cast<const IReader*>(&element) != nullptr) {
    return "Reader";
  }
  return "Algorithm";
}

// Saturated addition that does not overflow and exceed SIZE_MAX.
//
// From http://locklessinc.com/articles/sat_arithmetic/
size_t saturatedAdd(size_t a, size_t b) {
  size_t res = a + b;
  res |= -static_cast<int>(res < a);
  return res;
}

}  // namespace

Sequencer::Sequencer(const Sequencer::Config& cfg)
    : m_cfg(cfg),
      m_taskArena((m_cfg.numThreads < 0) ? tbb::task_arena::automatic
                                         : m_cfg.numThreads),
      m_logger(Acts::getDefaultLogger("Sequencer", m_cfg.logLevel)) {
#ifndef ACTS_EXAMPLES_NO_TBB
  if (m_cfg.numThreads == 1) {
#endif
    ACTS_INFO("Create Sequencer (single-threaded)");
#ifndef ACTS_EXAMPLES_NO_TBB
  } else {
    ROOT::EnableThreadSafety();
    ACTS_INFO("Create Sequencer with " << m_cfg.numThreads << " threads");
  }
#endif
}

void Sequencer::addContextDecorator(
    std::shared_ptr<IContextDecorator> decorator) {
  if (not decorator) {
    throw std::invalid_argument("Can not add empty/NULL context decorator");
  }
  m_decorators.push_back(std::move(decorator));
  ACTS_INFO("Added context decarator '" << m_decorators.back()->name() << "'");
}

void Sequencer::addReader(std::shared_ptr<IReader> reader) {
  if (not reader) {
    throw std::invalid_argument("Can not add empty/NULL reader");
  }
  m_readers.push_back(reader);
  addElement(std::move(reader));
}

void Sequencer::addAlgorithm(std::shared_ptr<IAlgorithm> algorithm) {
  if (not algorithm) {
    throw std::invalid_argument("Can not add empty/NULL algorithm");
  }

  addElement(std::move(algorithm));
}

void Sequencer::addWriter(std::shared_ptr<IWriter> writer) {
  if (not writer) {
    throw std::invalid_argument("Can not add empty/NULL writer");
  }
  addElement(std::move(writer));
}

void Sequencer::addElement(const std::shared_ptr<SequenceElement>& element) {
  if (not element) {
    throw std::invalid_argument("Can not add empty/NULL element");
  }

  m_sequenceElements.push_back(element);

  std::string elementType{getAlgorithmType(*element)};
  std::string elementTypeCapitalized = elementType;
  elementTypeCapitalized[0] = std::toupper(elementTypeCapitalized[0]);
  ACTS_INFO("Add " << elementType << " '" << element->name() << "'");

  if (!m_cfg.runDataFlowChecks) {
    return;
  }

  auto symbol = [&](const char* in) {
    std::string s = boost::core::demangle(in);
    size_t pos = 0;
    while (pos + 80 < s.size()) {
      ACTS_INFO("   " + s.substr(pos, pos + 80));
      pos += 80;
    }
    ACTS_INFO("   " + s.substr(pos));
  };

  bool valid = true;

  for (const auto* handle : element->readHandles()) {
    if (!handle->isInitialized()) {
      continue;
    }

    ACTS_INFO("<- " << handle->name() << " '" << handle->key() << "':");
    symbol(handle->typeInfo().name());

    if (auto it = m_whiteBoardState.find(handle->key());
        it != m_whiteBoardState.end()) {
      const auto& source = *it->second;
      if (!source.isCompatible(*handle)) {
        ACTS_ERROR("Adding "
                   << elementType << " " << element->name() << ":"
                   << "\n-> white board will contain key '" << handle->key()
                   << "'"
                   << "\nat this point in the sequence (source: "
                   << source.fullName() << "),"
                   << "\nbut the type will be\n"
                   << "'" << boost::core::demangle(source.typeInfo().name())
                   << "'"
                   << "\nand not\n"
                   << "'" << boost::core::demangle(handle->typeInfo().name())
                   << "'");
        valid = false;
      }
    } else {
      ACTS_ERROR("Adding " << elementType << " " << element->name() << ":"
                           << "\n-> white board will not contain key"
                           << " '" << handle->key()
                           << "' at this point in the sequence."
                           << "\n   Needed for read data handle '"
                           << handle->name() << "'")
      valid = false;
    }
  }

  if (valid) {  // only record outputs this if we're valid until here
    for (const auto* handle : element->writeHandles()) {
      if (!handle->isInitialized()) {
        continue;
      }

      ACTS_INFO("-> " << handle->name() << " '" << handle->key() << "':");
      symbol(handle->typeInfo().name());

      if (auto it = m_whiteBoardState.find(handle->key());
          it != m_whiteBoardState.end()) {
        const auto& source = *it->second;
        ACTS_ERROR("White board will already contain key '"
                   << handle->key() << "'. Source: '" << source.fullName()
                   << "' (cannot overwrite)");
        valid = false;
        break;
      }

      m_whiteBoardState.emplace(std::pair{handle->key(), handle});

      if (auto it = m_whiteboardObjectAliases.find(handle->key());
          it != m_whiteboardObjectAliases.end()) {
        ACTS_DEBUG("Key '" << handle->key() << "' aliased to '" << it->second
                           << "'");
        m_whiteBoardState[it->second] = handle;
      }
    }
  }

  if (!valid) {
    throw SequenceConfigurationException{};
  }
}

void Sequencer::addWhiteboardAlias(const std::string& aliasName,
                                   const std::string& objectName) {
  auto [it, success] =
      m_whiteboardObjectAliases.insert({objectName, aliasName});
  if (!success) {
    throw std::invalid_argument("Alias to '" + aliasName + "' -> '" +
                                objectName + "' already set");
  }

  if (auto oit = m_whiteBoardState.find(objectName);
      oit != m_whiteBoardState.end()) {
    m_whiteBoardState[aliasName] = oit->second;
  }
}

std::vector<std::string> Sequencer::listAlgorithmNames() const {
  std::vector<std::string> names;

  // WARNING this must be done in the same order as in the processing
  for (const auto& decorator : m_decorators) {
    names.push_back("Decorator:" + decorator->name());
  }
  for (const auto& algorithm : m_sequenceElements) {
    names.push_back(std::string(getAlgorithmType(*algorithm)) + ":" +
                    algorithm->name());
  }

  return names;
}

std::pair<std::size_t, std::size_t> Sequencer::determineEventsRange() const {
  constexpr auto kInvalidEventsRange = std::make_pair(SIZE_MAX, SIZE_MAX);

  // Note on skipping events:
  //
  // Previously, skipping events was only allowed when readers where
  // available, since only readers had a `.skip()` functionality. The
  // `.skip()` interface has been removed in favour of telling the readers the
  // event they are requested to read via the algorithm context. Skipping can
  // now also be used when no readers are configured, e.g. for generating only
  // a few specific events in a simulation setup.

  // determine intersection of event ranges available from readers
  size_t beg = 0u;
  size_t end = SIZE_MAX;
  for (const auto& reader : m_readers) {
    auto available = reader->availableEvents();
    beg = std::max(beg, available.first);
    end = std::min(end, available.second);
  }

  // since we use event ranges (and not just num events) they might not
  // overlap
  if (end < beg) {
    ACTS_ERROR("Available events ranges from readers do not overlap");
    return kInvalidEventsRange;
  }
  // configured readers without available events makes no sense
  // TODO could there be a use-case for zero events? run only setup functions?
  if (beg == end) {
    ACTS_ERROR("No events available");
    return kInvalidEventsRange;
  }
  // trying to skip too many events must be an error
  if (end <= saturatedAdd(beg, m_cfg.skip)) {
    ACTS_ERROR("Less events available than requested to skip");
    return kInvalidEventsRange;
  }
  // events range was not defined by either the readers or user command line.
  if ((beg == 0u) and (end == SIZE_MAX) and (!m_cfg.events.has_value())) {
    ACTS_ERROR("Could not determine number of events");
    return kInvalidEventsRange;
  }

  // take user selection into account
  auto begSelected = saturatedAdd(beg, m_cfg.skip);
  auto endSelected = end;
  if (m_cfg.events.has_value()) {
    auto endRequested = saturatedAdd(begSelected, m_cfg.events.value());
    endSelected = std::min(end, endRequested);
    if (end < endRequested) {
      ACTS_INFO("Restrict requested number of events to available ones");
    }
  }

  return {begSelected, endSelected};
}

// helpers for per-algorithm timing information
namespace {
using Clock = std::chrono::high_resolution_clock;
using Duration = Clock::duration;
using Timepoint = Clock::time_point;
using Seconds = std::chrono::duration<double>;
using NanoSeconds = std::chrono::duration<double, std::nano>;

// RAII-based stopwatch to time execution within a block
struct StopWatch {
  Timepoint start;
  Duration& store;

  StopWatch(Duration& s) : start(Clock::now()), store(s) {}
  ~StopWatch() { store += Clock::now() - start; }
};

// Convert duration to a printable string w/ reasonable unit.
template <typename D>
inline std::string asString(D duration) {
  double ns = std::chrono::duration_cast<NanoSeconds>(duration).count();
  if (1e9 < std::abs(ns)) {
    return std::to_string(ns / 1e9) + " s";
  } else if (1e6 < std::abs(ns)) {
    return std::to_string(ns / 1e6) + " ms";
  } else if (1e3 < std::abs(ns)) {
    return std::to_string(ns / 1e3) + " us";
  } else {
    return std::to_string(ns) + " ns";
  }
}

// Convert duration scaled to one event to a printable string.
template <typename D>
inline std::string perEvent(D duration, size_t numEvents) {
  return asString(duration / numEvents) + "/event";
}

// Store timing data
struct TimingInfo {
  std::string identifier;
  double time_total_s = 0;
  double time_perevent_s = 0;

  DFE_NAMEDTUPLE(TimingInfo, identifier, time_total_s, time_perevent_s);
};

void storeTiming(const std::vector<std::string>& identifiers,
                 const std::vector<Duration>& durations, std::size_t numEvents,
                 const std::string& path) {
  dfe::NamedTupleTsvWriter<TimingInfo> writer(path, 4);
  for (size_t i = 0; i < identifiers.size(); ++i) {
    TimingInfo info;
    info.identifier = identifiers[i];
    info.time_total_s =
        std::chrono::duration_cast<Seconds>(durations[i]).count();
    info.time_perevent_s = info.time_total_s / numEvents;
    writer.append(info);
  }
}
}  // namespace

int Sequencer::run() {
  // measure overall wall clock
  Timepoint clockWallStart = Clock::now();
  // per-algorithm time measures
  std::vector<std::string> names = listAlgorithmNames();
  std::vector<Duration> clocksAlgorithms(names.size(), Duration::zero());
  tbbWrap::queuing_mutex clocksAlgorithmsMutex;

  // processing only works w/ a well-known number of events
  // error message is already handled by the helper function
  std::pair<size_t, size_t> eventsRange = determineEventsRange();
  if ((eventsRange.first == SIZE_MAX) and (eventsRange.second == SIZE_MAX)) {
    return EXIT_FAILURE;
  }

  ACTS_INFO("Processing events [" << eventsRange.first << ", "
                                  << eventsRange.second << ")");
  ACTS_INFO("Starting event loop with " << m_cfg.numThreads << " threads");
  ACTS_INFO("  " << m_decorators.size() << " context decorators");
  ACTS_INFO("  " << m_sequenceElements.size() << " sequence elements");

  size_t nWriters = 0;
  size_t nReaders = 0;
  size_t nAlgorithms = 0;
  for (const auto& alg : m_sequenceElements) {
    if (dynamic_cast<const IWriter*>(alg.get()) != nullptr) {
      nWriters++;
    } else if (dynamic_cast<const IReader*>(alg.get()) != nullptr) {
      nReaders++;
    } else if (dynamic_cast<const IAlgorithm*>(alg.get()) != nullptr) {
      nAlgorithms++;
    } else {
      throw std::runtime_error{"Unknown sequence element type"};
    }
  }

  ACTS_INFO("  " << nReaders << " readers");
  ACTS_INFO("  " << nAlgorithms << " algorithms");
  ACTS_INFO("  " << nWriters << " writers");

  ACTS_VERBOSE("Initialize sequence elements");
  for (auto& alg : m_sequenceElements) {
    ACTS_VERBOSE("Initialize " << getAlgorithmType(*alg) << ": "
                               << alg->name());
    if (alg->initialize() != ProcessCode::SUCCESS) {
      ACTS_FATAL("Failed to initialize " << getAlgorithmType(*alg) << ": "
                                         << alg->name());
      throw std::runtime_error("Failed to process event data");
    }
  }

  // execute the parallel event loop
  std::atomic<size_t> nProcessedEvents = 0;
  size_t nTotalEvents = eventsRange.second - eventsRange.first;
  m_taskArena.execute([&] {
    tbbWrap::parallel_for(
        tbb::blocked_range<size_t>(eventsRange.first, eventsRange.second),
        [&](const tbb::blocked_range<size_t>& r) {
          std::vector<Duration> localClocksAlgorithms(names.size(),
                                                      Duration::zero());

          for (size_t event = r.begin(); event != r.end(); ++event) {
            ACTS_DEBUG("start processing event " << event);
            m_cfg.iterationCallback();
            // Use per-event store
            WhiteBoard eventStore(
                Acts::getDefaultLogger("EventStore#" + std::to_string(event),
                                       m_cfg.logLevel),
                m_whiteboardObjectAliases);
            // If we ever wanted to run algorithms in parallel, this needs to
            // be changed to Algorithm context copies
            AlgorithmContext context(0, event, eventStore);
            size_t ialgo = 0;

            /// Decorate the context
            for (auto& cdr : m_decorators) {
              StopWatch sw(localClocksAlgorithms[ialgo++]);
              ACTS_VERBOSE("Execute context decorator: " << cdr->name());
              if (cdr->decorate(++context) != ProcessCode::SUCCESS) {
                throw std::runtime_error("Failed to decorate event context");
              }
            }

            ACTS_VERBOSE("Execute sequence elements");

            for (auto& alg : m_sequenceElements) {
              StopWatch sw(localClocksAlgorithms[ialgo++]);
              ACTS_VERBOSE("Execute " << getAlgorithmType(*alg) << ": "
                                      << alg->name());
              if (alg->internalExecute(++context) != ProcessCode::SUCCESS) {
                ACTS_FATAL("Failed to execute " << getAlgorithmType(*alg)
                                                << ": " << alg->name());
                throw std::runtime_error("Failed to process event data");
              }
            }

            nProcessedEvents++;
            if (logger().level() <= Acts::Logging::DEBUG) {
              ACTS_DEBUG("finished event " << event);
            } else if (nTotalEvents <= 100) {
              ACTS_INFO("finished event " << event);
            } else if (nProcessedEvents % 100 == 0) {
              ACTS_INFO(nProcessedEvents << " / " << nTotalEvents
                                         << " events processed");
            }
          }

          // add timing info to global information
          {
            tbbWrap::queuing_mutex::scoped_lock lock(clocksAlgorithmsMutex);
            for (size_t i = 0; i < clocksAlgorithms.size(); ++i) {
              clocksAlgorithms[i] += localClocksAlgorithms[i];
            }
          }
        });
  });

  ACTS_VERBOSE("Finalize sequence elements");
  for (auto& alg : m_sequenceElements) {
    ACTS_VERBOSE("Finalize " << getAlgorithmType(*alg) << ": " << alg->name());
    if (alg->finalize() != ProcessCode::SUCCESS) {
      ACTS_FATAL("Failed to finalize " << getAlgorithmType(*alg) << ": "
                                       << alg->name());
      throw std::runtime_error("Failed to process event data");
    }
  }

  // summarize timing
  Duration totalWall = Clock::now() - clockWallStart;
  Duration totalReal = std::accumulate(
      clocksAlgorithms.begin(), clocksAlgorithms.end(), Duration::zero());
  size_t numEvents = eventsRange.second - eventsRange.first;
  ACTS_INFO("Processed " << numEvents << " events in " << asString(totalWall)
                         << " (wall clock)");
  ACTS_INFO("Average time per event: " << perEvent(totalReal, numEvents));
  ACTS_DEBUG("Average time per algorithm:");
  for (size_t i = 0; i < names.size(); ++i) {
    ACTS_DEBUG("  " << names[i] << ": "
                    << perEvent(clocksAlgorithms[i], numEvents));
  }

  if (!m_cfg.outputDir.empty()) {
    storeTiming(names, clocksAlgorithms, numEvents,
                joinPaths(m_cfg.outputDir, m_cfg.outputTimingFile));
  }

  return EXIT_SUCCESS;
}

}  // namespace ActsExamples
