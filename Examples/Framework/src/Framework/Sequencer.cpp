// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Framework/Sequencer.hpp"

#include <TROOT.h>
#include <algorithm>
#include <chrono>
#include <dfe/dfe_io_dsv.hpp>
#include <dfe/dfe_namedtuple.hpp>
#include <exception>
#include <numeric>
#include <tbb/tbb.h>

#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Utilities/Paths.hpp"

FW::Sequencer::Sequencer(const Sequencer::Config& cfg)
    : m_cfg(cfg),
      m_logger(Acts::getDefaultLogger("Sequencer", m_cfg.logLevel)) {
  // automatically determine the number of concurrent threads to use
  if (m_cfg.numThreads < 0) {
    m_cfg.numThreads = tbb::task_scheduler_init::default_num_threads();
  }
  ROOT::EnableThreadSafety();
}

void FW::Sequencer::addService(std::shared_ptr<IService> service) {
  if (not service) {
    throw std::invalid_argument("Can not add empty/NULL service");
  }
  m_services.push_back(std::move(service));
  ACTS_INFO("Added service '" << m_services.back()->name() << "'");
}

void FW::Sequencer::addContextDecorator(
    std::shared_ptr<IContextDecorator> decorator) {
  if (not decorator) {
    throw std::invalid_argument("Can not add empty/NULL context decorator");
  }
  m_decorators.push_back(std::move(decorator));
  ACTS_INFO("Added context decarator '" << m_decorators.back()->name() << "'");
}

void FW::Sequencer::addReader(std::shared_ptr<IReader> reader) {
  if (not reader) {
    throw std::invalid_argument("Can not add empty/NULL reader");
  }
  m_readers.push_back(std::move(reader));
  ACTS_INFO("Added reader '" << m_readers.back()->name() << "'");
}

void FW::Sequencer::addAlgorithm(std::shared_ptr<IAlgorithm> algorithm) {
  if (not algorithm) {
    throw std::invalid_argument("Can not add empty/NULL algorithm");
  }
  m_algorithms.push_back(std::move(algorithm));
  ACTS_INFO("Added algorithm '" << m_algorithms.back()->name() << "'");
}

void FW::Sequencer::addWriter(std::shared_ptr<IWriter> writer) {
  if (not writer) {
    throw std::invalid_argument("Can not add empty/NULL writer");
  }
  m_writers.push_back(std::move(writer));
  ACTS_INFO("Added writer '" << m_writers.back()->name() << "'");
}

std::vector<std::string> FW::Sequencer::listAlgorithmNames() const {
  std::vector<std::string> names;

  // WARNING this must be done in the same order as in the processing
  for (const auto& service : m_services) {
    names.push_back("Service:" + service->name());
  }
  for (const auto& decorator : m_decorators) {
    names.push_back("Decorator:" + decorator->name());
  }
  for (const auto& reader : m_readers) {
    names.push_back("Reader:" + reader->name());
  }
  for (const auto& algorithm : m_algorithms) {
    names.push_back("Algorithm:" + algorithm->name());
  }
  for (const auto& writer : m_writers) {
    names.push_back("Writer:" + writer->name());
  }

  return names;
}

namespace {
// Saturated addition that does not overflow and exceed SIZE_MAX.
//
// From http://locklessinc.com/articles/sat_arithmetic/
size_t saturatedAdd(size_t a, size_t b) {
  size_t res = a + b;
  res |= -(res < a);
  return res;
}
}  // namespace

std::pair<std::size_t, std::size_t> FW::Sequencer::determineEventsRange()
    const {
  constexpr auto kInvalidEventsRange = std::make_pair(SIZE_MAX, SIZE_MAX);

  // Note on skipping events:
  //
  // Previously, skipping events was only allowed when readers where available,
  // since only readers had a `.skip()` functionality. The `.skip()` interface
  // has been removed in favour of telling the readers the event they are
  // requested to read via the algorithm context.
  // Skipping can now also be used when no readers are configured, e.g. for
  // generating only a few specific events in a simulation setup.

  // determine intersection of event ranges available from readers
  size_t beg = 0u;
  size_t end = SIZE_MAX;
  for (const auto& reader : m_readers) {
    auto available = reader->availableEvents();
    beg = std::max(beg, available.first);
    end = std::min(end, available.second);
  }

  // since we use event ranges (and not just num events) they might not overlap
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
  if ((beg == 0u) and (end == SIZE_MAX) and (m_cfg.events == SIZE_MAX)) {
    ACTS_ERROR("Could not determine number of events");
    return kInvalidEventsRange;
  }

  // take user selection into account
  auto begSelected = saturatedAdd(beg, m_cfg.skip);
  auto endRequested = saturatedAdd(begSelected, m_cfg.events);
  auto endSelected = std::min(end, endRequested);
  if (end < endRequested) {
    ACTS_INFO("Restrict requested number of events to available ones");
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
  double time_total_s;
  double time_perevent_s;

  DFE_NAMEDTUPLE(TimingInfo, identifier, time_total_s, time_perevent_s);
};

void storeTiming(const std::vector<std::string>& identifiers,
                 const std::vector<Duration>& durations, std::size_t numEvents,
                 std::string path) {
  dfe::NamedTupleTsvWriter<TimingInfo> writer(std::move(path), 4);
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

int FW::Sequencer::run() {
  // measure overall wall clock
  Timepoint clockWallStart = Clock::now();
  // per-algorithm time measures
  std::vector<std::string> names = listAlgorithmNames();
  std::vector<Duration> clocksAlgorithms(names.size(), Duration::zero());
  tbb::queuing_mutex clocksAlgorithmsMutex;

  // processing only works w/ a well-known number of events
  // error message is already handled by the helper function
  std::pair<size_t, size_t> eventsRange = determineEventsRange();
  if ((eventsRange.first == SIZE_MAX) and (eventsRange.second == SIZE_MAX)) {
    return EXIT_FAILURE;
  }

  ACTS_INFO("Processing events [" << eventsRange.first << ", "
                                  << eventsRange.second << ")");
  ACTS_INFO("Starting event loop with " << m_cfg.numThreads << " threads");
  ACTS_INFO("  " << m_services.size() << " services");
  ACTS_INFO("  " << m_decorators.size() << " context decorators");
  ACTS_INFO("  " << m_readers.size() << " readers");
  ACTS_INFO("  " << m_algorithms.size() << " algorithms");
  ACTS_INFO("  " << m_writers.size() << " writers");

  // run start-of-run hooks
  for (auto& service : m_services) {
    names.push_back("Service:" + service->name() + ":startRun");
    clocksAlgorithms.push_back(Duration::zero());
    StopWatch sw(clocksAlgorithms.back());
    service->startRun();
  }

  // execute the parallel event loop
  tbb::task_scheduler_init init(m_cfg.numThreads);
  tbb::parallel_for(
      tbb::blocked_range<size_t>(eventsRange.first, eventsRange.second),
      [&](const tbb::blocked_range<size_t>& r) {
        std::vector<Duration> localClocksAlgorithms(names.size(),
                                                    Duration::zero());

        for (size_t event = r.begin(); event != r.end(); ++event) {
          // Use per-event store
          WhiteBoard eventStore(Acts::getDefaultLogger(
              "EventStore#" + std::to_string(event), m_cfg.logLevel));
          // If we ever wanted to run algorithms in parallel, this needs to be
          // changed to Algorithm context copies
          AlgorithmContext context(0, event, eventStore);
          size_t ialgo = 0;

          // Prepare event store w/ service information
          for (auto& service : m_services) {
            StopWatch sw(localClocksAlgorithms[ialgo++]);
            service->prepare(++context);
          }
          /// Decorate the context
          for (auto& cdr : m_decorators) {
            StopWatch sw(localClocksAlgorithms[ialgo++]);
            if (cdr->decorate(++context) != ProcessCode::SUCCESS) {
              throw std::runtime_error("Failed to decorate event context");
            }
          }
          // Read everything in
          for (auto& rdr : m_readers) {
            StopWatch sw(localClocksAlgorithms[ialgo++]);
            if (rdr->read(++context) != ProcessCode::SUCCESS) {
              throw std::runtime_error("Failed to read input data");
            }
          }
          // Execute all algorithms
          for (auto& alg : m_algorithms) {
            StopWatch sw(localClocksAlgorithms[ialgo++]);
            if (alg->execute(++context) != ProcessCode::SUCCESS) {
              throw std::runtime_error("Failed to process event data");
            }
          }
          // Write out results
          for (auto& wrt : m_writers) {
            StopWatch sw(localClocksAlgorithms[ialgo++]);
            if (wrt->write(++context) != ProcessCode::SUCCESS) {
              throw std::runtime_error("Failed to write output data");
            }
          }
          ACTS_INFO("finished event " << event);
        }

        // add timing info to global information
        {
          tbb::queuing_mutex::scoped_lock lock(clocksAlgorithmsMutex);
          for (size_t i = 0; i < clocksAlgorithms.size(); ++i) {
            clocksAlgorithms[i] += localClocksAlgorithms[i];
          }
        }
      });

  // run end-of-run hooks
  for (auto& wrt : m_writers) {
    names.push_back("Writer:" + wrt->name() + ":endRun");
    clocksAlgorithms.push_back(Duration::zero());
    StopWatch sw(clocksAlgorithms.back());
    if (wrt->endRun() != ProcessCode::SUCCESS) {
      return EXIT_FAILURE;
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
  storeTiming(names, clocksAlgorithms, numEvents,
              joinPaths(m_cfg.outputDir, "timing.tsv"));

  return EXIT_SUCCESS;
}
