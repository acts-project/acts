// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <exception>
#include <future>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string_view>
#include <thread>

#include <G4MTRunManager.hh>
#include <G4Run.hh>
#include <Randomize.hh>
#include <tbb/global_control.h>

#include "NativeWorkerExchange.hpp"
#include "NativeWorkerTransport.hpp"

namespace {
using namespace ActsExamples;
using namespace ActsTests::NativeWorkers;

constexpr int firstEvent = 17;
constexpr int nEvents = 4;
constexpr int badEvent = firstEvent + 2;

Request input(const AlgorithmContext& ctx) {
  const int id = static_cast<int>(ctx.eventNumber);
  const RandomNumbers rng({42});
  // Use a positive seed representable by the supported CLHEP engines.
  return {id, static_cast<long>(1 + rng.generateSeed(ctx) % 900000000),
          2. * (id - firstEvent)};
}

struct Track {
  double random;
  double intercept;
  double slope;
};

// Hold one track inside its first silicon step until reconstruction of a
// different event finishes. No ordering inference depends on wall-clock times.
struct Observations {
  explicit Observations(bool overlap) : testOverlap(overlap) {}

  void transportStep(int id) {
    std::unique_lock lock(mutex);
    workers.insert(std::this_thread::get_id());
    if (!testOverlap) {
      return;
    }
    if (id == firstEvent) {
      firstAtSensor = true;
      changed.notify_all();
    }
    if (!changed.wait_for(lock, std::chrono::seconds(10), [&] {
          return id == firstEvent ? !tracks.empty() : firstAtSensor;
        })) {
      throw std::runtime_error(
          "Timed out waiting for transport/reconstruction overlap");
    }
  }

  std::mutex mutex;
  std::condition_variable changed;
  const bool testOverlap;
  bool firstAtSensor = false;
  std::set<std::thread::id> workers;
  std::map<int, Track> tracks;
  bool reconstructedDuringRun = false;
  std::atomic<bool> beamReturned{false};
};

class Transport final : public IAlgorithm {
 public:
  Transport(Exchange& exchange, bool injectFailure)
      : IAlgorithm("Transport"),
        m_exchange(exchange),
        m_injectFailure(injectFailure) {
    m_hits.initialize("transport");
  }
  ProcessCode execute(const AlgorithmContext& ctx) const override {
    auto request = input(ctx);
    request.fail = m_injectFailure && request.eventId == badEvent;
    auto response = m_exchange.submit(request);
    try {
      m_hits(ctx, response.get());
    } catch (const TransportFailure&) {
      if (request.fail) {
        return ProcessCode::SKIP;
      }
      throw;
    }
    return ProcessCode::SUCCESS;
  }

 private:
  Exchange& m_exchange;
  bool m_injectFailure;
  WriteDataHandle<Result> m_hits{this, "TransportOutput"};
};

class Reconstruct final : public IAlgorithm {
 public:
  explicit Reconstruct(Observations& observed)
      : IAlgorithm("Reconstruct"), m_observed(observed) {
    m_hits.initialize("transport");
  }
  ProcessCode execute(const AlgorithmContext& ctx) const override {
    const auto& result = m_hits(ctx);
    if (result.eventId != static_cast<int>(ctx.eventNumber) ||
        result.hits.size() != 2) {
      throw std::runtime_error("Wrong transport result on the ACTS event");
    }
    const auto [x1, z1] = result.hits[0];
    const auto [x2, z2] = result.hits[1];
    if (std::abs(z1 - 29.75) > 1e-6 || std::abs(z2 - 59.75) > 1e-6) {
      throw std::runtime_error("Missing silicon entry crossings");
    }
    const double slope = (x2 - x1) / (z2 - z1);
    std::lock_guard lock(m_observed.mutex);
    if (!m_observed.beamReturned) {
      m_observed.reconstructedDuringRun = true;
    }
    m_observed.tracks.emplace(result.eventId,
                              Track{result.random, x1 - slope * z1, slope});
    m_observed.changed.notify_all();
    return ProcessCode::SUCCESS;
  }

 private:
  Observations& m_observed;
  ReadDataHandle<Result> m_hits{this, "TransportInput"};
};
}  // namespace

BOOST_AUTO_TEST_CASE(HandOffEventsToNativeWorkers) {
  const auto& suite = boost::unit_test::framework::master_test_suite();
  const std::string_view mode = suite.argc > 1 ? suite.argv[1] : "overlap";
  const bool earlyEnd = mode == "short";
  const int nWorkers = mode == "reference" ? 1 : 2;
  const int nCallers = earlyEnd ? 4 : 2;
  // Allow concurrent ACTS callers even on a single-core test runner.
  tbb::global_control concurrency(tbb::global_control::max_allowed_parallelism,
                                  nCallers);
  const int budget = earlyEnd ? 64 : nEvents;
  Observations observed(mode == "overlap");
  Exchange exchange(2);
  G4MTRunManager manager;
  manager.SetNumberOfThreads(nWorkers);
  BOOST_REQUIRE_EQUAL(manager.GetNumberOfThreads(), nWorkers);
  configure(manager, exchange, [&](int id) { observed.transportStep(id); });
  manager.Initialize();

  const auto engineState = G4Random::getTheEngine()->put();
  auto acts = std::async(std::launch::async, [&] {
    try {
      Sequencer::Config cfg;
      cfg.skip = firstEvent;
      cfg.events = nEvents;
      cfg.numThreads = nCallers;
      cfg.trackFpes = false;
      cfg.logLevel = Acts::Logging::WARNING;
      Sequencer sequencer(cfg);
      sequencer.addAlgorithm(std::make_shared<Transport>(exchange, earlyEnd));
      sequencer.addAlgorithm(std::make_shared<Reconstruct>(observed));
      const int status = sequencer.run();
      exchange.close();
      return status;
    } catch (...) {
      exchange.cancel(std::current_exception());
      throw;
    }
  });

  try {
    manager.BeamOn(budget);
  } catch (...) {
    exchange.cancel(std::current_exception());
    throw;
  }
  observed.beamReturned = true;
  // Wake callers even if the master consumed its event budget prematurely.
  exchange.cancel(
      std::make_exception_ptr(std::runtime_error("Geant4 run ended")));
  BOOST_REQUIRE_EQUAL(acts.get(), 0);
  BOOST_CHECK(observed.reconstructedDuringRun);
  BOOST_CHECK_LE(exchange.maxOutstanding(), 2u);
  if (observed.testOverlap) {
    BOOST_CHECK_EQUAL(observed.workers.size(), 2u);
    BOOST_CHECK_EQUAL(exchange.maxOutstanding(), 2u);
  }
  BOOST_REQUIRE_EQUAL(observed.tracks.size(),
                      nEvents - static_cast<int>(earlyEnd));

  // Both worker counts must reproduce the same per-event random draw and
  // fitted track. Compute the reference independently on the test thread.
  std::unique_ptr<CLHEP::HepRandomEngine> reference(
      CLHEP::HepRandomEngine::newEngine(engineState));
  BOOST_REQUIRE(reference);
  WhiteBoard board;
  for (int id = firstEvent; id < firstEvent + nEvents; ++id) {
    if (earlyEnd && id == badEvent) {
      BOOST_CHECK(!observed.tracks.contains(id));
      continue;
    }
    const auto request = input({1, static_cast<std::size_t>(id), board, 0});
    reference->setSeed(request.seed, 3);
    const double expected = reference->flat();
    const auto& track = observed.tracks.at(id);
    BOOST_CHECK_EQUAL(track.random, expected);
    BOOST_CHECK_SMALL(track.intercept - (request.x + expected), 1e-8);
    BOOST_CHECK_SMALL(track.slope, 1e-9);
  }
  BOOST_REQUIRE(manager.GetCurrentRun() != nullptr);
  const int slots = manager.GetCurrentRun()->GetNumberOfEvent();
  BOOST_TEST_MESSAGE("ACTS submitted "
                     << nEvents << ", reconstructed " << observed.tracks.size()
                     << "; Geant4 processed " << slots << " of " << budget
                     << " slots");
  BOOST_CHECK_GE(slots, nEvents);
  BOOST_CHECK_LE(slots, earlyEnd ? nEvents + nWorkers : nEvents);
}
