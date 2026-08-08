// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// Throughput benchmark for the synthetic event generator itself.
///
/// The seeding benchmarks time a seeder on a generated event. This one times
/// the generation, so that a change to the generator can be told apart from a
/// change to the seeder. What is timed is one event, single threaded; the
/// detector is built once outside the timed region, being per-job work.
///
/// The event multiplicity is reported next to the time, because the two have to
/// be read together: generating ten percent faster while producing ten percent
/// fewer space points has made nothing faster.
///
/// The cost is the helix propagation, at 156 surfaces and some 52k tracks per
/// event; `--reuse` measures no better than the default, so the container
/// allocation is not a measurable part of it. At a fixed space point count, six
/// times the surfaces come to only 12% more time, the surface loop being an
/// intersection with an early exit.

#include "Acts/Definitions/Units.hpp"
#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/EventCsvWriter.hpp"
#include "ActsFatras/Synthetic/EventGenerator.hpp"
#include "ActsFatras/Synthetic/SeedingTruth.hpp"
#include "ActsTests/CommonHelpers/BenchmarkTools.hpp"

#include <chrono>
#include <cstddef>
#include <iostream>
#include <string>

#include <boost/program_options.hpp>

#include "SyntheticEventOptions.hpp"

namespace po = boost::program_options;

using namespace Acts::UnitLiterals;
using namespace ActsTests;
using namespace ActsFatras::Synthetic;

int main(int argc, char* argv[]) {
  SyntheticEventOptions::Values shared;
  std::size_t numRuns = 20;
  std::size_t warmupMs = 2000;
  float truthPt = 1000.f;
  bool reuse = false;
  std::string dumpPrefix;

  try {
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message");
    SyntheticEventOptions::add(desc, shared);
    desc.add_options()("runs",
                       po::value<std::size_t>(&numRuns)->default_value(numRuns),
                       "number of events to time")(
        "warmup", po::value<std::size_t>(&warmupMs)->default_value(warmupMs),
        "warm-up time in ms before the timed events")(
        "truth-pt", po::value<float>(&truthPt)->default_value(truthPt),
        "momentum in MeV a primary has to reach to be counted as seedable")(
        "reuse", po::bool_switch(&reuse),
        "reuse the event storage between events instead of starting from an "
        "empty one, which is what an experiment with a per-event arena would "
        "see")("dump", po::value<std::string>(&dumpPrefix),
               "write one event to <prefix>_spacepoints.csv and "
               "<prefix>_particles.csv before timing anything");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help") > 0) {
      std::cout << desc << std::endl;
      return 0;
    }
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }

  DetectorLayout layout;
  EventConfig config;
  try {
    layout = SyntheticEventOptions::makeLayout(shared);
    config = SyntheticEventOptions::makeConfig(shared);
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }

  // microBenchmark takes quartiles over the runs, which needs at least two
  if (numRuns < 2) {
    std::cerr << "error: --runs needs at least two events" << std::endl;
    return 1;
  }

  EventGenerator generator(layout, config);

  // Written before anything is timed, so that a dump costs one event rather
  // than a full benchmark.
  if (!dumpPrefix.empty()) {
    const Event dumped = generator.generate();
    writeEventCsv(dumped, layout, dumpPrefix);
    generator.reset(config.seed);
    std::cout << "wrote " << dumpPrefix << "_spacepoints.csv and " << dumpPrefix
              << "_particles.csv" << std::endl;
  }

  std::cout << "layout=" << shared.layout
            << " surfaces=" << layout.surfaces.size()
            << " layers=" << layout.layers.size()
            << " pileup=" << config.generation.pileup
            << " electronRate=" << config.simulation.secondaries.electronRate
            << " nuclearRate=" << config.simulation.secondaries.nuclearRate
            << " seed=" << config.seed << " reuse=" << (reuse ? 1 : 0)
            << std::endl;

  // Reused across events when --reuse is given: `SpacePointContainer::clear`
  // keeps the columns and their capacity.
  Event reusedEvent;

  const auto generateOne = [&] {
    if (reuse) {
      generator.generate(reusedEvent);
      assumeRead(reusedEvent);
    } else {
      Event event = generator.generate();
      assumeRead(event);
    }
  };

  // Warm up by hand: microBenchmark's own warm-up runs a wall-clock-bounded
  // number of iterations, which would leave the generator at an unpredictable
  // point in its random sequence. Restarting the sequence afterwards keeps the
  // timed events at 0 to numRuns-1 while the caches have seen realistic sizes.
  {
    const auto warmupStart = std::chrono::steady_clock::now();
    const auto warmupTime = std::chrono::milliseconds(warmupMs);
    while (std::chrono::steady_clock::now() - warmupStart < warmupTime) {
      generateOne();
    }
  }
  generator.reset(config.seed);

  // one event per run, so the reported median and error are per event
  const auto result =
      microBenchmark(generateOne, 1, numRuns, std::chrono::milliseconds(0));

  // Replay the same events outside the timed region to describe them. Doing it
  // inside would put the truth bookkeeping into the measurement.
  generator.reset(config.seed);
  EventSummary total;
  Event event;
  for (std::size_t i = 0; i < numRuns; ++i) {
    generator.generate(event);
    const EventSummary summary = summarize(event, truthPt * 1_MeV / 1_GeV);
    total.spacePoints += summary.spacePoints;
    total.primaries += summary.primaries;
    total.secondaries += summary.secondaries;
    total.seedablePrimaries += summary.seedablePrimaries;
    total.primaryHits += summary.primaryHits;
    total.secondaryHits += summary.secondaryHits;
  }

  const auto perEvent = [&](std::size_t value) { return value / numRuns; };
  std::cout << "mean/event: spacePoints=" << perEvent(total.spacePoints)
            << " primaries=" << perEvent(total.primaries)
            << " secondaries=" << perEvent(total.secondaries)
            << " primaryHits=" << perEvent(total.primaryHits)
            << " secondaryHits=" << perEvent(total.secondaryHits)
            << " seedable=" << perEvent(total.seedablePrimaries) << std::endl;

  // one event per run, so the median run time is already the time per event
  const double msPerEvent = result.runTimeMedian().count() / 1e6;
  std::cout << "per event: " << msPerEvent << " ms -> " << 1000. / msPerEvent
            << " events/s, "
            << static_cast<double>(perEvent(total.spacePoints)) /
                   (msPerEvent * 1000.)
            << " Mspacepoints/s" << std::endl;
  std::cout << result << std::endl;

  return 0;
}
