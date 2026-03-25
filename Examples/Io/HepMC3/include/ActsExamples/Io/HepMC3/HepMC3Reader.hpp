// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Utilities/VertexGenerators.hpp"

#include <filesystem>
#include <memory>
#include <mutex>
#include <string>

namespace HepMC3 {
class GenEvent;
class Reader;
}  // namespace HepMC3

namespace ActsExamples {

struct MultiplicityGenerator;

/// HepMC3 event reader.
///
/// This reader supports reading events from one or more HepMC3 files, with
/// flexible control over how many physical events are merged into each logical
/// event.
///
/// ## Multiplicity Generators
///
/// Each input file must have an associated multiplicity generator that
/// determines how many physical events to read from that file for each logical
/// event:
/// - FixedMultiplicityGenerator: Always reads a fixed number of events
/// (deterministic)
/// - PoissonMultiplicityGenerator: Reads a Poisson-distributed number of events
/// (stochastic)
///
/// ## Usage Patterns
///
/// 1. Single file with default multiplicity:
///    ```cpp
///    HepMC3Reader::Config cfg;
///    cfg.inputPath = "events.hepmc3";
///    cfg.outputEvent = "hepmc3_event";
///    HepMC3Reader reader(cfg, Acts::Logging::INFO);
///    ```
///    Automatically uses FixedMultiplicityGenerator(n=1).
///
/// 2. Multiple files with fixed multiplicities (e.g., hard-scatter + pileup):
///    ```cpp
///    HepMC3Reader::Config cfg;
///    cfg.inputs = {
///      {.path = "signal.hepmc3",
///       .multiplicityGenerator =
///       std::make_shared<FixedMultiplicityGenerator>(1)},
///      {.path = "pileup.hepmc3",
///       .multiplicityGenerator =
///       std::make_shared<FixedMultiplicityGenerator>(50)}
///    };
///    cfg.outputEvent = "hepmc3_event";
///    HepMC3Reader reader(cfg, Acts::Logging::INFO);
///    ```
///
/// 3. With stochastic pileup multiplicity:
///    ```cpp
///    HepMC3Reader::Config cfg;
///    cfg.inputs = {
///      {.path = "signal.hepmc3",
///       .multiplicityGenerator =
///       std::make_shared<FixedMultiplicityGenerator>(1)},
///      {.path = "pileup.hepmc3",
///       .multiplicityGenerator =
///       std::make_shared<PoissonMultiplicityGenerator>(50.0)}
///    };
///    cfg.outputEvent = "hepmc3_event";
///    cfg.randomNumbers = rng; // Required for stochastic generators
///    HepMC3Reader reader(cfg, Acts::Logging::INFO);
///    ```
///
/// ## Event Skipping
///
/// Event skipping (via Sequencer's skip parameter) is only supported when ALL
/// multiplicity generators are Fixed. This is because skipping requires knowing
/// exactly how many physical events to skip, which is only deterministic with
/// FixedMultiplicityGenerator.
///
class HepMC3Reader final : public IReader {
 public:
  /// Input specification per file
  struct Input {
    /// Path to the HepMC3 file
    std::filesystem::path path;
    /// Multiplicity generator determining how many events to read per logical
    /// event. Must always be set. Use FixedMultiplicityGenerator for
    /// deterministic behavior, or PoissonMultiplicityGenerator for stochastic
    /// pileup simulation.
    std::shared_ptr<const MultiplicityGenerator> multiplicityGenerator;
  };

  struct Config {
    /// Input files to read. For each file, the multiplicity generator
    /// determines how many events are read per logical event. This can be used
    /// to read e.g. hard-scatter events from one file and pileup events from
    /// another. Mutually exclusive with inputPath.
    std::vector<Input> inputs;

    /// Convenience parameter for single-file reading. Mutually exclusive with
    /// inputs. If set, creates a single input with
    /// FixedMultiplicityGenerator(n=1).
    std::optional<std::filesystem::path> inputPath;

    /// The output collection
    std::string outputEvent;

    /// If true, print the event listing
    bool printListing = false;

    /// HepMC3 does not expose the number of events in the file, so we need to
    /// provide it here if known, otherwise the reader will have to read the
    /// whole file
    std::optional<std::size_t> numEvents = std::nullopt;

    /// If true, the reader will check if the read `HepMC3::GenEvent` has the
    /// same event number as the internal one.
    /// This will only be correct if the events were written in sequential order
    /// and numbered correctly.
    bool checkEventNumber = true;

    /// In multi-threading mode, the reader will need to buffer events to read
    /// them predictably and in order. This defines the maximum queue size being
    /// used. If this number is exceeded the reader will error out.
    std::size_t maxEventBufferSize = 128;

    /// The random number service. Required if:
    /// - vertexGenerator is set, OR
    /// - any non-Fixed multiplicityGenerator is used (e.g.,
    /// PoissonMultiplicityGenerator) Not required when using only
    /// FixedMultiplicityGenerator with no vertexGenerator.
    std::shared_ptr<const RandomNumbers> randomNumbers;
    /// Position generator that will be used to shift read events
    std::shared_ptr<PrimaryVertexPositionGenerator> vertexGenerator;
  };

  /// Construct the particle reader.
  ///
  /// @param [in] cfg The configuration object
  /// @param [in] lvl The logging level
  HepMC3Reader(const Config& cfg, Acts::Logging::Level lvl);

  ~HepMC3Reader() override;

  std::string name() const override;

  /// Return the available events range.
  std::pair<std::size_t, std::size_t> availableEvents() const override;

  /// Read out data from the input stream.
  ProcessCode read(const AlgorithmContext& ctx) override;

  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

  ProcessCode skip(std::size_t events) override;

 private:
  std::size_t determineNumEvents(const std::filesystem::path& path) const;

  static std::shared_ptr<HepMC3::GenEvent> makeEvent();

  ProcessCode readSingleFile(const AlgorithmContext& ctx,
                             std::shared_ptr<HepMC3::GenEvent>& outputEvent);

  ProcessCode readCached(
      const AlgorithmContext& ctx,
      std::vector<std::shared_ptr<HepMC3::GenEvent>>& events);

  ProcessCode readBuffer(
      const AlgorithmContext& ctx,
      std::vector<std::shared_ptr<HepMC3::GenEvent>>& outputEvents);

  ProcessCode readLogicalEvent(
      const AlgorithmContext& ctx,
      std::vector<std::shared_ptr<HepMC3::GenEvent>>& events);

  /// The configuration of this writer
  Config m_cfg;
  /// Number of events
  std::pair<std::size_t, std::size_t> m_eventsRange;
  /// The logger
  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const { return *m_logger; }

  WriteDataHandle<std::shared_ptr<HepMC3::GenEvent>> m_outputEvent{
      this, "OutputEvent"};

  std::vector<
      std::pair<std::size_t, std::vector<std::shared_ptr<HepMC3::GenEvent>>>>
      m_events;
  std::size_t m_nextEvent = 0;

  std::size_t m_maxEventBufferSize = 0;

  bool m_bufferError = false;

  std::mutex m_queueMutex;

  /// Precomputed flag: true if RNG is needed (vertex generator or non-Fixed
  /// multiplicity)
  bool m_needsRng = false;

  /// Precomputed flag: true if any input uses a non-Fixed multiplicity
  /// generator
  bool m_hasNonFixedMultiplicity = false;

  struct InputConfig {
    std::shared_ptr<HepMC3::Reader> reader;
    std::filesystem::path path;
    std::shared_ptr<const MultiplicityGenerator> multiplicityGenerator;
    std::size_t eventsRead = 0;  // Track events read from this input
  };

  std::vector<InputConfig> m_inputs;
};

}  // namespace ActsExamples
