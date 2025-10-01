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
class HepMC3Reader final : public IReader {
 public:
  struct Config {
    /// Input specification per file
    struct Input {
      /// Path to the HepMC3 file
      std::filesystem::path path;
      /// Fixed number of events to read per logical event (used if multiplicityGenerator is null)
      std::size_t numEvents = 1;
      /// Optional multiplicity generator for variable event sampling
      std::shared_ptr<const MultiplicityGenerator> multiplicityGenerator = nullptr;
    };

    /// Input files to read. For each file, a specific number of events is
    /// read per logical event. This can be used to read e.g. hard-scatter
    /// events from one file and pileup events from another.
    std::vector<Input> inputs;

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

    /// The random number service. Required if vertexGenerator or any multiplicityGenerator is set.
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
  ProcessCode read(const ActsExamples::AlgorithmContext& ctx) override;

  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

  ProcessCode skip(std::size_t events) override;

 private:
  std::size_t determineNumEvents(HepMC3::Reader& reader) const;

  static std::shared_ptr<HepMC3::GenEvent> makeEvent();

  ProcessCode readSingleFile(const ActsExamples::AlgorithmContext& ctx,
                             std::shared_ptr<HepMC3::GenEvent>& outputEvent);

  ProcessCode readCached(
      const ActsExamples::AlgorithmContext& ctx,
      std::vector<std::shared_ptr<HepMC3::GenEvent>>& events);

  ProcessCode readBuffer(
      const ActsExamples::AlgorithmContext& ctx,
      std::vector<std::shared_ptr<HepMC3::GenEvent>>& outputEvents);

  ProcessCode readLogicalEvent(
      const ActsExamples::AlgorithmContext& ctx,
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

  struct InputConfig {
    std::shared_ptr<HepMC3::Reader> reader;
    std::size_t numEvents;
    std::filesystem::path path;
  };

  std::vector<InputConfig> m_inputs;
};

}  // namespace ActsExamples
