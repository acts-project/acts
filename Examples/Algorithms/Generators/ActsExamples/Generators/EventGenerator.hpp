// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Utilities/MultiplicityGenerators.hpp"
#include "ActsExamples/Utilities/ParametricParticleGenerator.hpp"
#include "ActsExamples/Utilities/VertexGenerators.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace HepMC3 {
class GenEvent;
class GenVertex;
}  // namespace HepMC3

namespace ActsExamples {

/// Event generator based on separate particles and vertex generators.
///
/// This must be a reader and not just an algorithm since it might read in
/// pre-computed generator samples, e.g. via HepMC, and therefore has an
/// internal state that will be modified.
class EventGenerator final : public ActsExamples::IReader {
 public:
  /// Combined set of generator functions.
  ///
  /// Each generator creates a number of primary vertices (multiplicity),
  /// each with an separate vertex position and time (vertex), and a set of
  /// associated particles grouped into secondary vertices (process) anchored
  /// at the primary vertex position. The first group of particles generated
  /// by the process are the particles associated directly to the primary
  /// vertex.
  ///
  /// The process generator is responsible for defining all components of the
  /// particle barcode except the primary vertex. The primary vertex will be
  /// set/overwritten by the event generator.

  /// @brief Combined struct which contains all generator components
  struct Generator {
    std::shared_ptr<MultiplicityGenerator> multiplicity = nullptr;
    std::shared_ptr<PrimaryVertexPositionGenerator> vertex = nullptr;
    std::shared_ptr<ParticlesGenerator> particles = nullptr;
  };

  struct Config {
    /// Name of the output event collection.
    std::optional<std::string> outputEvent = "hepmc3_event";

    /// List of generators that should be used to generate the event.
    std::vector<Generator> generators;
    /// The random number service.
    std::shared_ptr<const RandomNumbers> randomNumbers;

    /// If true, print the listing of the generated event. This can be very
    /// verbose
    bool printListing = false;
  };

  EventGenerator(const Config& cfg, Acts::Logging::Level lvl);

  /// Name of the reader.
  std::string name() const final;
  /// Available events range. Always return
  /// [0,std::numeric_limits<std::size_t>::max()) since we generate them.
  std::pair<std::size_t, std::size_t> availableEvents() const final;
  /// Generate an event.
  ProcessCode read(const AlgorithmContext& ctx) final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  void handleVertex(const HepMC3::GenVertex& genVertex, SimVertex& vertex,
                    std::vector<SimVertex>& vertices,
                    std::vector<SimParticle>& particles,
                    std::size_t& nSecondaryVertices, std::size_t& nParticles,
                    std::vector<bool>& seenVertices);

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  WriteDataHandle<std::shared_ptr<HepMC3::GenEvent>> m_outputEvent{
      this, "OutputEvent"};
};

}  // namespace ActsExamples
