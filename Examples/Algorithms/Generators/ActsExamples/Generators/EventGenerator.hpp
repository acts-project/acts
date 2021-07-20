// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2018-03-13
/// @author Moritz Kiehn <msmk@cern.ch>

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <functional>
#include <memory>
#include <utility>
#include <vector>

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
  // using MultiplicityGenerator = std::function<size_t(RandomEngine&)>;
  // using VertexGenerator = std::function<Acts::Vector4(RandomEngine&)>;
  // using ParticlesGenerator =
  // std::function<SimParticleContainer(RandomEngine&)>;

  struct MultiplicityGenerator {
    virtual ~MultiplicityGenerator() = default;
    virtual size_t operator()(RandomEngine& rng) const = 0;
  };

  struct VertexGenerator {
    virtual ~VertexGenerator() = default;
    virtual Acts::Vector4 operator()(RandomEngine& rng) const = 0;
  };

  struct ParticlesGenerator {
    virtual ~ParticlesGenerator() = default;
    virtual SimParticleContainer operator()(RandomEngine& rng) = 0;
  };

  struct Generator {
    std::shared_ptr<MultiplicityGenerator> multiplicity = nullptr;
    std::shared_ptr<VertexGenerator> vertex = nullptr;
    std::shared_ptr<ParticlesGenerator> particles = nullptr;
  };

  struct Config {
    /// Name of the output particles collection.
    std::string outputParticles;
    /// List of generators that should be used to generate the event.
    std::vector<Generator> generators;
    /// The random number service.
    std::shared_ptr<const RandomNumbers> randomNumbers;
  };

  EventGenerator(const Config& cfg, Acts::Logging::Level lvl);

  /// Name of the reader.
  std::string name() const final;
  /// Available events range. Always return [0,SIZE_MAX) since we generate them.
  std::pair<size_t, size_t> availableEvents() const final;
  /// Generate an event.
  ProcessCode read(const AlgorithmContext& context) final;

  const Config& config() const { return m_cfg; }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace ActsExamples
