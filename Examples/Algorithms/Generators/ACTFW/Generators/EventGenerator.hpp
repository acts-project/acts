// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2018-03-13
/// @author Moritz Kiehn <msmk@cern.ch>

#pragma once

#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/IReader.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace FW {

/// Event generator based on separate process and vertex generators.
///
/// This must be a reader and not just an algorithm since it might read in
/// pre-computed generator samples, e.g. via HEPMC, and therefore has an
/// internal state that will be modified.
class EventGenerator final : public FW::IReader {
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
  using MultiplicityGenerator = std::function<size_t(RandomEngine&)>;
  using VertexGenerator = std::function<Acts::Vector4D(RandomEngine&)>;
  using ProcessGenerator = std::function<std::vector<SimVertex>(RandomEngine&)>;
  struct Generator {
    MultiplicityGenerator multiplicity = nullptr;
    VertexGenerator vertex = nullptr;
    ProcessGenerator process = nullptr;
  };

  struct Config {
    /// Name of the output collection
    std::string output;
    /// List of generators that should be used to generate the event
    std::vector<Generator> generators;
    /// Shuffle generated vertices to mix generator output
    bool shuffle = false;
    /// The random number service
    std::shared_ptr<const RandomNumbers> randomNumbers;
  };

  EventGenerator(const Config& cfg, Acts::Logging::Level lvl);

  /// Name of the reader.
  std::string name() const override final;
  /// Available events range. Always return [0,SIZE_MAX) since we generate them.
  std::pair<size_t, size_t> availableEvents() const final override;
  /// Generate an event.
  ProcessCode read(const AlgorithmContext& context) override final;

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace FW
