// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IReader.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <memory>
#include <string>

#include <edm4hep/MCParticleCollection.h>
#include <podio/ROOTFrameReader.h>

namespace ActsExamples {

/// Read particles from EDM4hep.
///
/// Inpersistent information:
/// - particle ID
/// - process
class EDM4hepReader final : public IReader {
 public:
  struct Config {
    /// Where to read input file from.
    std::string inputPath;
    /// Name of the particle collection in EDM4hep.
    std::string inputParticles = "MCParticles";
    /// Which particle collection to read into.
    std::string outputParticles;
  };

  /// Construct the particle reader.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  EDM4hepReader(const Config& config, Acts::Logging::Level level);

  std::string name() const final;

  /// Return the available events range.
  std::pair<std::size_t, std::size_t> availableEvents() const final;

  /// Read out data from the input stream.
  ProcessCode read(const ActsExamples::AlgorithmContext& ctx) final;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

  void processChildren(const edm4hep::MCParticle& particle, SimBarcode parentId,
                       SimParticleContainer::sequence_type& particles,
                       std::size_t& nSecondaryVertices,
                       std::size_t& maxGen) const;

  static void setSubParticleIds(
      SimParticleContainer::sequence_type::iterator begin,
      SimParticleContainer::sequence_type::iterator end);

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::pair<std::size_t, std::size_t> m_eventsRange;
  std::unique_ptr<const Acts::Logger> m_logger;

  podio::ROOTFrameReader m_reader;

  WriteDataHandle<SimParticleContainer> m_outputParticles{this,
                                                          "OutputParticles"};
};

}  // namespace ActsExamples
