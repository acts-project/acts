// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/Framework/IReader.hpp"

#include <memory>
#include <string>

#include "edm4hep/MCParticleCollection.h"
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

namespace ActsExamples {

/// Read in a simhit collection from EDM4hep.
///
/// Inpersistent information:
/// - particle ID
/// - after4 momentum
/// - hit index
/// - digitization channel
class EDM4hepSimHitReader final : public IReader {
 public:
  struct Config {
    /// Where to the read input file from.
    std::string inputPath;
    /// Name of the particle collection in EDM4hep.
    std::string inputParticles = "MCParticles";
    /// Which particle collection to read into.
    std::string outputParticles;
    /// Output simulated (truth) hits collection.
    std::string outputSimHits;
    /// DD4hep detector for cellID resolution.
    std::shared_ptr<DD4hep::DD4hepDetector> dd4hepDetector;
  };

  /// Construct the simhit reader.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  EDM4hepSimHitReader(const Config& config, Acts::Logging::Level level);

  std::string name() const final;

  /// Return the available events range.
  std::pair<size_t, size_t> availableEvents() const final;

  /// Read out data from the input stream.
  ProcessCode read(const ActsExamples::AlgorithmContext& ctx) final;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::pair<size_t, size_t> m_eventsRange;
  std::unique_ptr<const Acts::Logger> m_logger;

  podio::ROOTReader m_reader;
  podio::EventStore m_store;

  std::vector<std::string> m_collections;

  const edm4hep::MCParticleCollection* m_mcParticleCollection;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
