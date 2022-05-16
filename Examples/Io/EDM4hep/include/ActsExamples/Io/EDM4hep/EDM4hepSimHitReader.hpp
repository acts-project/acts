// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IReader.hpp"

#include <memory>
#include <string>

#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

namespace ActsExamples {

class EDM4hepSimHitReader final : public IReader {
 public:
  struct Config {
    /// Where to the read input file from.
    std::string inputPath;
    /// Output simulated (truth) hits collection.
    std::string outputSimHits;
  };

  /// Construct the simhit reader.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  EDM4hepSimHitReader(const Config& config, Acts::Logging::Level level);

  std::string name() const final override;

  /// Return the available events range.
  std::pair<size_t, size_t> availableEvents() const final override;

  /// Read out data from the input stream.
  ProcessCode read(const ActsExamples::AlgorithmContext& ctx) final override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::pair<size_t, size_t> m_eventsRange;
  std::unique_ptr<const Acts::Logger> m_logger;

  podio::ROOTReader m_reader;
  podio::EventStore m_store;

  std::vector<std::string> m_simHitCollections;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
