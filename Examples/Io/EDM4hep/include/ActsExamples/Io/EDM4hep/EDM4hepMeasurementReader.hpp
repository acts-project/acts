// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IReader.hpp"

#include <memory>
#include <string>

#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

namespace ActsExamples {

/// Read in a measurement cluster collection from EDM4hep.
///
/// Inpersistent information:
/// - hit index
/// - 1D local coords?
/// - segment path
///
/// Known issues:
/// - cluster channels are read from inappropriate fields
/// - local 2D coordinates and time are read from position
class EDM4hepMeasurementReader final : public IReader {
 public:
  struct Config {
    /// Where to read the input file from.
    std::string inputPath;
    /// Output measurement collection.
    std::string outputMeasurements;
    /// Output measurement to sim hit collection.
    std::string outputMeasurementSimHitsMap;
    /// Output source links collection.
    std::string outputSourceLinks;
    /// Output cluster collection (optional).
    std::string outputClusters;
  };

  /// Construct the cluster reader.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  EDM4hepMeasurementReader(const Config& config, Acts::Logging::Level level);

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

  const edm4hep::TrackerHitPlaneCollection* m_trackerHitPlaneCollection;
  const edm4hep::TrackerHitCollection* m_trackerHitRawCollection;
};

}  // namespace ActsExamples
