// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IReader.hpp"

#include <memory>
#include <string>
#include <unordered_map>

namespace Acts {
class Surface;
}

namespace ActsExamples {

/// Read in a planar cluster collection in comma-separated-value format.
///
/// This reads three files per event file in the configured input
/// directory. By default it reads file in the current working directory.
/// Files are assumed to be named using the following schema
///
///     event000000001-cells.csv
///     event000000001-hits.csv
///     event000000001-truth.csv
///     event000000002-cells.csv
///     event000000002-hits.csv
///     event000000002-truth.csv
///
/// and each line in the file corresponds to one hit/cluster.
class CsvPlanarClusterReader final : public IReader {
 public:
  struct Config {
    /// Where to read input files from.
    std::string inputDir;
    /// Output cluster collection.
    std::string outputClusters;
    /// For each cluster/ hit index the original hit id stored on file.
    std::string outputHitIds;
    /// Output hit-particles mapping collection.
    std::string outputMeasurementParticlesMap;
    /// Output simulated (truth) hits collection.
    std::string outputSimHits;
    /// Tracking geometry required to access global-to-local transforms.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
  };

  /// Construct the cluster reader.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  CsvPlanarClusterReader(const Config& config, Acts::Logging::Level level);

  std::string name() const override;

  /// Return the available events range.
  std::pair<size_t, size_t> availableEvents() const override;

  /// Read out data from the input stream.
  ProcessCode read(const ActsExamples::AlgorithmContext& ctx) override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::pair<size_t, size_t> m_eventsRange;
  std::unique_ptr<const Acts::Logger> m_logger;

  WriteDataHandle<GeometryIdMultimap<Acts::PlanarModuleCluster>>
      m_outputClusters{this, "OutputClusters"};
  WriteDataHandle<IndexMultimap<ActsFatras::Barcode>>
      m_outputMeasurementParticlesMap{this, "OutputMeasurementParticlesMap"};
  WriteDataHandle<SimHitContainer> m_outputSimHits{this, "OutputSimHits"};
  WriteDataHandle<std::vector<uint64_t>> m_outputHitIds{this, "OutputHitIds"};

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
