// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <limits>
#include <memory>
#include <string>

namespace ActsExamples {

/// Write out a planar cluster collection in comma-separated-value format.
///
/// This writes multiples file per event containing information about the
/// space points, local constituent cells, and hit-particle truth mapping
/// into the configured output directory. By default it writes to the
/// current working directory. Files are named using the following schema
///
///     event000000001-cells.csv
///     event000000001-hits.csv
///     event000000001-truth.csv
///     event000000002-cells.csv
///     event000000002-hits.csv
///     event000000002-truth.csv
///     ...
///
/// and each line in the file corresponds to one hit/cluster.
class CsvPlanarClusterWriter final
    : public WriterT<GeometryIdMultimap<Acts::PlanarModuleCluster>> {
 public:
  struct Config {
    /// Which cluster collection to write.
    std::string inputClusters;
    /// Which simulated (truth) hits collection to use.
    std::string inputSimHits;
    /// Where to place output files
    std::string outputDir;
    /// Number of decimal digits for floating point precision in output.
    size_t outputPrecision = std::numeric_limits<float>::max_digits10;
    /// Tracking geometry required to access global-to-local transforms.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
  };

  /// Construct the cluster writer.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  CsvPlanarClusterWriter(const Config& config, Acts::Logging::Level level);

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 protected:
  /// Type-specific write implementation.
  ///
  /// @param[in] ctx is the algorithm context
  /// @param[in] particles are the particle to be written
  ProcessCode writeT(
      const AlgorithmContext& ctx,
      const GeometryIdMultimap<Acts::PlanarModuleCluster>& clusters) override;

 private:
  Config m_cfg;

  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
};

}  // namespace ActsExamples
