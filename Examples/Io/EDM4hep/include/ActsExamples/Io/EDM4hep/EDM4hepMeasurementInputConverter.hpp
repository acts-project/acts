// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"
#include "ActsExamples/Io/Podio/PodioInputConverter.hpp"

#include <memory>
#include <string>

namespace ActsExamples {

class DD4hepDetector;

/// Read in a measurement collection from EDM4hep TrackerHitLocal format.
class EDM4hepMeasurementInputConverter final : public PodioInputConverter {
 public:
  struct Config {
    /// Where to read the input frame from.
    std::string inputFrame;
    /// Name of the input tracker hit local collection.
    std::string inputTrackerHitsLocal = "ActsTrackerHitsLocal";
    /// Output measurement collection.
    std::string outputMeasurements;
    /// Output measurement to sim hit collection.
    std::string outputMeasurementSimHitsMap;
    /// Output cluster collection (optional).
    std::string outputClusters;

    /// DD4hep detector for cellID to geometry identifier resolution.
    std::shared_ptr<DD4hepDetector> dd4hepDetector;
  };

  /// Construct the cluster reader.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  explicit EDM4hepMeasurementInputConverter(
      const Config& config,
      std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Read out data from the input stream.
  ProcessCode convert(const AlgorithmContext& ctx,
                      const podio::Frame& frame) const override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  EDM4hepUtil::MapGeometryIdFrom m_geometryMapper;

  WriteDataHandle<MeasurementContainer> m_outputMeasurements{
      this, "OutputMeasurements"};

  WriteDataHandle<IndexMultimap<Index>> m_outputMeasurementSimHitsMap{
      this, "OutputMeasurementSimHitsMap"};

  WriteDataHandle<ClusterContainer> m_outputClusters{this, "OutputClusters"};
};

}  // namespace ActsExamples
