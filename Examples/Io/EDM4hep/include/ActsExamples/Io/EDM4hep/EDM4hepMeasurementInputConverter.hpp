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
#include "ActsExamples/Io/EDM4hep/EDM4hepInputConverter.hpp"

#include <string>

namespace podio {
class Frame;
}

namespace ActsExamples {

/// Read in a measurement cluster collection as EDM4hep from a @c podio::Frame.
///
/// Inpersistent information:
/// - hit index
/// - 1D local coords?
/// - segment path
///
/// Known issues:
/// - cluster channels are read from inappropriate fields
/// - local 2D coordinates and time are read from position
class EDM4hepMeasurementInputConverter final : public EDM4hepInputConverter {
 public:
  struct Config {
    /// Where to read the input frame from.
    std::string inputFrame;
    /// Output measurement collection.
    std::string outputMeasurements;
    /// Output measurement to sim hit collection.
    std::string outputMeasurementSimHitsMap;
    /// Output cluster collection (optional).
    std::string outputClusters;
  };

  /// Construct the cluster reader.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  EDM4hepMeasurementInputConverter(const Config& config,
                                   Acts::Logging::Level level);

  /// Read out data from the input stream.
  ProcessCode convert(const AlgorithmContext& ctx,
                      const podio::Frame& frame) const override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<podio::Frame> m_inputFrame{this, "InputFrame"};

  WriteDataHandle<MeasurementContainer> m_outputMeasurements{
      this, "OutputMeasurements"};

  WriteDataHandle<IndexMultimap<Index>> m_outputMeasurementSimHitsMap{
      this, "OutputMeasurementSimHitsMap"};

  WriteDataHandle<ClusterContainer> m_outputClusters{this, "OutputClusters"};
};

}  // namespace ActsExamples
