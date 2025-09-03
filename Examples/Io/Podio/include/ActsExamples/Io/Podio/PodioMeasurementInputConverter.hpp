// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Io/Podio/PodioInputConverter.hpp"
#include "ActsPlugins/EDM4hep/EDM4hepUtil.hpp"

namespace ActsExamples {

class PodioMeasurementInputConverter : public PodioInputConverter {
 public:
  struct Config {
    std::string inputFrame;
    /// Input measurement collection name in podio
    std::string inputMeasurements = "ActsMeasurements";
    /// Output measurement collection
    // std::string outputMeasurements;
  };

  /// constructor
  /// @param config is the configuration object
  /// @param level is the output logging level
  explicit PodioMeasurementInputConverter(
      const Config& config, Acts::Logging::Level level = Acts::Logging::INFO);

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

  ProcessCode convert(const AlgorithmContext& ctx,
                      const podio::Frame& frame) const final;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
