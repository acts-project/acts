// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Io/Podio/PodioCollectionDataHandle.hpp"
#include "ActsExamples/Io/Podio/PodioOutputConverter.hpp"
#include "ActsPlugins/EDM4hep/EDM4hepUtil.hpp"

#include <edm4hep/SimTrackerHit.h>

namespace podio {
class CollectionBase;
}

namespace ActsPodioEdm {
class MeasurementCollection;
}

namespace ActsExamples {

class MeasurementContainer;

/// Write out a track collection to EDM4hep objects
class PodioMeasurementOutputConverter : public PodioOutputConverter {
 public:
  struct Config {
    /// Input measurement component
    std::string inputMeasurements;
    /// Output podio measurement collection
    std::string outputMeasurements = "ActsMeasurements";
    /// Input simhit association (optional)
    std::optional<std::string> inputSimHitAssociation = std::nullopt;
    /// Input collection to map measured hits to simulated hits (optional)
    std::optional<std::string> inputMeasurementSimHitsMap = std::nullopt;
  };

  explicit PodioMeasurementOutputConverter(
      const Config& config,
      std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

  std::vector<std::string> collections() const final;

 protected:
  /// @brief Write method called by the base class
  /// @param [in] context is the algorithm context for consistency
  ProcessCode execute(const AlgorithmContext& context) const final;

 private:
  Config m_cfg;

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};

  ReadDataHandle<ActsPlugins::EDM4hepUtil::SimHitAssociation>
      m_inputSimHitAssociation{this, "InputSimHitAssociation"};

  ReadDataHandle<IndexMultimap<Index>> m_inputMeasurementSimHitsMap{
      this, "InputMeasurementSimHitsMap"};

  PodioCollectionWriteHandle<ActsPodioEdm::MeasurementCollection>
      m_outputMeasurements{this, "OutputMeasurements"};
};

}  // namespace ActsExamples
