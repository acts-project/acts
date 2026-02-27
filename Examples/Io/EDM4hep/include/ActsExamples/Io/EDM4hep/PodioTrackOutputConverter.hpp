// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/Podio/CollectionBaseWriteHandle.hpp"
#include "ActsExamples/Io/Podio/PodioOutputConverter.hpp"

#include <memory>
#include <string>

namespace ActsPlugins::PodioUtil {
class ConversionHelper;
}

namespace podio {
class CollectionBase;
}

namespace ActsPodioEdm {
class MutableTrackerHitLocal;
}

namespace ActsExamples {

/// Write out a track container to PODIO track collections
///
/// This algorithm reads a ConstTrackContainer from the whiteboard and
/// writes it to PodioTrackContainer format, which is the Acts native
/// EDM4hep format. This preserves all track states and dynamic columns,
/// unlike the standard EDM4hep format conversion.
class PodioTrackOutputConverter : public PodioOutputConverter {
 public:
  struct Config {
    /// Input track container
    std::string inputTracks;
    /// Output track collection in podio format
    std::string outputTracks = "tracks";
    /// Input measurement collection
    std::string inputMeasurements;
    /// DD4hep detector
    std::shared_ptr<DD4hepDetector> detector;
  };

  /// Constructor
  /// @param config is the configuration object
  /// @param logger is the logger
  explicit PodioTrackOutputConverter(
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
  void writeMeasurement(const Acts::GeometryContext& gctx,
                        const ConstVariableBoundMeasurementProxy& meas,
                        ActsPodioEdm::MutableTrackerHitLocal& to) const;

  const Acts::Surface* surfaceByIdentifier(
      Acts::GeometryIdentifier geometryId) const;

  static dd4hep::CellID cellIdFromSurface(const Acts::Surface& surface);

  Config m_cfg;

  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};

  CollectionBaseWriteHandle m_outputTracks{this, "OutputTracks"};
  CollectionBaseWriteHandle m_outputTrackStates{this, "OutputTrackStates"};
  CollectionBaseWriteHandle m_outputMeasurements{this, "OutputMeasurements"};
  CollectionBaseWriteHandle m_outputTrackStateHitLinks{
      this, "OutputTrackStateHitLinks"};
};

}  // namespace ActsExamples
