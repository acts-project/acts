// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/Podio/PodioCollectionDataHandle.hpp"
#include "ActsExamples/Io/Podio/PodioOutputConverter.hpp"

#include <string>

namespace edm4hep {
class TrackerHitPlaneCollection;
class TrackerHit3DCollection;
}  // namespace edm4hep

namespace ActsExamples {

/// Write out a measurement cluster collection to EDM4hep objects
///
/// Inpersistent information:
/// - hit index
/// - 1D local coords?
/// - segment path
///
/// Known issues:
/// - cluster channels are written to inappropriate fields
/// - local 2D coordinates and time are written to position
class EDM4hepMeasurementOutputConverter final : public PodioOutputConverter {
 public:
  struct Config {
    /// Which measurement collection to write.
    std::string inputMeasurements;
    /// Which cluster collection to write (optional)
    std::string inputClusters;
    /// Name of the output tracker hit plane collection.
    std::string outputTrackerHitsPlane = "ActsTrackerHitsPlane";
    /// Name of the output tracker hit raw collection.
    std::string outputTrackerHitsRaw = "ActsTrackerHitsRaw";
  };

  /// Constructor with
  /// @param config configuration struct
  /// @param level logging level
  explicit EDM4hepMeasurementOutputConverter(
      const Config& config,
      std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

  std::vector<std::string> collections() const override;

 private:
  /// This implementation converts the measurements to EDM4hep.
  ///
  /// @param ctx The Algorithm context with per event information
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  Config m_cfg;

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};

  ReadDataHandle<ClusterContainer> m_inputClusters{this, "InputClusters"};

  PodioCollectionWriteHandle<edm4hep::TrackerHitPlaneCollection>
      m_outputTrackerHitsPlane{this, "OutputTrackerHitsPlane"};
  PodioCollectionWriteHandle<edm4hep::TrackerHit3DCollection>
      m_outputTrackerHitsRaw{this, "OutputTrackerHitsRaw"};
};

}  // namespace ActsExamples
