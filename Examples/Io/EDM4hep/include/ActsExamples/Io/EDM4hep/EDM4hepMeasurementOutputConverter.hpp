// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/Podio/PodioCollectionDataHandle.hpp"
#include "ActsExamples/Io/Podio/PodioOutputConverter.hpp"
#include "ActsPlugins/EDM4hep/EDM4hepUtil.hpp"
#include "ActsPodioEdm/TrackerHitLocalCollection.h"
#include "ActsPodioEdm/TrackerHitLocalSimTrackerHitLinkCollection.h"

#include <memory>
#include <string>

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
    /// Name of the output tracker hit raw collection.
    std::string outputTrackerHitsLocal;

    /// Optional sim hit linking. All three fields must be set together.
    /// Input sim hit association (internal index ↔ edm4hep hit).
    std::optional<std::string> inputSimHitAssociation;
    /// Input map from measurement index to internal sim hit index.
    std::optional<std::string> inputMeasurementSimHitsMap;
    /// Name of the output TrackerHitLocalSimTrackerHitLink collection.
    std::optional<std::string> outputSimHitLinks;

    /// Tracking geometry for surface lookup (local-to-global transform).
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
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

  ReadDataHandle<ActsPlugins::EDM4hepUtil::SimHitAssociation>
      m_inputSimHitAssociation{this, "InputSimHitAssociation"};

  ReadDataHandle<IndexMultimap<Index>> m_inputMeasurementSimHitsMap{
      this, "InputMeasurementSimHitsMap"};

  PodioCollectionWriteHandle<ActsPodioEdm::TrackerHitLocalCollection>
      m_outputTrackerHitsLocal{this, "OutputTrackerHitsLocal"};

  PodioCollectionWriteHandle<
      ActsPodioEdm::TrackerHitLocalSimTrackerHitLinkCollection>
      m_outputSimHitLinks{this, "OutputSimHitLinks"};
};

}  // namespace ActsExamples
