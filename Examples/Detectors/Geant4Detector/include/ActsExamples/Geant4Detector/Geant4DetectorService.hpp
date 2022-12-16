// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/ProtoDetector.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/BareService.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"

namespace Acts {
class TrackingGeometry;
namespace Experimental {
class Detector;
}  // namespace Experimental
}  // namespace Acts

namespace ActsExamples {
namespace Geant4 {

/// @class Geant4DetectorService
///
class Geant4DetectorService final : public BareService {
 public:
  /// Nested configuration struct
  struct Config {
    /// Name of the service to be built
    std::string name = "";
    /// The input file name
    std::string gdmlFile = "";
    /// The detector name
    std::string detectorName = "";
    /// The proto detector description
    Acts::ProtoDetector protoDetector;
    /// Indicate if the detector instance should be built
    bool buildDetector = true;
    /// Indicate if the tracking geometry instance should be built
    bool buildTrackingGeometry = true;
    /// Name for sensitive seleciton
    std::string sensitiveSelectionName = "sens_vol";
    /// Name for passive selection
    std::string passiveSelectionName = "pass_vol";
    /// Logging level of the master tool
    Acts::Logging::Level logLevel = Acts::Logging::INFO;
    /// Logging level of the child tools
    Acts::Logging::Level toolLogLevel = Acts::Logging::INFO;
  };

  Geant4DetectorService(const Config& cfg);
  ~Geant4DetectorService() final = default;

  /// Start-of-run hook to be called before any events are processed.
  ///
  void startRun() final;

  /// Create and return a detector from the detector construction
  ///
  /// @return a detector shared instance
  std::shared_ptr<Acts::Experimental::Detector> detector() const;

  /// Create and return a tracking eometry from the detector construction
  ///
  /// @return a detector shared instance
  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry() const;

 private:
  Config m_cfg;
  std::unique_ptr<GdmlDetectorConstruction> m_g4DetectorConstruction = nullptr;
  std::shared_ptr<Acts::Experimental::Detector> m_detector = nullptr;
  std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeometry = nullptr;
};

inline std::shared_ptr<Acts::Experimental::Detector>
Geant4DetectorService::detector() const {
  return m_detector;
}

inline std::shared_ptr<const Acts::TrackingGeometry>
Geant4DetectorService::trackingGeometry() const {
  return m_trackingGeometry;
}

}  // namespace Geant4
}  // namespace ActsExamples
