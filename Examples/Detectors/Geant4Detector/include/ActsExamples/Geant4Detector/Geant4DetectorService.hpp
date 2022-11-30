// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/BareService.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"

namespace Acts {
namespace Experimental {
struct Detector {};
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
    std::string serviceName = "";

    std::string gdmlFile = "";
    std::string detectorName = "";
  };

  Geant4DetectorService(const Config& cfg,
                        Acts::Logging::Level level = Acts::Logging::INFO);
  ~Geant4DetectorService() final = default;

  /// Start-of-run hook to be called before any events are processed.
  ///
  void startRun() final;

  /// Create and return a detector from the detector construction
  ///
  /// @return a detector shared instance
  std::shared_ptr<Acts::Experimental::Detector> detector() const;

 private:
  Config m_cfg;
  std::unique_ptr<GdmlDetectorConstruction> m_g4DetectorConstruction = nullptr;
  std::shared_ptr<Acts::Experimental::Detector> m_detector = nullptr;

};

inline std::shared_ptr<Acts::Experimental::Detector> Geant4DetectorService::detector() const {
    return m_detector;
}


}  // namespace Geant4
}  // namespace ActsExamples
