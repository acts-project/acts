// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4Detector/Geant4Detector.hpp"

#include "Acts/Geometry/Detector.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"

#include <boost/program_options.hpp>

auto ActsExamples::Geant4::Geant4Detector::constructDetector(
    const ActsExamples::Geant4::Geant4DetectorService::Config& cfg)
    -> std::pair<DetectorPtr, ContextDecorators> {
  auto g4DetectorService = ActsExamples::Geant4::Geant4DetectorService(cfg);

  g4DetectorService.startRun();

  return std::make_pair<DetectorPtr, ContextDecorators>(nullptr, {});
}

auto ActsExamples::Geant4::Geant4Detector::constructTrackingGeometry(
    const ActsExamples::Geant4::Geant4DetectorService::Config& cfg)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  auto g4DetectorService = ActsExamples::Geant4::Geant4DetectorService(cfg);

  g4DetectorService.startRun();

  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(
      g4DetectorService.trackingGeometry(), {});
}