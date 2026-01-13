// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DetectorCommons/Detector.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"

namespace ActsExamples {

Detector::Detector(std::unique_ptr<const Acts::Logger> logger)
    : m_logger(std::move(logger)),
      m_nominalGeometryContext(
          Acts::GeometryContext::dangerouslyDefaultConstruct()) {}

Detector::~Detector() = default;

std::vector<std::shared_ptr<IContextDecorator>> Detector::contextDecorators()
    const {
  return m_contextDecorators;
}

std::unique_ptr<G4VUserDetectorConstruction>
Detector::buildGeant4DetectorConstruction(
    const Geant4ConstructionOptions& /*options*/) const {
  throw std::runtime_error("Geant4 detector construction is not available.");
}

const Acts::GeometryContext& Detector::nominalGeometryContext() const {
  return m_nominalGeometryContext;
}

std::shared_ptr<const Acts::TrackingGeometry> Detector::trackingGeometry()
    const {
  if (m_trackingGeometry == nullptr) {
    throw std::runtime_error("Tracking geometry is not built");
  }
  return m_trackingGeometry;
}

const Acts::Logger& Detector::logger() const {
  return *m_logger;
}

}  // namespace ActsExamples
