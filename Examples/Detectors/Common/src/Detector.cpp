// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DetectorCommons/Detector.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace ActsExamples::DetectorCommons {

Detector::Detector(std::unique_ptr<const Acts::Logger> logger)
    : m_logger(std::move(logger)) {}

Detector::Detector(Detector&&) = default;

Detector::~Detector() = default;

Detector& Detector::operator=(Detector&&) = default;

std::tuple<std::shared_ptr<const Acts::TrackingGeometry>,
           Detector::ContextDecorators, Detector::DetectorStore>
Detector::trackingGeometry() {
  if (m_trackingGeometry == nullptr) {
    Acts::GeometryContext gctx = buildGeometryContext();
    buildTrackingGeometry(gctx);
  }
  return {m_trackingGeometry, m_contextDecorators, m_detectorStore};
}

std::tuple<std::shared_ptr<Acts::Experimental::Detector>,
           Detector::ContextDecorators, Detector::DetectorStore>
Detector::detector() {
  if (m_detector == nullptr) {
    buildDetector();
  }
  return {m_detector, m_contextDecorators, m_detectorStore};
}

void Detector::drop() {
  m_trackingGeometry.reset();
  m_detector.reset();
  m_contextDecorators.clear();
  m_detectorStore.clear();
}

Acts::GeometryContext Detector::buildGeometryContext() const {
  return Acts::GeometryContext();
}

void Detector::buildDetector() {}

}  // namespace ActsExamples::DetectorCommons
