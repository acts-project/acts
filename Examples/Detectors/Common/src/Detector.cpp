// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/DetectorCommons/Detector.hpp"

#include "Acts/Utilities/Logger.hpp"

namespace ActsExamples::DetectorCommons {

Detector::Detector(std::unique_ptr<const Acts::Logger> logger)
    : m_logger(std::move(logger)) {}

std::tuple<Detector::TrackingGeometryPtr, Detector::ContextDecorators,
           Detector::DetectorStore>
Detector::trackingGeometry() {
  if (m_trackingGeometry == nullptr) {
    buildTrackingGeometry();
  }
  return {m_trackingGeometry, m_contextDecorators, m_detectorStore};
}

std::tuple<Detector::DetectorPtr, Detector::ContextDecorators,
           Detector::DetectorStore>
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

void Detector::buildDetector() {}

}  // namespace ActsExamples::DetectorCommons
