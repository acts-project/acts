// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GenericDetector/GenericDetector.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DetectorCommons/Detector.hpp"
#include "ActsExamples/GenericDetector/BuildGenericDetector.hpp"

namespace ActsExamples::Generic {

GenericDetector::GenericDetector(const Config& cfg)
    : DetectorCommons::Detector(
          Acts::getDefaultLogger("GenericDetector", cfg.logLevel)),
      m_cfg(cfg) {}

void GenericDetector::buildTrackingGeometry() {
  Acts::GeometryContext gctx;
  std::vector<std::vector<std::shared_ptr<DetectorElement>>> detectorStore;
  m_trackingGeometry = ActsExamples::Generic::buildDetector<DetectorElement>(
      gctx, detectorStore, m_cfg.buildLevel, m_cfg.materialDecorator,
      m_cfg.buildProto, m_cfg.surfaceLogLevel, m_cfg.layerLogLevel,
      m_cfg.volumeLogLevel);
  for (auto& something : detectorStore) {
    for (auto& element : something) {
      m_detectorStore.push_back(std::move(element));
    }
  }
}

}  // namespace ActsExamples::Generic
