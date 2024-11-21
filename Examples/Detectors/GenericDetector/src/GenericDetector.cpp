// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GenericDetector/GenericDetector.hpp"

#include "ActsExamples/DetectorCommons/DetectorBase.hpp"
#include "ActsExamples/GenericDetector/BuildGenericDetector.hpp"
#include "ActsExamples/GenericDetector/GenericDetectorElement.hpp"

namespace ActsExamples {

GenericDetector::GenericDetector(const Config& cfg)
    : DetectorBase(Acts::getDefaultLogger("GenericDetector", cfg.logLevel)),
      m_cfg(cfg) {}

Gen1GeometryHolder GenericDetector::buildGen1Geometry() {
  Gen1GeometryHolder result;

  std::vector<std::vector<std::shared_ptr<GenericDetectorElement>>>
      detectorStore;
  result.trackingGeometry =
      ActsExamples::Generic::buildDetector<GenericDetectorElement>(
          result.geometryContext, detectorStore, m_cfg.buildLevel,
          m_cfg.materialDecorator, m_cfg.buildProto, m_cfg.surfaceLogLevel,
          m_cfg.layerLogLevel, m_cfg.volumeLogLevel);

  for (const auto& something : detectorStore) {
    for (const auto& element : something) {
      result.detectorStore.push_back(element);
    }
  }

  return result;
}

}  // namespace ActsExamples
