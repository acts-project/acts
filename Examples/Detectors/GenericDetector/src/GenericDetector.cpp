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

GenericDetectorFactory::GenericDetectorFactory(const Config& cfg)
    : DetectorFactoryBase(
          Acts::getDefaultLogger("GenericDetectorFactory", cfg.logLevel)),
      m_cfg(cfg) {}

std::shared_ptr<DetectorBase> GenericDetectorFactory::buildDetector() const {
  Acts::GeometryContext geometryContext;
  std::vector<std::shared_ptr<const Acts::DetectorElementBase>> detectorStore;
  std::shared_ptr<const Acts::TrackingGeometry> gen1Geometry;
  std::shared_ptr<Acts::Experimental::Detector> gen2Geometry;
  std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>
      contextDecorators;

  geometryContext = Acts::GeometryContext();

  std::vector<std::vector<std::shared_ptr<GenericDetectorElement>>>
      specificDetectorStore;
  gen1Geometry = Generic::buildDetector<GenericDetectorElement>(
      geometryContext, specificDetectorStore, m_cfg.buildLevel,
      m_cfg.materialDecorator, m_cfg.buildProto, m_cfg.surfaceLogLevel,
      m_cfg.layerLogLevel, m_cfg.volumeLogLevel);

  for (const auto& something : specificDetectorStore) {
    for (const auto& element : something) {
      detectorStore.push_back(element);
    }
  }

  return std::make_shared<PreConstructedDetector>(
      std::move(geometryContext), std::move(detectorStore),
      std::move(gen1Geometry), std::move(gen2Geometry),
      std::move(contextDecorators));
}

}  // namespace ActsExamples
