// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GenericDetector/GenericDetector.hpp"

#include "ActsExamples/GenericDetector/BuildGenericDetector.hpp"
#include "ActsExamples/GenericDetector/GenericDetectorElement.hpp"

namespace ActsExamples {

GenericDetector::GenericDetector(const Config& cfg)
    : Detector(Acts::getDefaultLogger("GenericDetector", cfg.logLevel)),
      m_cfg(cfg) {
  m_nominalGeometryContext = Acts::GeometryContext();

  auto detectorElementFactory =
      [this](std::shared_ptr<const Acts::Transform3> transform,
             std::shared_ptr<const Acts::PlanarBounds> bounds, double thickness,
             std::shared_ptr<const Acts::ISurfaceMaterial> material)
      -> std::shared_ptr<GenericDetectorElement> {
    auto id =
        static_cast<GenericDetectorElement::Identifier>(m_detectorStore.size());
    auto detElem = std::make_shared<GenericDetectorElement>(
        id, std::move(transform), std::move(bounds), thickness,
        std::move(material));
    m_detectorStore.push_back(detElem);
    return detElem;
  };

  m_trackingGeometry = Generic::buildDetector(
      m_nominalGeometryContext, detectorElementFactory, m_cfg.buildLevel,
      m_cfg.materialDecorator, m_cfg.buildProto, m_cfg.surfaceLogLevel,
      m_cfg.layerLogLevel, m_cfg.volumeLogLevel, m_cfg.gen3);
}

}  // namespace ActsExamples
