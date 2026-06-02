// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GenericDetector/AlignedGenericDetector.hpp"

namespace ActsExamples {

AlignedGenericDetector::AlignedGenericDetector(const Config& cfg)
    : GenericDetector(cfg, GenericDetector::NoBuildTag{}) {
  m_nominalGeometryContext =
      Acts::GeometryContext::dangerouslyDefaultConstruct();

  // Set the detector element factory
  auto alignedDetectorElementFactory =
      [&](const Acts::Transform3& transform,
          std::shared_ptr<const Acts::PlanarBounds> bounds, double thickness,
          std::shared_ptr<const Acts::ISurfaceMaterial> material) {
        auto id = static_cast<GenericDetectorElement::Identifier>(
            m_detectorStore.size());
        auto detElem = std::make_shared<AlignedGenericDetectorElement>(
            id, transform, std::move(bounds), thickness, std::move(material));
        m_detectorStore.push_back(detElem);
        return detElem;
      };
  buildTrackingGeometry(alignedDetectorElementFactory);
}

}  // namespace ActsExamples
