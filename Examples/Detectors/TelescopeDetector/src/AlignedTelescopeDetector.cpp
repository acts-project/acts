// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TelescopeDetector/AlignedTelescopeDetector.hpp"

#include "ActsExamples/TelescopeDetector/BuildTelescopeDetector.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetectorElement.hpp"

namespace ActsExamples {

AlignedTelescopeDetector::AlignedTelescopeDetector(const Config& cfg)
    : TelescopeDetector(cfg, NoBuildTag{}) {
  m_nominalGeometryContext =
      Acts::GeometryContext::dangerouslyDefaultConstruct();

  // Set the detector element factory
  auto alignedDetectorElementFactory =
      [](const Acts::Transform3& transform,
         std::variant<std::shared_ptr<const Acts::PlanarBounds>,
                      std::shared_ptr<const Acts::DiscBounds>>
             bounds,
         double thickness,
         std::shared_ptr<const Acts::ISurfaceMaterial> material,
         std::vector<std::shared_ptr<const Acts::SurfacePlacementBase>>&
             detStore) {
        auto ID =
            static_cast<TelescopeDetectorElement::Identifier>(detStore.size());
        std::shared_ptr<AlignedTelescopeDetectorElement> detElem;
        if (bounds.index() == 0) {
          detElem = std::make_shared<AlignedTelescopeDetectorElement>(
              ID, std::make_shared<Acts::Transform3>(transform),
              std::get<std::shared_ptr<const Acts::PlanarBounds>>(bounds),
              thickness, std::move(material));
        } else {
          detElem = std::make_shared<AlignedTelescopeDetectorElement>(
              ID, std::make_shared<Acts::Transform3>(transform),
              std::get<std::shared_ptr<const Acts::DiscBounds>>(bounds),
              thickness, std::move(material));
        }
        detStore.push_back(detElem);
        return detElem;
      };
  m_trackingGeometry = buildTelescopeDetector(
      m_nominalGeometryContext, alignedDetectorElementFactory, m_detectorStore,
      cfg.positions, cfg.stereos, cfg.offsets, cfg.bounds, cfg.thickness,
      static_cast<TelescopeSurfaceType>(cfg.surfaceType),
      static_cast<Acts::AxisDirection>(cfg.binValue));
}

}  // namespace ActsExamples
