// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TelescopeDetector/TelescopeDetectorElement.hpp"

#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

namespace ActsExamples {

TelescopeDetectorElement::TelescopeDetectorElement(
    std::shared_ptr<const Acts::Transform3> transform,
    std::shared_ptr<const Acts::PlanarBounds> pBounds, double thickness,
    std::shared_ptr<const Acts::ISurfaceMaterial> material)
    : m_elementTransform(std::move(transform)),
      m_elementThickness(thickness),
      m_elementPlanarBounds(std::move(pBounds)),
      m_elementDiscBounds(nullptr),
      m_elementMaterial(std::move(material)) {}

TelescopeDetectorElement::TelescopeDetectorElement(
    std::shared_ptr<const Acts::Transform3> transform,
    std::shared_ptr<const Acts::DiscBounds> dBounds, double thickness,
    std::shared_ptr<const Acts::ISurfaceMaterial> material)
    : m_elementTransform(std::move(transform)),
      m_elementThickness(thickness),
      m_elementPlanarBounds(nullptr),
      m_elementDiscBounds(std::move(dBounds)),
      m_elementMaterial(std::move(material)) {}

std::shared_ptr<Acts::Surface> TelescopeDetectorElement::createSurface() {
  std::shared_ptr<Acts::Surface> surf;
  if (m_elementPlanarBounds) {
    surf = Acts::Surface::makeShared<Acts::PlaneSurface>(*m_elementTransform,
                                                         m_elementPlanarBounds);
  } else {
    surf = Acts::Surface::makeShared<Acts::DiscSurface>(*m_elementTransform,
                                                        m_elementDiscBounds);
  }
  assignSurface(surf);
  surf->assignThickness(m_elementThickness);
  if (m_elementMaterial) {
    surf->assignSurfaceMaterial(m_elementMaterial);
  }
  return surf;
}

}  // namespace ActsExamples
