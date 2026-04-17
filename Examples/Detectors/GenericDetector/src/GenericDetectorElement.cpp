// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GenericDetector/GenericDetectorElement.hpp"

#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

#include <utility>

namespace ActsExamples {

GenericDetectorElement::GenericDetectorElement(
    const Identifier identifier, const Acts::Transform3& transform,
    std::shared_ptr<const Acts::PlanarBounds> pBounds, double thickness,
    std::shared_ptr<const Acts::ISurfaceMaterial> material)
    : m_elementIdentifier(identifier),
      m_elementTransform(transform),
      m_elementThickness(thickness),
      m_elementPlanarBounds(std::move(pBounds)),
      m_elementDiscBounds(nullptr),
      m_elementMaterial(std::move(material)) {}

GenericDetectorElement::GenericDetectorElement(
    const Identifier identifier, const Acts::Transform3& transform,
    std::shared_ptr<const Acts::DiscBounds> dBounds, double thickness,
    std::shared_ptr<const Acts::ISurfaceMaterial> material)
    : m_elementIdentifier(identifier),
      m_elementTransform(transform),
      m_elementThickness(thickness),
      m_elementPlanarBounds(nullptr),
      m_elementDiscBounds(std::move(dBounds)),
      m_elementMaterial(std::move(material)) {}

std::shared_ptr<Acts::Surface> GenericDetectorElement::createSurface() {
  std::shared_ptr<Acts::Surface> surf;
  if (m_elementPlanarBounds) {
    surf = Acts::Surface::makeShared<Acts::PlaneSurface>(m_elementTransform,
                                                         m_elementPlanarBounds);
  } else {
    surf = Acts::Surface::makeShared<Acts::DiscSurface>(m_elementTransform,
                                                        m_elementDiscBounds);
  }
  assignSurface(surf);
  surf->assignThickness(m_elementThickness);
  if (m_elementMaterial) {
    surf->assignSurfaceMaterial(m_elementMaterial);
  }
  return surf;
}

std::pair<std::shared_ptr<GenericDetectorElement>,
          std::shared_ptr<Acts::Surface>>
GenericDetectorElement::create(
    const Identifier identifier, const Acts::Transform3& transform,
    std::shared_ptr<const Acts::PlanarBounds> pBounds, double thickness,
    std::shared_ptr<const Acts::ISurfaceMaterial> material) {
  auto el = std::make_shared<GenericDetectorElement>(
      identifier, transform, std::move(pBounds), thickness,
      std::move(material));
  auto surf = el->createSurface();
  return {std::move(el), std::move(surf)};
}

std::pair<std::shared_ptr<GenericDetectorElement>,
          std::shared_ptr<Acts::Surface>>
GenericDetectorElement::create(
    const Identifier identifier, const Acts::Transform3& transform,
    std::shared_ptr<const Acts::DiscBounds> dBounds, double thickness,
    std::shared_ptr<const Acts::ISurfaceMaterial> material) {
  auto el = std::make_shared<GenericDetectorElement>(
      identifier, transform, std::move(dBounds), thickness,
      std::move(material));
  auto surf = el->createSurface();
  return {std::move(el), std::move(surf)};
}

const Acts::Transform3& GenericDetectorElement::localToGlobalTransform(
    const Acts::GeometryContext& /*gctx*/) const {
  return m_elementTransform;
}

const Acts::Transform3& GenericDetectorElement::nominalTransform() const {
  return m_elementTransform;
}

double GenericDetectorElement::thickness() const {
  return m_elementThickness;
}

GenericDetectorElement::Identifier GenericDetectorElement::identifier() const {
  return m_elementIdentifier;
}

}  // namespace ActsExamples
