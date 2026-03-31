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

#include <cassert>
#include <utility>

namespace ActsExamples {

GenericDetectorElement::GenericDetectorElement(
    const Identifier identifier, const Acts::Transform3& transform,
    std::shared_ptr<const Acts::PlanarBounds> pBounds, double thickness,
    std::shared_ptr<const Acts::ISurfaceMaterial> material)
    : m_elementIdentifier(identifier),
      m_elementTransform(transform),
      m_elementThickness(thickness),
      m_elementPlanarBounds(pBounds),
      m_elementDiscBounds(nullptr),
      m_elementMaterial(material) {
  // Backward-compat: create surface immediately via deprecated raw-ref path
  // so that surface() works without calling createSurface() first.
  auto surf =
      Acts::Surface::makeShared<Acts::PlaneSurface>(std::move(pBounds), *this);
  surf->assignSurfaceMaterial(std::move(material));
  surf->assignThickness(thickness);
  m_elementSurface = surf;
  m_legacySurface = std::move(surf);
}

GenericDetectorElement::GenericDetectorElement(
    const Identifier identifier, const Acts::Transform3& transform,
    std::shared_ptr<const Acts::DiscBounds> dBounds, double thickness,
    std::shared_ptr<const Acts::ISurfaceMaterial> material)
    : m_elementIdentifier(identifier),
      m_elementTransform(transform),
      m_elementThickness(thickness),
      m_elementPlanarBounds(nullptr),
      m_elementDiscBounds(dBounds),
      m_elementMaterial(material) {
  auto surf =
      Acts::Surface::makeShared<Acts::DiscSurface>(std::move(dBounds), *this);
  surf->assignSurfaceMaterial(std::move(material));
  surf->assignThickness(thickness);
  m_elementSurface = surf;
  m_legacySurface = std::move(surf);
}

std::shared_ptr<Acts::Surface> GenericDetectorElement::createSurface() const {
  std::shared_ptr<Acts::Surface> surf;
  if (m_elementPlanarBounds) {
    surf = Acts::Surface::makeShared<Acts::PlaneSurface>(m_elementPlanarBounds,
                                                         shared_from_this());
  } else {
    surf = Acts::Surface::makeShared<Acts::DiscSurface>(m_elementDiscBounds,
                                                        shared_from_this());
  }
  surf->assignThickness(m_elementThickness);
  if (m_elementMaterial) {
    surf->assignSurfaceMaterial(m_elementMaterial);
  }
  // Replace legacy surface with the new owning surface
  m_legacySurface.reset();
  m_elementSurface = surf;
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

const Acts::Surface& GenericDetectorElement::surface() const {
  auto s = m_elementSurface.lock();
  assert(s && "GenericDetectorElement: call createSurface() before surface()");
  return *s;
}

Acts::Surface& GenericDetectorElement::surface() {
  auto s = m_elementSurface.lock();
  assert(s && "GenericDetectorElement: call createSurface() before surface()");
  return *s;
}

double GenericDetectorElement::thickness() const {
  return m_elementThickness;
}

GenericDetectorElement::Identifier GenericDetectorElement::identifier() const {
  return m_elementIdentifier;
}

}  // namespace ActsExamples
