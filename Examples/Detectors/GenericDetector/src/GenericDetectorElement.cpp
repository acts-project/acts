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
    : Acts::DetectorElementBase(),
      m_elementIdentifier(identifier),
      m_elementTransform(transform),
      m_elementSurface(
          Acts::Surface::makeShared<Acts::PlaneSurface>(pBounds, *this)),
      m_elementThickness(thickness),
      m_elementPlanarBounds(std::move(pBounds)),
      m_elementDiscBounds(nullptr) {
  m_elementSurface->assignSurfaceMaterial(std::move(material));
}

GenericDetectorElement::GenericDetectorElement(
    const Identifier identifier, const Acts::Transform3& transform,
    std::shared_ptr<const Acts::DiscBounds> dBounds, double thickness,
    std::shared_ptr<const Acts::ISurfaceMaterial> material)
    : Acts::DetectorElementBase(),
      m_elementIdentifier(identifier),
      m_elementTransform(transform),
      m_elementSurface(
          Acts::Surface::makeShared<Acts::DiscSurface>(dBounds, *this)),
      m_elementThickness(thickness),
      m_elementPlanarBounds(nullptr),
      m_elementDiscBounds(std::move(dBounds)) {
  m_elementSurface->assignSurfaceMaterial(std::move(material));
}

const Acts::Transform3& GenericDetectorElement::transform(
    const Acts::GeometryContext& /*gctx*/) const {
  return m_elementTransform;
}

const Acts::Transform3& GenericDetectorElement::nominalTransform() const {
  return m_elementTransform;
}

const Acts::Surface& GenericDetectorElement::surface() const {
  return *m_elementSurface;
}

Acts::Surface& GenericDetectorElement::surface() {
  return *m_elementSurface;
}

double GenericDetectorElement::thickness() const {
  return m_elementThickness;
}

GenericDetectorElement::Identifier GenericDetectorElement::identifier() const {
  return m_elementIdentifier;
}

}  // namespace ActsExamples
