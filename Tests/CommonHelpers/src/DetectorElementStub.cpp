// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsTests/CommonHelpers/DetectorElementStub.hpp"

#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/SurfacePlacementBase.hpp"
#include "ActsTests/CommonHelpers/LineSurfaceStub.hpp"

using namespace Acts;

ActsTests::DetectorElementStub::DetectorElementStub(const Transform3& transform)
    : m_elementTransform(transform) {}

ActsTests::DetectorElementStub::DetectorElementStub(
    const Transform3& transform, std::shared_ptr<const CylinderBounds> cBounds,
    double thickness, std::shared_ptr<const ISurfaceMaterial> material)
    : m_elementTransform(transform),
      m_elementThickness(thickness),
      m_cylinderBounds(std::move(cBounds)),
      m_deferredMaterial(std::move(material)) {
  auto surf = Surface::makeShared<CylinderSurface>(transform, m_cylinderBounds);
  surf->assignThickness(m_elementThickness);
  surf->assignSurfaceMaterial(m_deferredMaterial);
  m_elementSurface = surf;
  m_legacySurface = std::move(surf);
}

ActsTests::DetectorElementStub::DetectorElementStub(
    const Transform3& transform, std::shared_ptr<const PlanarBounds> pBounds,
    double thickness, std::shared_ptr<const ISurfaceMaterial> material)
    : m_elementTransform(transform),
      m_elementThickness(thickness),
      m_planarBounds(std::move(pBounds)),
      m_deferredMaterial(std::move(material)) {
  auto surf = Surface::makeShared<PlaneSurface>(transform, m_planarBounds);
  surf->assignThickness(m_elementThickness);
  surf->assignSurfaceMaterial(m_deferredMaterial);
  m_elementSurface = surf;
  m_legacySurface = std::move(surf);
}

ActsTests::DetectorElementStub::DetectorElementStub(
    const Transform3& transform, std::shared_ptr<const LineBounds> lBounds,
    double thickness, std::shared_ptr<const ISurfaceMaterial> material)
    : m_elementTransform(transform),
      m_elementThickness(thickness),
      m_lineBounds(std::move(lBounds)),
      m_deferredMaterial(std::move(material)) {
  auto surf = Surface::makeShared<LineSurfaceStub>(transform, m_lineBounds);
  surf->assignThickness(m_elementThickness);
  surf->assignSurfaceMaterial(m_deferredMaterial);
  m_elementSurface = surf;
  m_legacySurface = std::move(surf);
}

std::shared_ptr<Acts::Surface> ActsTests::DetectorElementStub::createSurface() {
  std::shared_ptr<Acts::Surface> surf;

  if (m_cylinderBounds) {
    surf = Surface::makeShared<CylinderSurface>(m_elementTransform,
                                                m_cylinderBounds);
  } else if (m_planarBounds) {
    surf = Surface::makeShared<PlaneSurface>(m_elementTransform, m_planarBounds);
  } else if (m_lineBounds) {
    surf =
        Surface::makeShared<LineSurfaceStub>(m_elementTransform, m_lineBounds);
  }

  if (surf != nullptr) {
    assignSurface(surf);
    surf->assignThickness(m_elementThickness);
    surf->assignSurfaceMaterial(m_deferredMaterial);
  }

  m_legacySurface = surf;
  m_elementSurface = surf;
  return surf;
}

const Acts::Surface& ActsTests::DetectorElementStub::surface() const {
  auto s = m_elementSurface.lock();
  assert(s && "DetectorElementStub: call createSurface() before surface()");
  return *s;
}

Acts::Surface& ActsTests::DetectorElementStub::surface() {
  auto s = m_elementSurface.lock();
  assert(s && "DetectorElementStub: call createSurface() before surface()");
  return *s;
}
