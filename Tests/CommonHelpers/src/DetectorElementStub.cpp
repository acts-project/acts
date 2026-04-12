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

#include <cassert>

using namespace Acts;

ActsTests::DetectorElementStub::DetectorElementStub(const Transform3& transform)
    : m_elementTransform(transform) {}

ActsTests::DetectorElementStub::DetectorElementStub(
    const Transform3& transform, std::shared_ptr<const CylinderBounds> cBounds,
    double thickness, std::shared_ptr<const ISurfaceMaterial> material)
    : m_elementTransform(transform),
      m_elementThickness(thickness),
      m_cylinderBounds(cBounds),
      m_deferredMaterial(material) {
  // Backward-compat: create surface immediately using deprecated raw-ref path
  // so that surface() works without calling createSurface() first.
  auto surf = Surface::makeShared<CylinderSurface>(std::move(cBounds), *this);
  surf->assignThickness(thickness);
  surf->assignSurfaceMaterial(std::move(material));
  assert(surf->surfacePlacement() == this);
  assert(surf->isSensitive() == isSensitive());
  m_elementSurface = surf;
  m_legacySurface = std::move(surf);
}

ActsTests::DetectorElementStub::DetectorElementStub(
    const Transform3& transform, std::shared_ptr<const PlanarBounds> pBounds,
    double thickness, std::shared_ptr<const ISurfaceMaterial> material)
    : m_elementTransform(transform),
      m_elementThickness(thickness),
      m_planarBounds(pBounds),
      m_deferredMaterial(material) {
  auto surf = Surface::makeShared<PlaneSurface>(std::move(pBounds), *this);
  surf->assignThickness(thickness);
  surf->assignSurfaceMaterial(std::move(material));
  assert(surf->surfacePlacement() == this);
  assert(surf->isSensitive() == isSensitive());
  m_elementSurface = surf;
  m_legacySurface = std::move(surf);
}

ActsTests::DetectorElementStub::DetectorElementStub(
    const Transform3& transform, std::shared_ptr<const LineBounds> lBounds,
    double thickness, std::shared_ptr<const ISurfaceMaterial> material)
    : m_elementTransform(transform),
      m_elementThickness(thickness),
      m_lineBounds(lBounds),
      m_deferredMaterial(material) {
  auto surf = Surface::makeShared<LineSurfaceStub>(std::move(lBounds), *this);
  surf->assignThickness(thickness);
  surf->assignSurfaceMaterial(std::move(material));
  assert(surf->surfacePlacement() == this);
  assert(surf->isSensitive() == isSensitive());
  m_elementSurface = surf;
  m_legacySurface = std::move(surf);
}

std::shared_ptr<Acts::Surface> ActsTests::DetectorElementStub::createSurface()
    const {
  std::shared_ptr<Acts::Surface> surf;

  if (m_cylinderBounds) {
    surf = Surface::makeShared<CylinderSurface>(m_cylinderBounds,
                                                shared_from_this());
  } else if (m_planarBounds) {
    surf =
        Surface::makeShared<PlaneSurface>(m_planarBounds, shared_from_this());
  } else if (m_lineBounds) {
    surf =
        Surface::makeShared<LineSurfaceStub>(m_lineBounds, shared_from_this());
  }

  if (surf != nullptr) {
    surf->assignThickness(m_elementThickness);
    surf->assignSurfaceMaterial(m_deferredMaterial);
    assert(surf->surfacePlacement() == this);
    assert(surf->isSensitive() == isSensitive());
  }

  // Replace the weak_ptr so surface() returns the new owning surface
  m_legacySurface.reset();
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
