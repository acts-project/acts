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
    : DetectorElementBase(), m_elementTransform(transform) {}

ActsTests::DetectorElementStub::DetectorElementStub(
    const Transform3& transform, std::shared_ptr<const CylinderBounds> cBounds,
    double thickness, std::shared_ptr<const ISurfaceMaterial> material)
    : DetectorElementBase(),
      m_elementTransform(transform),
      m_elementThickness(thickness) {
  m_elementSurface =
      Surface::makeShared<CylinderSurface>(std::move(cBounds), *this);
  m_elementSurface->assignSurfaceMaterial(std::move(material));
}

ActsTests::DetectorElementStub::DetectorElementStub(
    const Transform3& transform, std::shared_ptr<const PlanarBounds> pBounds,
    double thickness, std::shared_ptr<const ISurfaceMaterial> material)
    : DetectorElementBase(),
      m_elementTransform(transform),
      m_elementThickness(thickness) {
  m_elementSurface =
      Surface::makeShared<PlaneSurface>(std::move(pBounds), *this);
  m_elementSurface->assignSurfaceMaterial(std::move(material));
}

ActsTests::DetectorElementStub::DetectorElementStub(
    const Transform3& transform, std::shared_ptr<const LineBounds> lBounds,
    double thickness, std::shared_ptr<const ISurfaceMaterial> material)
    : DetectorElementBase(),
      m_elementTransform(transform),
      m_elementThickness(thickness) {
  m_elementSurface =
      Surface::makeShared<LineSurfaceStub>(std::move(lBounds), *this);
  m_elementSurface->assignSurfaceMaterial(std::move(material));
}
