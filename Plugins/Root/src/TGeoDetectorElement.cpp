// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Root/TGeoDetectorElement.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsPlugins/Root/TGeoSurfaceConverter.hpp"

#include <cassert>
#include <utility>

#include "RtypesCore.h"
#include "TGeoBoolNode.h"

using namespace Acts;

using Line2D = Eigen::Hyperplane<double, 2>;

namespace ActsPlugins {

TGeoDetectorElement::TGeoDetectorElement(
    const Identifier& identifier, const TGeoNode& tGeoNode,
    const TGeoMatrix& tGeoMatrix, TGeoAxes axes, double scalor,
    std::shared_ptr<const ISurfaceMaterial> material)
    : m_detElement(&tGeoNode),
      m_identifier(identifier),
      m_deferredMaterial(material) {
  const Double_t* translation = tGeoMatrix.GetTranslation();
  const Double_t* rotation = tGeoMatrix.GetRotationMatrix();

  auto sensor = m_detElement->GetVolume();
  auto tgShape = sensor->GetShape();

  auto [cBounds, cTransform, cThickness] =
      TGeoSurfaceConverter::cylinderComponents(*tgShape, rotation, translation,
                                               axes, scalor);
  if (cBounds != nullptr) {
    m_transform = cTransform;
    m_bounds = cBounds;
    m_thickness = cThickness;
    // Backward-compat: create surface immediately via deprecated raw-ref path
    auto surf = Surface::makeShared<CylinderSurface>(cBounds, *this);
    surf->assignSurfaceMaterial(std::move(material));
    surf->assignThickness(m_thickness);
    m_surface = surf;
    m_legacySurface = std::move(surf);
    return;
  }

  auto [dBounds, dTransform, dThickness] = TGeoSurfaceConverter::discComponents(
      *tgShape, rotation, translation, axes, scalor);
  if (dBounds != nullptr) {
    m_bounds = dBounds;
    m_transform = dTransform;
    m_thickness = dThickness;
    auto surf = Surface::makeShared<DiscSurface>(dBounds, *this);
    surf->assignSurfaceMaterial(std::move(material));
    surf->assignThickness(m_thickness);
    m_surface = surf;
    m_legacySurface = std::move(surf);
    return;
  }

  auto [pBounds, pTransform, pThickness] =
      TGeoSurfaceConverter::planeComponents(*tgShape, rotation, translation,
                                            axes, scalor);
  if (pBounds != nullptr) {
    m_bounds = pBounds;
    m_transform = pTransform;
    m_thickness = pThickness;
    auto surf = Surface::makeShared<PlaneSurface>(pBounds, *this);
    surf->assignSurfaceMaterial(std::move(material));
    surf->assignThickness(m_thickness);
    m_surface = surf;
    m_legacySurface = std::move(surf);
  }
}

TGeoDetectorElement::TGeoDetectorElement(
    const Identifier& identifier, const TGeoNode& tGeoNode,
    const Transform3& tgTransform,
    const std::shared_ptr<const PlanarBounds>& tgBounds, double tgThickness)
    : m_detElement(&tGeoNode),
      m_transform(tgTransform),
      m_identifier(identifier),
      m_bounds(tgBounds),
      m_thickness(tgThickness) {
  auto surf = Surface::makeShared<PlaneSurface>(tgBounds, *this);
  surf->assignThickness(m_thickness);
  m_surface = surf;
  m_legacySurface = std::move(surf);
}

TGeoDetectorElement::TGeoDetectorElement(
    const Identifier& identifier, const TGeoNode& tGeoNode,
    const Transform3& tgTransform,
    const std::shared_ptr<const DiscBounds>& tgBounds, double tgThickness)
    : m_detElement(&tGeoNode),
      m_transform(tgTransform),
      m_identifier(identifier),
      m_bounds(tgBounds),
      m_thickness(tgThickness) {
  auto surf = Surface::makeShared<DiscSurface>(tgBounds, *this);
  surf->assignThickness(m_thickness);
  m_surface = surf;
  m_legacySurface = std::move(surf);
}

TGeoDetectorElement::~TGeoDetectorElement() = default;

std::shared_ptr<Acts::Surface> TGeoDetectorElement::createSurface() const {
  std::shared_ptr<Acts::Surface> surf;

  if (auto cBounds =
          std::dynamic_pointer_cast<const CylinderBounds>(m_bounds)) {
    surf = Surface::makeShared<CylinderSurface>(cBounds, shared_from_this());
  } else if (auto dBounds =
                 std::dynamic_pointer_cast<const DiscBounds>(m_bounds)) {
    surf = Surface::makeShared<DiscSurface>(dBounds, shared_from_this());
  } else if (auto pBounds =
                 std::dynamic_pointer_cast<const PlanarBounds>(m_bounds)) {
    surf = Surface::makeShared<PlaneSurface>(pBounds, shared_from_this());
  }

  if (surf != nullptr) {
    surf->assignThickness(m_thickness);
    if (m_deferredMaterial) {
      surf->assignSurfaceMaterial(m_deferredMaterial);
    }
  }

  // Replace legacy surface with the new owning surface
  m_legacySurface.reset();
  m_surface = surf;
  return surf;
}

const Acts::Surface& TGeoDetectorElement::surface() const {
  auto s = m_surface.lock();
  assert(s && "TGeoDetectorElement: call createSurface() before surface()");
  return *s;
}

Acts::Surface& TGeoDetectorElement::surface() {
  auto s = m_surface.lock();
  assert(s && "TGeoDetectorElement: call createSurface() before surface()");
  return *s;
}

const Transform3& TGeoDetectorElement::nominalTransform() const {
  return m_transform;
}

}  // namespace ActsPlugins
