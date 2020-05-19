// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderVolumeBounds.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Visualization/IVisualization.hpp"

#include <cmath>
#include <iostream>

Acts::CylinderVolumeBounds::CylinderVolumeBounds(
    const CylinderBounds& cBounds, double thickness) noexcept(false)
    : VolumeBounds() {
  double cR = cBounds.get(CylinderBounds::eR);
  if (thickness <= 0. or (cR - 0.5 * thickness) < 0.) {
    throw(std::invalid_argument(
        "CylinderVolumeBounds: invalid extrusion thickness."));
  }
  m_values[eMinR] = cR - 0.5 * thickness;
  m_values[eMaxR] = cR + 0.5 * thickness;
  m_values[eHalfLengthZ] = cBounds.get(CylinderBounds::eHalfLengthZ);
  m_values[eHalfPhiSector] = cBounds.get(CylinderBounds::eHalfPhiSector);
  m_values[eAveragePhi] = cBounds.get(CylinderBounds::eAveragePhi);
  buildSurfaceBounds();
}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(
    const RadialBounds& rBounds, double thickness) noexcept(false)
    : VolumeBounds() {
  if (thickness <= 0.) {
    throw(std::invalid_argument(
        "CylinderVolumeBounds: invalid extrusion thickness."));
  }
  m_values[eMinR] = rBounds.get(RadialBounds::eMinR);
  m_values[eMaxR] = rBounds.get(RadialBounds::eMaxR);
  m_values[eHalfLengthZ] = 0.5 * thickness;
  m_values[eHalfPhiSector] = rBounds.get(RadialBounds::eHalfPhiSector);
  m_values[eAveragePhi] = rBounds.get(RadialBounds::eAveragePhi);
  buildSurfaceBounds();
}

Acts::OrientedSurfaces Acts::CylinderVolumeBounds::orientedSurfaces(
    const Transform3D* transformPtr) const {
  OrientedSurfaces oSurfaces;
  oSurfaces.reserve(6);

  // Set the transform
  Transform3D transform =
      (transformPtr == nullptr) ? Transform3D::Identity() : (*transformPtr);
  auto trfShared = std::make_shared<Transform3D>(transform);

  // [0] Bottom Disc (negative z)
  auto dSurface = Surface::makeShared<DiscSurface>(
      std::make_shared<const Transform3D>(
          transform * Translation3D(0., 0., -get(eHalfLengthZ))),
      m_discBounds);
  oSurfaces.push_back(OrientedSurface(std::move(dSurface), forward));
  // [1] Top Disc (positive z)
  dSurface = Surface::makeShared<DiscSurface>(
      std::make_shared<const Transform3D>(
          transform * Translation3D(0., 0., get(eHalfLengthZ))),
      m_discBounds);
  oSurfaces.push_back(OrientedSurface(std::move(dSurface), backward));

  // [2] Outer Cylinder
  auto cSurface =
      Surface::makeShared<CylinderSurface>(trfShared, m_outerCylinderBounds);
  oSurfaces.push_back(OrientedSurface(std::move(cSurface), backward));

  // [3] Inner Cylinder (optional)
  if (m_innerCylinderBounds != nullptr) {
    cSurface =
        Surface::makeShared<CylinderSurface>(trfShared, m_innerCylinderBounds);
    oSurfaces.push_back(OrientedSurface(std::move(cSurface), forward));
  }

  // [4] & [5] - Sectoral planes (optional)
  if (m_sectorPlaneBounds != nullptr) {
    // sectorPlane 1 (negative phi)
    const Transform3D* sp1Transform = new Transform3D(
        transform * AngleAxis3D(-get(eHalfPhiSector), Vector3D(0., 0., 1.)) *
        Translation3D(0.5 * (get(eMinR) + get(eMaxR)), 0., 0.) *
        AngleAxis3D(M_PI / 2, Vector3D(1., 0., 0.)));
    auto pSurface = Surface::makeShared<PlaneSurface>(
        std::shared_ptr<const Transform3D>(sp1Transform), m_sectorPlaneBounds);
    oSurfaces.push_back(OrientedSurface(std::move(pSurface), forward));
    // sectorPlane 2 (positive phi)
    const Transform3D* sp2Transform = new Transform3D(
        transform * AngleAxis3D(get(eHalfPhiSector), Vector3D(0., 0., 1.)) *
        Translation3D(0.5 * (get(eMinR) + get(eMaxR)), 0., 0.) *
        AngleAxis3D(-M_PI / 2, Vector3D(1., 0., 0.)));
    pSurface = Surface::makeShared<PlaneSurface>(
        std::shared_ptr<const Transform3D>(sp2Transform), m_sectorPlaneBounds);
    oSurfaces.push_back(OrientedSurface(std::move(pSurface), backward));
  }
  return oSurfaces;
}

void Acts::CylinderVolumeBounds::buildSurfaceBounds() {
  if (get(eMinR) > s_epsilon) {
    m_innerCylinderBounds = std::make_shared<const CylinderBounds>(
        get(eMinR), get(eHalfLengthZ), get(eHalfPhiSector), get(eAveragePhi));
  }
  m_outerCylinderBounds = std::make_shared<const CylinderBounds>(
      get(eMaxR), get(eHalfLengthZ), get(eHalfPhiSector), get(eAveragePhi));
  m_discBounds = std::make_shared<const RadialBounds>(
      get(eMinR), get(eMaxR), get(eHalfPhiSector), get(eAveragePhi));

  if (std::abs(get(eHalfPhiSector) - M_PI) > s_epsilon) {
    m_sectorPlaneBounds = std::make_shared<const RectangleBounds>(
        0.5 * (get(eMaxR) - get(eMinR)), get(eHalfLengthZ));
  }
}

std::ostream& Acts::CylinderVolumeBounds::toStream(std::ostream& sl) const {
  return dumpT<std::ostream>(sl);
}

Acts::Volume::BoundingBox Acts::CylinderVolumeBounds::boundingBox(
    const Transform3D* trf, const Vector3D& envelope,
    const Volume* entity) const {
  double xmax, xmin, ymax, ymin;
  xmax = get(eMaxR);

  if (get(eHalfPhiSector) > M_PI / 2.) {
    // more than half
    ymax = xmax;
    ymin = -xmax;
    xmin = xmax * std::cos(get(eHalfPhiSector));
  } else {
    // less than half
    ymax = get(eMaxR) * std::sin(get(eHalfPhiSector));
    ymin = -ymax;
    // in this case, xmin is given by the inner radius
    xmin = get(eMinR) * std::cos(get(eHalfPhiSector));
  }

  Vector3D vmin(xmin, ymin, -get(eHalfLengthZ));
  Vector3D vmax(xmax, ymax, get(eHalfLengthZ));

  // this is probably not perfect, but at least conservative
  Volume::BoundingBox box{entity, vmin - envelope, vmax + envelope};
  return trf == nullptr ? box : box.transformed(*trf);
}
