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
#include <cmath>
#include <iostream>
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/IVisualization.hpp"

const double Acts::CylinderVolumeBounds::s_numericalStable = 10e-2;

Acts::CylinderVolumeBounds::CylinderVolumeBounds()
    : VolumeBounds(), m_valueStore(4, 0.) {}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(double radius, double halez)
    : VolumeBounds(), m_valueStore(4, 0.) {
  m_valueStore.at(bv_innerRadius) = 0.;
  m_valueStore.at(bv_outerRadius) = std::abs(radius);
  m_valueStore.at(bv_halfPhiSector) = M_PI;
  m_valueStore.at(bv_halfZ) = std::abs(halez);
}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(double rinner, double router,
                                                 double halez)
    : VolumeBounds(), m_valueStore(4, 0.) {
  m_valueStore.at(bv_innerRadius) = std::abs(rinner);
  m_valueStore.at(bv_outerRadius) = std::abs(router);
  m_valueStore.at(bv_halfPhiSector) = M_PI;
  m_valueStore.at(bv_halfZ) = std::abs(halez);
}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(double rinner, double router,
                                                 double haphi, double halez)
    : VolumeBounds(), m_valueStore(4, 0.) {
  m_valueStore.at(bv_innerRadius) = std::abs(rinner);
  m_valueStore.at(bv_outerRadius) = std::abs(router);
  m_valueStore.at(bv_halfPhiSector) = std::abs(haphi);
  m_valueStore.at(bv_halfZ) = std::abs(halez);
}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(const CylinderBounds& cBounds,
                                                 double thickness)
    : VolumeBounds(), m_valueStore(4, 0.) {
  double cR = cBounds.r();
  m_valueStore.at(bv_innerRadius) = cR - 0.5 * thickness;
  m_valueStore.at(bv_outerRadius) = cR + 0.5 * thickness;
  m_valueStore.at(bv_halfPhiSector) = cBounds.halfPhiSector();
  m_valueStore.at(bv_halfZ) = cBounds.halflengthZ();
}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(const RadialBounds& rBounds,
                                                 double thickness)
    : VolumeBounds(), m_valueStore(4, 0.) {
  m_valueStore.at(bv_innerRadius) = rBounds.rMin();
  m_valueStore.at(bv_outerRadius) = rBounds.rMax();
  m_valueStore.at(bv_halfPhiSector) = rBounds.halfPhiSector();
  m_valueStore.at(bv_halfZ) = 0.5 * thickness;
}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(
    const CylinderVolumeBounds& cylbo)
    : VolumeBounds(), m_valueStore(cylbo.m_valueStore) {}

Acts::CylinderVolumeBounds::~CylinderVolumeBounds() = default;

Acts::CylinderVolumeBounds& Acts::CylinderVolumeBounds::operator=(
    const CylinderVolumeBounds& cylbo) {
  if (this != &cylbo) {
    m_valueStore = cylbo.m_valueStore;
  }
  return *this;
}

std::vector<std::shared_ptr<const Acts::Surface>>
Acts::CylinderVolumeBounds::decomposeToSurfaces(
    const Transform3D* transformPtr) const {
  std::vector<std::shared_ptr<const Surface>> rSurfaces;
  rSurfaces.reserve(6);

  // set the transform
  Transform3D transform =
      (transformPtr == nullptr) ? Transform3D::Identity() : (*transformPtr);
  auto trfShared = std::make_shared<Transform3D>(transform);

  const Transform3D* tTransform = nullptr;
  Vector3D cylCenter(transform.translation());

  std::shared_ptr<const DiscBounds> dBounds = discBounds();
  // bottom Disc (negative z)
  tTransform =
      new Transform3D(transform * AngleAxis3D(M_PI, Vector3D(1., 0., 0.)) *
                      Translation3D(Vector3D(0., 0., halflengthZ())));
  rSurfaces.push_back(Surface::makeShared<DiscSurface>(
      std::shared_ptr<const Transform3D>(tTransform), dBounds));
  // top Disc (positive z)
  tTransform = new Transform3D(transform *
                               Translation3D(Vector3D(0., 0., halflengthZ())));
  rSurfaces.push_back(Surface::makeShared<DiscSurface>(
      std::shared_ptr<const Transform3D>(tTransform), dBounds));

  // outer Cylinder - shares the transform
  rSurfaces.push_back(
      Surface::makeShared<CylinderSurface>(trfShared, outerCylinderBounds()));

  // innermost Cylinder
  if (innerRadius() > s_numericalStable) {
    rSurfaces.push_back(
        Surface::makeShared<CylinderSurface>(trfShared, innerCylinderBounds()));
  }

  // the cylinder is sectoral
  if (std::abs(halfPhiSector() - M_PI) > s_numericalStable) {
    std::shared_ptr<const PlanarBounds> sp12Bounds = sectorPlaneBounds();
    // sectorPlane 1 (negative phi)
    const Transform3D* sp1Transform = new Transform3D(
        transform * AngleAxis3D(-halfPhiSector(), Vector3D(0., 0., 1.)) *
        Translation3D(Vector3D(mediumRadius(), 0., 0.)) *
        AngleAxis3D(M_PI / 2, Vector3D(1., 0., 0.)));
    rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
        std::shared_ptr<const Transform3D>(sp1Transform), sp12Bounds));
    // sectorPlane 2 (positive phi)
    const Transform3D* sp2Transform = new Transform3D(
        transform * AngleAxis3D(halfPhiSector(), Vector3D(0., 0., 1.)) *
        Translation3D(Vector3D(mediumRadius(), 0., 0.)) *
        AngleAxis3D(-M_PI / 2, Vector3D(1., 0., 0.)));
    rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
        std::shared_ptr<const Transform3D>(sp2Transform), sp12Bounds));
  }
  return rSurfaces;
}

std::shared_ptr<const Acts::CylinderBounds>
Acts::CylinderVolumeBounds::innerCylinderBounds() const {
  return std::make_shared<const CylinderBounds>(
      m_valueStore.at(bv_innerRadius), m_valueStore.at(bv_halfPhiSector),
      m_valueStore.at(bv_halfZ));
}

std::shared_ptr<const Acts::CylinderBounds>
Acts::CylinderVolumeBounds::outerCylinderBounds() const {
  return std::make_shared<const CylinderBounds>(
      m_valueStore.at(bv_outerRadius), m_valueStore.at(bv_halfPhiSector),
      m_valueStore.at(bv_halfZ));
}

std::shared_ptr<const Acts::DiscBounds> Acts::CylinderVolumeBounds::discBounds()
    const {
  return std::shared_ptr<const DiscBounds>(new RadialBounds(
      m_valueStore.at(bv_innerRadius), m_valueStore.at(bv_outerRadius),
      m_valueStore.at(bv_halfPhiSector)));
}

std::shared_ptr<const Acts::PlanarBounds>
Acts::CylinderVolumeBounds::sectorPlaneBounds() const {
  return std::shared_ptr<const PlanarBounds>(new RectangleBounds(
      0.5 * (m_valueStore.at(bv_outerRadius) - m_valueStore.at(bv_innerRadius)),
      m_valueStore.at(bv_halfZ)));
}

std::ostream& Acts::CylinderVolumeBounds::toStream(std::ostream& sl) const {
  return dumpT<std::ostream>(sl);
}

Acts::Volume::BoundingBox Acts::CylinderVolumeBounds::boundingBox(
    const Transform3D* trf, const Vector3D& envelope,
    const Volume* entity) const {
  double xmax, xmin, ymax, ymin;
  xmax = outerRadius();

  if (halfPhiSector() > M_PI / 2.) {
    // more than half
    ymax = outerRadius();
    ymin = -outerRadius();
    xmin = outerRadius() * std::cos(halfPhiSector());
  } else {
    // less than half
    ymax = outerRadius() * std::sin(halfPhiSector());
    ymin = -ymax;
    // in this case, xmin is given by the inner radius
    xmin = innerRadius() * std::cos(halfPhiSector());
  }

  Vector3D vmin(xmin, ymin, -halflengthZ());
  Vector3D vmax(xmax, ymax, halflengthZ());

  // this is probably not perfect, but at least conservative
  Volume::BoundingBox box{entity, vmin - envelope, vmax + envelope};
  return trf == nullptr ? box : box.transformed(*trf);
}
