// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderVolumeBounds.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Volumes/CylinderVolumeBounds.hpp"
#include <cmath>
#include <iostream>
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"

const double Acts::CylinderVolumeBounds::s_numericalStable = 10e-2;

Acts::CylinderVolumeBounds::CylinderVolumeBounds()
  : VolumeBounds(), m_valueStore(4, 0.)
{
}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(double radius, double halez)
  : VolumeBounds(), m_valueStore(4, 0.)
{
  m_valueStore.at(bv_innerRadius)   = 0.;
  m_valueStore.at(bv_outerRadius)   = std::abs(radius);
  m_valueStore.at(bv_halfPhiSector) = M_PI;
  m_valueStore.at(bv_halfZ)         = std::abs(halez);
}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(double rinner,
                                                 double router,
                                                 double halez)
  : VolumeBounds(), m_valueStore(4, 0.)
{
  m_valueStore.at(bv_innerRadius)   = std::abs(rinner);
  m_valueStore.at(bv_outerRadius)   = std::abs(router);
  m_valueStore.at(bv_halfPhiSector) = M_PI;
  m_valueStore.at(bv_halfZ)         = std::abs(halez);
}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(double rinner,
                                                 double router,
                                                 double haphi,
                                                 double halez)
  : VolumeBounds(), m_valueStore(4, 0.)
{
  m_valueStore.at(bv_innerRadius)   = std::abs(rinner);
  m_valueStore.at(bv_outerRadius)   = std::abs(router);
  m_valueStore.at(bv_halfPhiSector) = std::abs(haphi);
  m_valueStore.at(bv_halfZ)         = std::abs(halez);
}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(
    const CylinderVolumeBounds& cylbo)
  : VolumeBounds(), m_valueStore(cylbo.m_valueStore)
{
}

Acts::CylinderVolumeBounds::~CylinderVolumeBounds() = default;

Acts::CylinderVolumeBounds&
Acts::CylinderVolumeBounds::operator=(const CylinderVolumeBounds& cylbo)
{
  if (this != &cylbo) {
    m_valueStore = cylbo.m_valueStore;
  }
  return *this;
}

std::vector<std::shared_ptr<const Acts::Surface>>
Acts::CylinderVolumeBounds::decomposeToSurfaces(
    std::shared_ptr<const Transform3D> transformPtr) const
{
  std::vector<std::shared_ptr<const Surface>> rSurfaces;
  rSurfaces.reserve(6);

  // set the transform
  Transform3D transform = (transformPtr == nullptr) ? Transform3D::Identity()
                                                    : (*transformPtr.get());
  const Transform3D* tTransform = nullptr;
  Vector3D           cylCenter(transform.translation());

  std::shared_ptr<const DiscBounds> dBounds = discBounds();
  // bottom Disc (negative z)
  tTransform
      = new Transform3D(transform * AngleAxis3D(M_PI, Vector3D(1., 0., 0.))
                        * Translation3D(Vector3D(0., 0., halflengthZ())));
  rSurfaces.push_back(Surface::makeShared<DiscSurface>(
      std::shared_ptr<const Transform3D>(tTransform), dBounds));
  // top Disc (positive z)
  tTransform = new Transform3D(
      transform * Translation3D(Vector3D(0., 0., halflengthZ())));
  rSurfaces.push_back(Surface::makeShared<DiscSurface>(
      std::shared_ptr<const Transform3D>(tTransform), dBounds));

  // outer Cylinder - shares the transform
  rSurfaces.push_back(Surface::makeShared<CylinderSurface>(
      transformPtr, outerCylinderBounds()));

  // innermost Cylinder
  if (innerRadius() > s_numericalStable) {
    rSurfaces.push_back(Surface::makeShared<CylinderSurface>(
        transformPtr, innerCylinderBounds()));
  }

  // the cylinder is sectoral
  if (std::abs(halfPhiSector() - M_PI) > s_numericalStable) {
    std::shared_ptr<const PlanarBounds> sp12Bounds = sectorPlaneBounds();
    // sectorPlane 1 (negative phi)
    const Transform3D* sp1Transform = new Transform3D(
        transform * AngleAxis3D(-halfPhiSector(), Vector3D(0., 0., 1.))
        * Translation3D(Vector3D(mediumRadius(), 0., 0.))
        * AngleAxis3D(M_PI / 2, Vector3D(1., 0., 0.)));
    rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
        std::shared_ptr<const Transform3D>(sp1Transform), sp12Bounds));
    // sectorPlane 2 (positive phi)
    const Transform3D* sp2Transform = new Transform3D(
        transform * AngleAxis3D(halfPhiSector(), Vector3D(0., 0., 1.))
        * Translation3D(Vector3D(mediumRadius(), 0., 0.))
        * AngleAxis3D(-M_PI / 2, Vector3D(1., 0., 0.)));
    rSurfaces.push_back(Surface::makeShared<PlaneSurface>(
        std::shared_ptr<const Transform3D>(sp2Transform), sp12Bounds));
  }
  return rSurfaces;
}

std::shared_ptr<const Acts::CylinderBounds>
Acts::CylinderVolumeBounds::innerCylinderBounds() const
{
  return std::make_shared<const CylinderBounds>(
      m_valueStore.at(bv_innerRadius),
      m_valueStore.at(bv_halfPhiSector),
      m_valueStore.at(bv_halfZ));
}

std::shared_ptr<const Acts::CylinderBounds>
Acts::CylinderVolumeBounds::outerCylinderBounds() const
{
  return std::make_shared<const CylinderBounds>(
      m_valueStore.at(bv_outerRadius),
      m_valueStore.at(bv_halfPhiSector),
      m_valueStore.at(bv_halfZ));
}

std::shared_ptr<const Acts::DiscBounds>
Acts::CylinderVolumeBounds::discBounds() const
{
  return std::shared_ptr<const DiscBounds>(
      new RadialBounds(m_valueStore.at(bv_innerRadius),
                       m_valueStore.at(bv_outerRadius),
                       m_valueStore.at(bv_halfPhiSector)));
}

std::shared_ptr<const Acts::PlanarBounds>
Acts::CylinderVolumeBounds::sectorPlaneBounds() const
{
  return std::shared_ptr<const PlanarBounds>(new RectangleBounds(
      0.5 * (m_valueStore.at(bv_outerRadius) - m_valueStore.at(bv_innerRadius)),
      m_valueStore.at(bv_halfZ)));
}

std::ostream&
Acts::CylinderVolumeBounds::dump(std::ostream& sl) const
{
  return dumpT<std::ostream>(sl);
}
