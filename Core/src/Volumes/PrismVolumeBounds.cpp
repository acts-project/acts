// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PrismVolumeBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Volumes/PrismVolumeBounds.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Surfaces/TriangleBounds.hpp"
#include "ACTS/Volumes/Volume.hpp"
#include <iomanip>
#include <iostream>
#include <math.h>

Acts::PrismVolumeBounds::PrismVolumeBounds()
  : VolumeBounds(), m_xyVtx(), m_halfZ(), m_baseBounds(nullptr), m_ordering(-1)
{
}

Acts::PrismVolumeBounds::PrismVolumeBounds(
    std::vector< Vector2D > xyVtx,
    float halez)
  : VolumeBounds(), m_xyVtx(xyVtx), m_halfZ(halez), m_baseBounds(new TriangleBounds(xyVtx)), m_ordering(-1)
{
}


Acts::PrismVolumeBounds::PrismVolumeBounds(const Acts::PrismVolumeBounds& trabo)
  : VolumeBounds()
  , m_halfZ(trabo.m_halfZ)
  , m_baseBounds(trabo.m_baseBounds)
  , m_ordering(trabo.m_ordering)
{
  m_xyVtx.resize(trabo.m_xyVtx.size());
  m_xyVtx.assign(trabo.m_xyVtx.begin(), trabo.m_xyVtx.end());
}

Acts::PrismVolumeBounds::~PrismVolumeBounds()
{
}

Acts::PrismVolumeBounds&
Acts::PrismVolumeBounds::operator=(const Acts::PrismVolumeBounds& trabo)
{
  if (this != &trabo) {
    m_halfZ       = trabo.m_halfZ;
    m_xyVtx       = trabo.m_xyVtx;
    m_baseBounds  = trabo.m_baseBounds;
    m_ordering    = trabo.m_ordering;
  }
  return *this;
}

const std::vector<const Acts::Surface*>
Acts::PrismVolumeBounds::decomposeToSurfaces(
    std::shared_ptr<Acts::Transform3D> transformPtr) const
{
  std::vector<const Surface*> rSurfaces;

  // the transform
  Transform3D transform = (transformPtr == nullptr)
      ? Transform3D::Identity()
      : (*transformPtr.get());
  Transform3D* tTransform = 0;
  // face surfaces xy
  //  (1) - at positive local z
  tTransform = new Transform3D(
      transform * Translation3D(Vector3D(0., 0., m_halfZ)));
  
  
  PlaneSurface* xyPlane
      = new PlaneSurface(std::shared_ptr<Transform3D>(tTransform),
                         m_baseBounds);
  rSurfaces.push_back(xyPlane);
  //  (2) - at negative local z
  tTransform = new Transform3D(
      transform * Translation3D(Vector3D(0., 0., -m_halfZ))
      * AngleAxis3D(M_PI, Vector3D(1., 0., 0.)));
  PlaneSurface* xymPlane
      = new PlaneSurface(std::shared_ptr<Transform3D>(tTransform),
                         std::shared_ptr<const PlanarBounds>(new TriangleBounds(mirror_xyVtx())));
  rSurfaces.push_back(xymPlane);
  // loop over xy vertices
  //  (3)
  for (unsigned int iv = 0; iv < m_xyVtx.size(); iv++) {
    if (iv != m_xyVtx.size() - 1)
      rSurfaces.push_back(sideSurf(transform, iv, iv + 1));
    else
      rSurfaces.push_back(sideSurf(transform, iv, 0));
  }

  return rSurfaces;
}

// faces in xy
PlaneSurface*
PrismVolumeBounds::sideSurf(Transform3D transform,
                            unsigned int      iv1,
                            unsigned int      iv2) const
{
  PlaneSurface* plane = 0;

  double xdif  = m_xyVtx.at(iv2).first - m_xyVtx.at(iv1).first;
  double ydif  = m_xyVtx.at(iv2).second - m_xyVtx.at(iv1).second;
  double xsize = sqrt(xdif * xdif + ydif * ydif);

  double ori = ordering() > 0 ? 1. : -1.;

  Vector3D pos(0.5 * (m_xyVtx.at(iv1).first + m_xyVtx.at(iv2).first),
                     0.5 * (m_xyVtx.at(iv1).second + m_xyVtx.at(iv2).second),
                     0.);
  double phi                   = ori * ydif < 0 ? M_PI / 2 : -M_PI / 2;
  if (ori > 0 && ydif > 0) phi = M_PI / 2;
  if (fabs(xdif) > 1e-6) {
    phi = atan(ydif / xdif);
    if (xdif < 0) phi += M_PI;
  }

  Transform3D* tTransform = new Transform3D(
      transform * Translation3D(pos)
      * AngleAxis3D(phi, Vector3D(0., 0., 1.))
      * AngleAxis3D(-ori * 0.5 * M_PI, Vector3D(1., 0., 0.)));
  plane
      = new PlaneSurface(std::shared_ptr<Transform3D>(tTransform),
                         std::shared_ptr<const PlanarBounds>
                           (new RectangleBounds(0.5 * xsize, m_halfZ)));

  // orientation
  int ivr                                     = (iv1 == 0 || iv2 == 0) ? 1 : 0;
  if (ivr == 1 && (iv1 == 1 || iv2 == 1)) ivr = 2;

  double         ox = m_xyVtx.at(ivr).first - pos[0];
  double         oy = m_xyVtx.at(ivr).second - pos[1];
  Vector3D d(ox, oy, 0.);

  // protect against wrong orientation
  if (d.dot(plane->normal()) > 0.) {
    delete plane;
    tTransform = new Transform3D(
        transform * Translation3D(pos)
        * AngleAxis3D(phi + M_PI, Vector3D(0., 0., 1.))
        * AngleAxis3D(-ori * 0.5 * M_PI, Vector3D(1., 0., 0.)));
    plane = new PlaneSurface(
        std::shared_ptr<Transform3D>(tTransform),
        std::shared_ptr<const PlanrBounds>(new RectangleBounds(0.5 * xsize, m_halfZ)));
  }

  return plane;
}

bool
PrismVolumeBounds::inside(const Vector3D& pos, double tol) const
{
  if (fabs(pos.z()) > m_halfZ + tol) return false;
  // xy plane
  Vector2D locp(pos.x(), pos.y());
  return (m_baseBounds->inside(locp, BoundaryCheck(true, true, tol, tol)));
}

std::vector<std::pair<TDD_real_t, TDD_real_t>>
PrismVolumeBounds::mirror_xyVtx() const
{
  // prepare the mirrored vertices
  std::vector<std::pair<TDD_real_t, TDD_real_t>> mirrored;
  mirrored.resize(m_xyVtx.size());
  // flip the second coordinate
  for (unsigned int i = 0; i < m_xyVtx.size(); i++)
    mirrored.at(i)    = std::pair<TDD_real_t, TDD_real_t>(m_xyVtx.at(i).first,
                                                       -m_xyVtx.at(i).second);
  // return the mirroed ones
  return mirrored;
}

int
PrismVolumeBounds::ordering() const
{
  if (m_ordering > -1) return m_ordering;

  m_ordering = 1;

  double yd2 = m_xyVtx.at(2).second - m_xyVtx.at(1).second;
  double yd0 = m_xyVtx.at(0).second - m_xyVtx.at(1).second;
  double xd2 = m_xyVtx.at(2).first - m_xyVtx.at(1).first;
  double xd0 = m_xyVtx.at(0).first - m_xyVtx.at(1).first;
  double ph2 = yd2 < 0 ? -M_PI / 2 : M_PI / 2;
  if (fabs(xd2) > 1e-6) {
    ph2 = atan(yd2 / xd2);
    if (xd2 < 0) ph2 += M_PI;
  }
  double ph0 = yd0 < 0 ? -M_PI / 2 : M_PI / 2;
  if (fabs(xd0) > 1e-6) {
    ph0 = atan(yd0 / xd0);
    if (xd0 < 0) ph0 += M_PI;
  }
  if (ph0 < 0) ph0 += 2 * M_PI;
  if (ph2 < 0) ph2 += 2 * M_PI;

  if ((ph0 > ph2 && ph0 - ph2 < M_PI) || (ph2 - ph0) > M_PI) m_ordering = 0;

  return m_ordering;
}

// ostream operator overload
std::ostream&
PrismVolumeBounds::dump(std::ostream& sl) const
{
  return dumpT<std::ostream>(sl);
}
