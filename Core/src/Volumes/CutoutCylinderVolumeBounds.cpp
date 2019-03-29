// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Volumes/CutoutCylinderVolumeBounds.hpp"
#include <memory>
#include <vector>
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/IVisualization.hpp"
#include "Acts/Volumes/BoundarySurfaceFace.hpp"
#include "Acts/Volumes/Volume.hpp"
#include "Acts/Volumes/VolumeBounds.hpp"

Acts::VolumeBounds*
Acts::CutoutCylinderVolumeBounds::clone() const
{
  return new CutoutCylinderVolumeBounds(*this);
}

bool
Acts::CutoutCylinderVolumeBounds::inside(const Acts::Vector3D& gpos,
                                         double                tol) const
{
  // first check whether we are in the outer envelope at all (ignore r_med)
  using VectorHelpers::perp;
  using VectorHelpers::phi;
  double ros = perp(gpos);

  bool insideR = (ros >= m_rmin - tol) && (ros <= m_rmax + tol);
  bool insideZ = std::abs(gpos.z()) <= m_dz1 + tol;

  if (!insideR || !insideZ) {
    return false;
  }

  // we're inside the outer volume, but we might be in inside the
  // cutout section in the middle
  bool insideRInner = ros <= m_rmed - tol;
  bool insideZInner = std::abs(gpos.z()) < m_dz2 - tol;

  return !insideRInner || !insideZInner;  // we are not, inside bounds
}

std::vector<std::shared_ptr<const Acts::Surface>>
Acts::CutoutCylinderVolumeBounds::decomposeToSurfaces(
    const Transform3D* transform) const
{
  std::vector<std::shared_ptr<const Acts::Surface>> surfaces;

  // transform copy
  std::shared_ptr<const Transform3D> trf;
  if (transform != nullptr) {
    trf = std::make_shared<const Transform3D>(*transform);
  } else {
    trf = std::make_shared<const Transform3D>(Transform3D::Identity());
  }

  if (m_rmin == 0.) {
    surfaces.resize(6);  // exactly six surfaces (no choke inner cover)
  } else {
    surfaces.resize(8);  // exactly eight surfaces
  }

  // outer cylinder envelope
  auto outer = Surface::makeShared<CylinderSurface>(trf, m_rmax, m_dz1);
  // surfaces.push_back(outer);
  surfaces.at(tubeOuterCover) = outer;

  // inner (small) cylinder envelope
  auto ctr_inner = Surface::makeShared<CylinderSurface>(trf, m_rmed, m_dz2);
  // surfaces.push_back(ctr_inner);
  surfaces.at(tubeInnerCover) = ctr_inner;

  // z position of the pos and neg choke points
  double hlChoke = (m_dz1 - m_dz2) * 0.5;
  double zChoke  = m_dz2 + hlChoke;

  if (m_rmin > 0.) {
    auto posChokeTrf = std::make_shared<const Transform3D>(
        *trf * Translation3D(Vector3D(0, 0, zChoke)));
    auto posInner
        = Surface::makeShared<CylinderSurface>(posChokeTrf, m_rmin, hlChoke);
    // surfaces.push_back(posInner);
    surfaces.at(index7) = posInner;

    auto negChokeTrf = std::make_shared<const Transform3D>(
        *trf * Translation3D(Vector3D(0, 0, -zChoke)));
    auto negInner
        = Surface::makeShared<CylinderSurface>(negChokeTrf, m_rmin, hlChoke);
    // surfaces.push_back(negInner);
    surfaces.at(index6) = negInner;
  }

  // outer disks
  auto posOutDiscTrf = std::make_shared<const Transform3D>(
      *trf * Translation3D(Vector3D(0, 0, m_dz1)));
  auto posOutDisc
      = Surface::makeShared<DiscSurface>(posOutDiscTrf, m_rmin, m_rmax);
  // surfaces.push_back(posOutDisc);
  surfaces.at(positiveFaceXY) = posOutDisc;

  auto negOutDiscTrf = std::make_shared<const Transform3D>(
      *trf * Translation3D(Vector3D(0, 0, -m_dz1))
      * AngleAxis3D(M_PI, Vector3D::UnitX()));
  auto negOutDisc
      = Surface::makeShared<DiscSurface>(negOutDiscTrf, m_rmin, m_rmax);
  // surfaces.push_back(negOutDisc);
  surfaces.at(negativeFaceXY) = negOutDisc;

  // inner disks
  auto posInDiscTrf = std::make_shared<const Transform3D>(
      *trf * Translation3D(Vector3D(0, 0, m_dz2)));
  auto posInDisc
      = Surface::makeShared<DiscSurface>(posInDiscTrf, m_rmin, m_rmed);
  // surfaces.push_back(posInDisc);
  surfaces.at(index5) = posInDisc;

  auto negInDiscTrf = std::make_shared<const Transform3D>(
      *trf * Translation3D(Vector3D(0, 0, -m_dz2))
      * AngleAxis3D(M_PI, Vector3D::UnitX()));
  auto negInDisc
      = Surface::makeShared<DiscSurface>(negInDiscTrf, m_rmin, m_rmed);
  // surfaces.push_back(negInDisc);
  surfaces.at(index4) = negInDisc;

  return surfaces;
}
std::ostream&
Acts::CutoutCylinderVolumeBounds::toStream(std::ostream& sl) const
{
  sl << "Acts::CutoutCylinderVolumeBounds(\n";
  sl << "rmin = " << m_rmin << " rmed = " << m_rmed << " rmax = " << m_rmax
     << "\n";
  sl << "dz1 = " << m_dz1 << " dz2 = " << m_dz2;
  return sl;
}

void
Acts::CutoutCylinderVolumeBounds::draw(IVisualization&    helper,
                                       const Transform3D& transform) const
{
  std::vector<std::shared_ptr<const Acts::Surface>> surfaces
      = decomposeToSurfaces(&transform);
  for (const auto& srf : surfaces) {
    auto cyl  = dynamic_cast<const CylinderSurface*>(srf.get());
    auto disc = dynamic_cast<const DiscSurface*>(srf.get());
    if (cyl != nullptr) {
      cyl->polyhedronRepresentation(50).draw(helper);
    } else {
      disc->polyhedronRepresentation(50).draw(helper);
    }
  }
}
