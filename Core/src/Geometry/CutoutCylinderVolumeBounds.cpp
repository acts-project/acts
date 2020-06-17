// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CutoutCylinderVolumeBounds.hpp"
#include "Acts/Geometry/BoundarySurfaceFace.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Visualization/IVisualization.hpp"

#include <memory>
#include <vector>

bool Acts::CutoutCylinderVolumeBounds::inside(const Acts::Vector3D& gpos,
                                              double tol) const {
  // first check whether we are in the outer envelope at all (ignore r_med)
  using VectorHelpers::perp;
  using VectorHelpers::phi;
  double ros = perp(gpos);

  bool insideR = (ros >= get(eMinR) - tol) && (ros <= get(eMaxR) + tol);
  bool insideZ = std::abs(gpos.z()) <= get(eHalfLengthZ) + tol;

  if (!insideR || !insideZ) {
    return false;
  }

  // we're inside the outer volume, but we might be in inside the
  // cutout section in the middle
  bool insideRInner = ros <= get(eMedR) - tol;
  bool insideZInner = std::abs(gpos.z()) < get(eHalfLengthZcutout) - tol;

  return !insideRInner || !insideZInner;  // we are not, inside bounds
}

Acts::OrientedSurfaces Acts::CutoutCylinderVolumeBounds::orientedSurfaces(
    const Transform3D* transform) const {
  OrientedSurfaces oSurfaces;

  // Transform copy
  std::shared_ptr<const Transform3D> trf;
  if (transform != nullptr) {
    trf = std::make_shared<const Transform3D>(*transform);
  } else {
    trf = std::make_shared<const Transform3D>(Transform3D::Identity());
  }

  if (get(eMinR) == 0.) {
    oSurfaces.resize(6);  // exactly six surfaces (no choke inner cover)
  } else {
    oSurfaces.resize(8);  // exactly eight surfaces
  }

  // Outer cylinder envelope
  auto outer = Surface::makeShared<CylinderSurface>(trf, m_outerCylinderBounds);
  oSurfaces.at(tubeOuterCover) = OrientedSurface(std::move(outer), backward);

  // Inner (cutout) cylinder envelope
  auto cutoutInner =
      Surface::makeShared<CylinderSurface>(trf, m_cutoutCylinderBounds);
  oSurfaces.at(tubeInnerCover) =
      OrientedSurface(std::move(cutoutInner), forward);

  // z position of the pos and neg choke points
  double hlChoke = (get(eHalfLengthZ) - get(eHalfLengthZcutout)) * 0.5;
  double zChoke = get(eHalfLengthZcutout) + hlChoke;

  if (m_innerCylinderBounds != nullptr) {
    auto posChokeTrf = std::make_shared<const Transform3D>(
        *trf * Translation3D(Vector3D(0, 0, zChoke)));
    auto posInner = Surface::makeShared<CylinderSurface>(posChokeTrf,
                                                         m_innerCylinderBounds);
    oSurfaces.at(index7) = OrientedSurface(std::move(posInner), forward);

    auto negChokeTrf = std::make_shared<const Transform3D>(
        *trf * Translation3D(Vector3D(0, 0, -zChoke)));
    auto negInner = Surface::makeShared<CylinderSurface>(negChokeTrf,
                                                         m_innerCylinderBounds);
    oSurfaces.at(index6) = OrientedSurface(std::move(negInner), forward);
  }

  // Two Outer disks
  auto posOutDiscTrf = std::make_shared<const Transform3D>(
      *trf * Translation3D(Vector3D(0, 0, get(eHalfLengthZ))));
  auto posOutDisc =
      Surface::makeShared<DiscSurface>(posOutDiscTrf, m_outerDiscBounds);
  oSurfaces.at(positiveFaceXY) =
      OrientedSurface(std::move(posOutDisc), backward);

  auto negOutDiscTrf = std::make_shared<const Transform3D>(
      *trf * Translation3D(Vector3D(0, 0, -get(eHalfLengthZ))));

  auto negOutDisc =
      Surface::makeShared<DiscSurface>(negOutDiscTrf, m_outerDiscBounds);
  oSurfaces.at(negativeFaceXY) =
      OrientedSurface(std::move(negOutDisc), forward);

  // Two Inner disks
  auto posInDiscTrf = std::make_shared<const Transform3D>(
      *trf * Translation3D(Vector3D(0, 0, get(eHalfLengthZcutout))));
  auto posInDisc =
      Surface::makeShared<DiscSurface>(posInDiscTrf, m_innerDiscBounds);
  oSurfaces.at(index5) = OrientedSurface(std::move(posInDisc), forward);

  auto negInDiscTrf = std::make_shared<const Transform3D>(
      *trf * Translation3D(Vector3D(0, 0, -get(eHalfLengthZcutout))));
  auto negInDisc =
      Surface::makeShared<DiscSurface>(negInDiscTrf, m_innerDiscBounds);
  oSurfaces.at(index4) = OrientedSurface(std::move(negInDisc), backward);

  return oSurfaces;
}

Acts::Volume::BoundingBox Acts::CutoutCylinderVolumeBounds::boundingBox(
    const Acts::Transform3D* trf, const Acts::Vector3D& envelope,
    const Acts::Volume* entity) const {
  Vector3D vmin, vmax;

  // no phi sector is possible, so this is just the outer size of
  // the cylinder

  vmax = {get(eMaxR), get(eMaxR), get(eHalfLengthZ)};
  vmin = {-get(eMaxR), -get(eMaxR), -get(eHalfLengthZ)};

  Acts::Volume::BoundingBox box(entity, vmin - envelope, vmax + envelope);
  // transform at the very end, if required
  return trf == nullptr ? box : box.transformed(*trf);
}

std::ostream& Acts::CutoutCylinderVolumeBounds::toStream(
    std::ostream& sl) const {
  sl << "Acts::CutoutCylinderVolumeBounds(\n";
  sl << "rmin = " << get(eMinR) << " rmed = " << get(eMedR)
     << " rmax = " << get(eMaxR) << "\n";
  sl << "dz1 = " << get(eHalfLengthZ) << " dz2 = " << get(eHalfLengthZcutout);
  return sl;
}

void Acts::CutoutCylinderVolumeBounds::buildSurfaceBounds() {
  if (get(eMinR) > s_epsilon) {
    double hlChoke = (get(eHalfLengthZ) - get(eHalfLengthZcutout)) * 0.5;
    m_innerCylinderBounds =
        std::make_shared<CylinderBounds>(get(eMinR), hlChoke);
  }

  m_cutoutCylinderBounds =
      std::make_shared<CylinderBounds>(get(eMedR), get(eHalfLengthZcutout));

  m_outerCylinderBounds =
      std::make_shared<CylinderBounds>(get(eMaxR), get(eHalfLengthZ));

  m_innerDiscBounds = std::make_shared<RadialBounds>(get(eMinR), get(eMedR));

  m_outerDiscBounds = std::make_shared<RadialBounds>(get(eMinR), get(eMaxR));
}
