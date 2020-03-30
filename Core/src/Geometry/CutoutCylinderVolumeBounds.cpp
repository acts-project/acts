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
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
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

std::vector<std::shared_ptr<const Acts::Surface>>
Acts::CutoutCylinderVolumeBounds::decomposeToSurfaces(
    const Transform3D* transform) const {
  std::vector<std::shared_ptr<const Acts::Surface>> surfaces;

  // transform copy
  std::shared_ptr<const Transform3D> trf;
  if (transform != nullptr) {
    trf = std::make_shared<const Transform3D>(*transform);
  } else {
    trf = std::make_shared<const Transform3D>(Transform3D::Identity());
  }

  if (get(eMinR) == 0.) {
    surfaces.resize(6);  // exactly six surfaces (no choke inner cover)
  } else {
    surfaces.resize(8);  // exactly eight surfaces
  }

  // outer cylinder envelope
  auto outerBounds =
      std::make_shared<CylinderBounds>(get(eMaxR), get(eHalfLengthZ));
  auto outer = Surface::makeShared<CylinderSurface>(trf, outerBounds);
  surfaces.at(tubeOuterCover) = outer;

  // inner (small) cylinder envelope
  auto ctr_innerBounds =
      std::make_shared<CylinderBounds>(get(eMedR), get(eHalfLengthZcutout));
  auto ctr_inner = Surface::makeShared<CylinderSurface>(trf, ctr_innerBounds);
  surfaces.at(tubeInnerCover) = ctr_inner;

  // z position of the pos and neg choke points
  double hlChoke = (get(eHalfLengthZ) - get(eHalfLengthZcutout)) * 0.5;
  double zChoke = get(eHalfLengthZcutout) + hlChoke;

  if (get(eMinR) > 0.) {
    auto posChokeTrf = std::make_shared<const Transform3D>(
        *trf * Translation3D(Vector3D(0, 0, zChoke)));
    auto posInner =
        Surface::makeShared<CylinderSurface>(posChokeTrf, get(eMinR), hlChoke);
    surfaces.at(index7) = posInner;

    auto negChokeTrf = std::make_shared<const Transform3D>(
        *trf * Translation3D(Vector3D(0, 0, -zChoke)));
    auto negInner =
        Surface::makeShared<CylinderSurface>(negChokeTrf, get(eMinR), hlChoke);
    surfaces.at(index6) = negInner;
  }

  // outer disks
  auto posOutDiscTrf = std::make_shared<const Transform3D>(
      *trf * Translation3D(Vector3D(0, 0, get(eHalfLengthZ))));
  auto posOutDisc =
      Surface::makeShared<DiscSurface>(posOutDiscTrf, get(eMinR), get(eMaxR));
  surfaces.at(positiveFaceXY) = posOutDisc;

  auto negOutDiscTrf = std::make_shared<const Transform3D>(
      *trf * Translation3D(Vector3D(0, 0, -get(eHalfLengthZ))) *
      AngleAxis3D(M_PI, Vector3D::UnitX()));
  auto negOutDisc =
      Surface::makeShared<DiscSurface>(negOutDiscTrf, get(eMinR), get(eMaxR));
  surfaces.at(negativeFaceXY) = negOutDisc;

  // inner disks
  auto posInDiscTrf = std::make_shared<const Transform3D>(
      *trf * Translation3D(Vector3D(0, 0, get(eHalfLengthZcutout))));
  auto posInDisc =
      Surface::makeShared<DiscSurface>(posInDiscTrf, get(eMinR), get(eMedR));
  surfaces.at(index5) = posInDisc;

  auto negInDiscTrf = std::make_shared<const Transform3D>(
      *trf * Translation3D(Vector3D(0, 0, -get(eHalfLengthZcutout))) *
      AngleAxis3D(M_PI, Vector3D::UnitX()));
  auto negInDisc =
      Surface::makeShared<DiscSurface>(negInDiscTrf, get(eMinR), get(eMedR));
  surfaces.at(index4) = negInDisc;

  return surfaces;
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
