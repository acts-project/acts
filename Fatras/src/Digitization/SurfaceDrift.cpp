// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Digitization/SurfaceDrift.hpp"

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsFatras/Digitization/DigitizationError.hpp"

#include <cmath>

namespace ActsFatras {

namespace {

/// Express the hit position and direction in the readout-local 3D frame.
///
/// The returned position has its third component (the surface-normal
/// coordinate) set to 0 — the readout reference plane — while the direction is
/// a unit vector whose third component is the projection onto the surface
/// normal. The first two components are the in-plane readout coordinates:
///   - Plane / Disc : Cartesian local (x, y)
///   - Cylinder     : unrolled (rPhi, z)
///
/// Returns false if the local frame cannot be established.
bool toLocalFrame(const Acts::GeometryContext& gctx,
                  const Acts::Surface& surface, const Acts::Vector3& pos,
                  const Acts::Vector3& dir, Acts::Vector3& pos3Local,
                  Acts::Vector3& dir3Local) {
  const auto invTransform = surface.localToGlobalTransform(gctx).inverse();

  if (surface.type() == Acts::Surface::SurfaceType::Cylinder) {
    // Position/direction in the cylinder-axis frame (axis aligned with local z)
    const Acts::Vector3 posCyl = invTransform * pos;
    const Acts::Vector3 dirCyl = invTransform.linear() * dir.normalized();

    const double phi = std::atan2(posCyl.y(), posCyl.x());
    const double R = std::hypot(posCyl.x(), posCyl.y());
    const double cphi = std::cos(phi);
    const double sphi = std::sin(phi);

    // Orthonormal readout-local basis at the hit:
    //   tangential = (-sinφ, cosφ, 0), axial = (0, 0, 1), radial = (cosφ, sinφ,
    //   0)
    const double dsRPhi = -sphi * dirCyl.x() + cphi * dirCyl.y();
    const double dsZ = dirCyl.z();
    const double dsR = cphi * dirCyl.x() + sphi * dirCyl.y();

    // (rPhi, z) readout position, normal coordinate = 0
    pos3Local = Acts::Vector3(R * phi, posCyl.z(), 0.);
    // (tangential, axial, radial); radial is the surface-normal component
    dir3Local = Acts::Vector3(dsRPhi, dsZ, dsR);
    return true;
  }

  // Plane / Disc: the Cartesian local frame, surface normal = local z.
  pos3Local = invTransform * pos;
  dir3Local = invTransform.linear() * dir.normalized();
  return true;
}

}  // namespace

Acts::Result<std::tuple<SurfaceDrift::Segment2D, SurfaceDrift::Segment3D>>
SurfaceDrift::toReadout(const Acts::GeometryContext& gctx,
                        const Acts::Surface& surface, double thickness,
                        const Acts::Vector3& pos, const Acts::Vector3& dir,
                        const Acts::Vector3& driftDir) const {
  Acts::Vector3 pos3Local = Acts::Vector3::Zero();
  Acts::Vector3 seg3Local = Acts::Vector3::Zero();
  if (!toLocalFrame(gctx, surface, pos, dir, pos3Local, seg3Local)) {
    return DigitizationError::DriftError;
  }

  // seg3Local is a unit direction; its z-component is the projection onto the
  // surface normal. If it (nearly) vanishes the track is parallel to the
  // surface and the drift is undefined.
  if (std::abs(seg3Local.z()) < Acts::s_epsilon) {
    return DigitizationError::DriftError;
  }

  // Scale the unit vector to the thickness of the module
  const double scale = thickness / seg3Local.z();
  seg3Local *= scale;
  // The drift direction is in the local frame, so we need to transform it
  const Acts::Vector3 entry = pos3Local - 0.5 * seg3Local;
  const Acts::Vector3 exit = pos3Local + 0.5 * seg3Local;
  Acts::Vector3 driftedEntry = entry;
  Acts::Vector3 driftedExit = exit;
  // Apply a Lorentz drift if the drift direction has both an in-plane and a
  // perpendicular (normal) component. driftDir is expressed in the same
  // readout-local frame as seg3Local.
  if (Acts::VectorHelpers::perp(driftDir) > Acts::s_epsilon &&
      std::abs(driftDir.z()) > Acts::s_epsilon) {
    // Apply the drift to the entry and exit points
    const double driftScale = 0.5 * thickness / driftDir.z();
    driftedEntry += driftScale * driftDir;
    driftedExit -= driftScale * driftDir;
  }

  // The drifted 2D readout segment and the undrifted 3D segment (whose length
  // is the true chord through the depletion zone, used downstream to rescale
  // the 2D activation path to 3D).
  return std::tuple{
      Segment2D{driftedEntry.segment<2>(0), driftedExit.segment<2>(0)},
      Segment3D{entry, exit}};
}

}  // namespace ActsFatras
