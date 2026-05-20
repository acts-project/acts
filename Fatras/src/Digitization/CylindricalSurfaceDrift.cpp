// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Digitization/CylindricalSurfaceDrift.hpp"

#include <cmath>

namespace {

/// Build the readout-local frame at the hit point on the cylinder.
///
/// Returns the 2-D position (rPhi, z), the unsigned tangential and axial
/// components of the direction (ds_rphi, ds_z) and the absolute value of the
/// radial component (|ds_r|) — all expressed in the surface's intrinsic
/// (axis-aligned) frame.
struct LocalDecomposition {
  Acts::Vector2 pos2D;
  double dsRPhi;
  double dsZ;
  double absDsR;
};

LocalDecomposition decompose(const Acts::Transform3& invTransform,
                             const Acts::Vector3& pos,
                             const Acts::Vector3& dir) {
  // Position in cylinder-axis frame: cylinder axis aligned with local z,
  // radial direction in the local xy-plane.
  const Acts::Vector3 posCyl = invTransform * pos;
  // Direction in cylinder-axis frame.
  const Acts::Vector3 dirCyl = invTransform.linear() * dir;

  const double phi = std::atan2(posCyl.y(), posCyl.x());
  const double R = std::hypot(posCyl.x(), posCyl.y());

  const double cphi = std::cos(phi);
  const double sphi = std::sin(phi);

  // Tangential unit vector at the hit (in the cylinder-axis frame).
  // Axial unit vector is (0, 0, 1).
  // Radial unit vector is (cosφ, sinφ, 0).
  const double dsRPhi = -sphi * dirCyl.x() + cphi * dirCyl.y();
  const double dsZ = dirCyl.z();
  const double dsR = cphi * dirCyl.x() + sphi * dirCyl.y();

  return {Acts::Vector2(R * phi, posCyl.z()), dsRPhi, dsZ, std::abs(dsR)};
}

}  // namespace

ActsFatras::CylindricalSurfaceDrift::Segment2D
ActsFatras::CylindricalSurfaceDrift::toReadout(
    const Acts::GeometryContext& gctx, const Acts::Surface& surface,
    double thickness, const Acts::Vector3& pos, const Acts::Vector3& dir,
    const Acts::Vector3& driftDir) const {
  const auto invTransform = surface.localToGlobalTransform(gctx).inverse();
  const auto decomp = decompose(invTransform, pos, dir);

  // Scale the unit direction to span the depletion thickness in the radial
  // direction. If the track is nearly tangent to the cylinder we clamp the
  // step length to the thickness itself to avoid blow-up.
  constexpr double kMinAbsRadial = 1e-12;
  const double scale =
      (decomp.absDsR > kMinAbsRadial) ? (thickness / decomp.absDsR) : thickness;
  const Acts::Vector2 seg2D(decomp.dsRPhi * scale, decomp.dsZ * scale);

  // Entry/exit before applying any Lorentz drift.
  Acts::Vector2 entry = decomp.pos2D - 0.5 * seg2D;
  Acts::Vector2 exit = decomp.pos2D + 0.5 * seg2D;

  // Optional Lorentz drift. The drift direction is interpreted in the readout
  // local 3-D frame: x = tangential (rPhi), y = axial (z), z = radial (depth).
  // Same convention as the planar case: driftDir.z() picks the readout side.
  if (!driftDir.segment<2>(0).isApprox(Acts::Vector2(0., 0.))) {
    const double absDriftZ = std::abs(driftDir.z());
    const double driftScale =
        (absDriftZ > kMinAbsRadial) ? (thickness / absDriftZ) : thickness;
    const Acts::Vector2 driftShift = driftDir.segment<2>(0) * driftScale;

    if (driftDir.z() > 0.) {
      entry += driftShift;
    }
    if (driftDir.z() < 0.) {
      exit += driftShift;
    }
  }

  return {entry, exit};
}
