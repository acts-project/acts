// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Digitization/CylindricalSurfaceMask.hpp"

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "ActsFatras/Digitization/DigitizationError.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace {

/// Liang–Barsky clipping of a segment against an axis-aligned rectangle.
///
/// On success returns a {tEnter, tExit} pair in [0, 1] describing the
/// parameter range of the clipped segment along (end - start). If the segment
/// lies completely outside the rectangle the returned pair has tEnter > tExit.
struct ClipResult {
  double tEnter;
  double tExit;
  bool valid;
};

ClipResult clipLiangBarsky(double x0, double y0, double x1, double y1,
                           double xmin, double xmax, double ymin, double ymax) {
  const double dx = x1 - x0;
  const double dy = y1 - y0;

  double tEnter = 0.0;
  double tExit = 1.0;

  const double p[4] = {-dx, dx, -dy, dy};
  const double q[4] = {x0 - xmin, xmax - x0, y0 - ymin, ymax - y0};

  for (int i = 0; i < 4; ++i) {
    if (p[i] == 0.0) {
      // Segment is parallel to this clipping edge.
      if (q[i] < 0.0) {
        return {0.0, 0.0, false};
      }
      continue;
    }
    const double t = q[i] / p[i];
    if (p[i] < 0.0) {
      tEnter = std::max(tEnter, t);
    } else {
      tExit = std::min(tExit, t);
    }
  }

  if (tEnter > tExit) {
    return {tEnter, tExit, false};
  }
  return {tEnter, tExit, true};
}

}  // namespace

Acts::Result<ActsFatras::CylindricalSurfaceMask::Segment2D>
ActsFatras::CylindricalSurfaceMask::apply(const Acts::Surface& surface,
                                          const Segment2D& segment) const {
  if (surface.type() != Acts::Surface::Cylinder ||
      surface.bounds().type() != Acts::SurfaceBounds::eCylinder) {
    return DigitizationError::UndefinedSurface;
  }

  const auto& cBounds =
      static_cast<const Acts::CylinderBounds&>(surface.bounds());

  const double R = cBounds.get(Acts::CylinderBounds::eR);
  const double halfZ = cBounds.get(Acts::CylinderBounds::eHalfLengthZ);
  const double halfPhi = cBounds.get(Acts::CylinderBounds::eHalfPhiSector);
  const double avgPhi = cBounds.get(Acts::CylinderBounds::eAveragePhi);

  // (rPhi, z) bounding box of the cylinder bounds.
  const double rPhiMin = R * (avgPhi - halfPhi);
  const double rPhiMax = R * (avgPhi + halfPhi);
  const double zMin = -halfZ;
  const double zMax = halfZ;

  // Local rPhi may be reported in the (-π R, π R] principal branch by the
  // drift step. If the sector is closed (full cylinder) any rPhi is valid and
  // we only need to clip in z. Otherwise unwrap by adding/subtracting 2π·R so
  // that the segment lies in the same branch as the sector midpoint.
  Segment2D seg = segment;
  const bool isClosed = std::abs(halfPhi - std::numbers::pi) < Acts::s_epsilon;

  if (!isClosed) {
    const double midRPhi = R * avgPhi;
    const double twoPiR = 2.0 * std::numbers::pi * R;
    for (auto& p : seg) {
      while (p[0] - midRPhi > std::numbers::pi * R) {
        p[0] -= twoPiR;
      }
      while (midRPhi - p[0] > std::numbers::pi * R) {
        p[0] += twoPiR;
      }
    }
  }

  // Fast exit: both endpoints inside.
  auto inside = [&](const Acts::Vector2& p) {
    const bool inZ = (p[1] >= zMin) && (p[1] <= zMax);
    if (isClosed) {
      return inZ;
    }
    return inZ && (p[0] >= rPhiMin) && (p[0] <= rPhiMax);
  };
  if (inside(seg[0]) && inside(seg[1])) {
    return seg;
  }

  // If the cylinder is closed in phi, the only clipping needed is in z.
  // Liang–Barsky with rPhi range = (-∞, +∞) reduces to the z-only case, but
  // implementing it explicitly avoids special-casing the rPhi parameter.
  const double useRPhiMin =
      isClosed ? -std::numeric_limits<double>::infinity() : rPhiMin;
  const double useRPhiMax =
      isClosed ? std::numeric_limits<double>::infinity() : rPhiMax;

  const auto clip = clipLiangBarsky(seg[0][0], seg[0][1], seg[1][0], seg[1][1],
                                    useRPhiMin, useRPhiMax, zMin, zMax);
  if (!clip.valid) {
    return DigitizationError::MaskingError;
  }

  const Acts::Vector2 d = seg[1] - seg[0];
  return Segment2D{seg[0] + clip.tEnter * d, seg[0] + clip.tExit * d};
}
