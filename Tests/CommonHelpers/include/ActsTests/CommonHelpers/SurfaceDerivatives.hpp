// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

namespace ActsTests {

/// Check the alignment derivative chain of @p surface against central finite
/// differences.
///
/// @p nominal must be the transform of @p surface, and @p position must lie on
/// the surface. @p rawDirection is normalized before use. A translation step
/// is translate(delta), which adds R*delta. A rotation step is a right-multiply
/// by AngleAxis3. The track ray stays fixed while the surface moves. The
/// checked derivatives are:
///
///   ds/da   alignmentToPathDerivative
///   db/df   freeToBoundJacobian
///   db/da   alignmentToBoundDerivative
///
/// A failure names the derivative and the first few mismatching entries.
/// @p makeSurface rebuilds the surface on a perturbed transform and returns a
/// pointer to it.
template <typename SurfaceFactory>
void checkSurfaceDerivatives(const Acts::GeometryContext& gctx,
                             const Acts::Surface& surface,
                             const Acts::Transform3& nominal,
                             const Acts::Vector3& position,
                             const Acts::Vector3& rawDirection,
                             SurfaceFactory makeSurface) {
  using namespace Acts;
  constexpr double kEps = 1e-6;
  constexpr double kTol = 1e-4;
  constexpr double kOffSurface = 1.0;
  constexpr double kOnSurface = 1e-3;
  constexpr double kTime = 1.7;
  constexpr double kQOverP = -0.4;
  const Vector3 direction = rawDirection.normalized();

  std::ostringstream errors;
  errors << std::scientific << std::setprecision(3);

  auto report = [&](const char* name, const auto& analytic,
                    const auto& numeric) {
    if (analytic.rows() != numeric.rows() ||
        analytic.cols() != numeric.cols()) {
      errors << name << " shape mismatch\n";
      return;
    }
    double maxAbs = 0.;
    std::ostringstream entries;
    int shown = 0;
    for (int row = 0; row < analytic.rows(); ++row) {
      for (int col = 0; col < analytic.cols(); ++col) {
        const double a = analytic(row, col);
        const double n = numeric(row, col);
        const double err = (std::isfinite(a) && std::isfinite(n))
                               ? std::abs(a - n)
                               : std::numeric_limits<double>::infinity();
        maxAbs = std::max(maxAbs, err);
        if (err > kTol && shown < 4) {
          entries << "    (" << row << "," << col << ") analytic=" << a
                  << " numeric=" << n << "\n";
          ++shown;
        }
      }
    }
    if (maxAbs > kTol) {
      errors << name << " max|diff|=" << maxAbs << "\n" << entries.str();
    }
  };

  auto perturb = [&](int axis, double delta, bool rotation) {
    Transform3 out = nominal;
    if (rotation) {
      out *= AngleAxis3(delta, Vector3::Unit(axis));
    } else {
      Vector3 deltaLocal = Vector3::Zero();
      deltaLocal[axis] = delta;
      out.translate(deltaLocal);
    }
    return out;
  };

  auto boundFromLocal = [&](const Vector2& local) {
    BoundVector bound = BoundVector::Zero();
    bound[eBoundLoc0] = local[0];
    bound[eBoundLoc1] = local[1];
    bound[eBoundPhi] = VectorHelpers::phi(direction);
    bound[eBoundTheta] = VectorHelpers::theta(direction);
    bound[eBoundQOverP] = kQOverP;
    bound[eBoundTime] = kTime;
    return bound;
  };

  auto boundOnRay = [&](const Surface& shifted, BoundVector& bound) {
    const auto intersection =
        shifted
            .intersect(gctx, position, direction, BoundaryTolerance::Infinite())
            .closest();
    if (!intersection.isValid()) {
      return false;
    }
    const Vector3 at = position + intersection.pathLength() * direction;
    auto local = shifted.globalToLocal(gctx, at, direction, kOnSurface);
    if (!local.ok()) {
      return false;
    }
    bound = boundFromLocal(local.value());
    return true;
  };

  const bool singularDirection =
      std::hypot(direction.x(), direction.y()) < 1e-8;

  FreeVector pathDerivative = FreeVector::Zero();
  pathDerivative.segment<3>(eFreePos0) = direction;
  const AlignmentToPathMatrix alignToPath =
      surface.alignmentToPathDerivative(gctx, position, direction);
  const FreeToBoundMatrix usedJac =
      surface.freeToBoundJacobian(gctx, position, direction);
  const AlignmentToBoundMatrix alignToBound =
      surface.alignmentToBoundDerivative(gctx, position, direction,
                                         pathDerivative);

  // ds/da
  {
    Eigen::Matrix<double, 1, eAlignmentSize> numeric =
        Eigen::Matrix<double, 1, eAlignmentSize>::Zero();
    bool ok = true;
    for (int axis = 0; axis < 3 && ok; ++axis) {
      for (const bool rotation : {false, true}) {
        const int col =
            (rotation ? eAlignmentRotation0 : eAlignmentCenter0) + axis;
        const auto plus = makeSurface(perturb(axis, kEps, rotation));
        const auto minus = makeSurface(perturb(axis, -kEps, rotation));
        const auto up = plus->intersect(gctx, position, direction,
                                        BoundaryTolerance::Infinite())
                            .closest();
        const auto um = minus
                            ->intersect(gctx, position, direction,
                                        BoundaryTolerance::Infinite())
                            .closest();
        if (!up.isValid() || !um.isValid()) {
          ok = false;
          break;
        }
        numeric(0, col) = (up.pathLength() - um.pathLength()) / (2. * kEps);
      }
    }
    if (!ok) {
      errors << "ds/da alignmentToPathDerivative evaluation failed\n";
    } else {
      report("ds/da alignmentToPathDerivative", alignToPath, numeric);
    }
  }

  // db/df. The bound angles are singular for a direction along z.
  if (singularDirection) {
    errors << "db/df freeToBoundJacobian direction is parallel to the z axis\n";
  } else {
    FreeVector nominalFree = FreeVector::Zero();
    nominalFree.segment<3>(eFreePos0) = position;
    nominalFree[eFreeTime] = kTime;
    nominalFree.segment<3>(eFreeDir0) = direction;
    nominalFree[eFreeQOverP] = kQOverP;
    FreeToBoundMatrix numeric = FreeToBoundMatrix::Zero();
    bool ok = true;
    for (int col = 0; col < static_cast<int>(eFreeSize) && ok; ++col) {
      FreeVector plus = nominalFree;
      FreeVector minus = nominalFree;
      plus[col] += kEps;
      minus[col] -= kEps;
      auto boundPlus =
          transformFreeToBoundParameters(plus, surface, gctx, kOffSurface);
      auto boundMinus =
          transformFreeToBoundParameters(minus, surface, gctx, kOffSurface);
      if (!boundPlus.ok() || !boundMinus.ok()) {
        ok = false;
        break;
      }
      numeric.col(col) = (boundPlus.value() - boundMinus.value()) / (2. * kEps);
    }
    if (!ok) {
      errors << "db/df freeToBoundJacobian evaluation failed\n";
    } else {
      report("db/df freeToBoundJacobian", usedJac, numeric);
    }
  }

  // db/da. The hit is the intersection of the fixed ray with the moved surface.
  {
    AlignmentToBoundMatrix numeric = AlignmentToBoundMatrix::Zero();
    bool ok = true;
    for (int axis = 0; axis < 3 && ok; ++axis) {
      for (const bool rotation : {false, true}) {
        const int col =
            (rotation ? eAlignmentRotation0 : eAlignmentCenter0) + axis;
        BoundVector plus = BoundVector::Zero();
        BoundVector minus = BoundVector::Zero();
        const auto movedPlus = makeSurface(perturb(axis, kEps, rotation));
        const auto movedMinus = makeSurface(perturb(axis, -kEps, rotation));
        if (!boundOnRay(*movedPlus, plus) || !boundOnRay(*movedMinus, minus)) {
          ok = false;
          break;
        }
        numeric.col(col) = (plus - minus) / (2. * kEps);
      }
    }
    if (!ok) {
      errors << "db/da alignmentToBoundDerivative evaluation failed\n";
    } else {
      report("db/da alignmentToBoundDerivative", alignToBound, numeric);
    }
  }

  BOOST_CHECK_MESSAGE(errors.str().empty(), errors.str());
}

}  // namespace ActsTests
