// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <format>

namespace Acts {
template <Experimental::CompositeSpacePoint SpacePoint_t>
std::string toString(const SpacePoint_t& measurement) {
  if constexpr (hasPrintOperator<SpacePoint_t>) {
    std::ostringstream sstr{};
    sstr << measurement;
    return sstr.str();
  } else {
    using ResidualIdx =
        Experimental::detail::CompSpacePointAuxiliaries::ResidualIdx;
    if (measurement.isStraw()) {
      return std::format(
          "straw SP @ {:} with r: {:.3f}+-{:.3f} & wire : {:} ",
          toString(measurement.localPosition()), measurement.driftRadius(),
          std::sqrt(
              measurement.covariance()[toUnderlying(ResidualIdx::bending)]),
          toString(measurement.sensorDirection()));
    } else {
      return std::format(
          "strip SP @ {:} with normal: {:}, strip dir: {:}, to next {:}, "
          "measures loc0/loc1/time: {:}/{:}/{:}",
          toString(measurement.localPosition()),
          toString(measurement.planeNormal()),
          toString(measurement.sensorDirection()),
          toString(measurement.toNextSensor()),
          measurement.measuresLoc0() ? "yes" : "no",
          measurement.measuresLoc1() ? "yes" : "no",
          measurement.hasTime() ? "yes" : "no");
    }
  }
}
}  // namespace Acts

namespace Acts::Experimental::detail {

template <CompositeSpacePoint SpacePoint_t>
CompSpacePointAuxiliaries::Vector CompSpacePointAuxiliaries::extrapolateToPlane(
    const Line_t& line, const SpacePoint_t& hit) {
  return extrapolateToPlane(line.position(), line.direction(), hit);
}
template <CompositeSpacePoint SpacePoint_t>
CompSpacePointAuxiliaries::Vector CompSpacePointAuxiliaries::extrapolateToPlane(
    const Vector& pos, const Vector& dir, const SpacePoint_t& hit) {
  using namespace Acts::PlanarHelper;
  const auto planeIsect =
      intersectPlane(pos, dir, hit.planeNormal(), hit.localPosition());
  return planeIsect.position();
}
template <CompositeSpacePoint SpacePoint_t>
double CompSpacePointAuxiliaries::chi2Term(const Line_t& line,
                                           const SpacePoint_t& hit) {
  return chi2Term(line.position(), line.direction(), hit);
}
template <CompositeSpacePoint SpacePoint_t>
double CompSpacePointAuxiliaries::chi2Term(const Vector& pos, const Vector& dir,
                                           const SpacePoint_t& hit) {
  double chiSq{0.};
  using namespace Acts::detail::LineHelper;
  constexpr auto bendIdx = toUnderlying(ResidualIdx::bending);
  constexpr auto nonBendIdx = toUnderlying(ResidualIdx::nonBending);
  if (hit.isStraw()) {
    const double dist = Acts::abs(
        signedDistance(pos, dir, hit.localPosition(), hit.sensorDirection()));
    chiSq = Acts::square(dist - hit.driftRadius()) / hit.covariance()[bendIdx];
    if (hit.measuresLoc0()) {
      auto closePointOnStraw =
          lineIntersect(pos, dir, hit.localPosition(), hit.sensorDirection());
      chiSq += Acts::square(closePointOnStraw.pathLength()) /
               hit.covariance()[nonBendIdx];
    }

  } else {
    const Vector distOnPlane =
        (extrapolateToPlane(pos, dir, hit) - hit.localPosition());
    const Vector& b1{hit.toNextSensor()};
    const Vector& b2{hit.sensorDirection()};
    Vector2 dist{distOnPlane.dot(b1), distOnPlane.dot(b2)};
    if (hit.measuresLoc0() && hit.measuresLoc1()) {
      /// Check whether the two vectors are orthogonal
      constexpr double tolerance = 1.e-8;
      const double stripAngle = b1.dot(b2);
      if (stripAngle > tolerance) {
        const double invDist = 1. / (1. - square(stripAngle));
        SquareMatrix<2> stereoDecomp{invDist * SquareMatrix<2>::Identity()};
        stereoDecomp(1, 0) = stereoDecomp(0, 1) = -stripAngle * invDist;
        dist = stereoDecomp * dist;
      }
      chiSq = Acts::square(dist[0]) / hit.covariance()[bendIdx] +
              Acts::square(dist[1]) / hit.covariance()[nonBendIdx];
    } else if (hit.measuresLoc0()) {
      chiSq = Acts::square(dist[0]) / hit.covariance()[nonBendIdx];
    } else {
      chiSq = Acts::square(dist[0]) / hit.covariance()[bendIdx];
      if (hit.covariance()[nonBendIdx] >
          std::numeric_limits<double>::epsilon()) {
        chiSq += square(dist[1]) / hit.covariance()[nonBendIdx];
      }
    }
  }
  return chiSq;
}
template <CompositeSpacePoint SpacePoint_t>
double CompSpacePointAuxiliaries::chi2Term(
    const Line_t& line, const Acts::Transform3& localToGlobal, const double t0,
    const SpacePoint_t& hit) {
  return chi2Term(line.position(), line.direction(), localToGlobal, t0, hit);
}
template <CompositeSpacePoint SpacePoint_t>
double CompSpacePointAuxiliaries::chi2Term(const Vector& pos, const Vector& dir,
                                           const Transform3& localToGlobal,
                                           const double t0,
                                           const SpacePoint_t& hit) {
  return chi2Term(
      pos, dir,
      t0 + (localToGlobal * extrapolateToPlane(pos, dir, hit)).norm() /
               PhysicalConstants::c,
      hit);
}

template <CompositeSpacePoint SpacePoint_t>
double CompSpacePointAuxiliaries::chi2Term(const Line_t& line, const double t0,
                                           const SpacePoint_t& hit) {
  return chi2Term(line.position(), line.direction(), t0, hit);
}
template <CompositeSpacePoint SpacePoint_t>
double CompSpacePointAuxiliaries::chi2Term(const Vector& pos, const Vector& dir,
                                           const double t0,
                                           const SpacePoint_t& hit) {
  double chiSq = chi2Term(pos, dir, hit);
  if (!hit.hasTime() || hit.isStraw()) {
    return chiSq;
  }
  chiSq += Acts::square(hit.time() - t0) /
           hit.covariance()[toUnderlying(ResidualIdx::time)];
  return chiSq;
}

template <CompositeSpacePoint Point_t>
int CompSpacePointAuxiliaries::strawSign(const Line_t& line,
                                         const Point_t& strawSp) {
  return strawSign(line.position(), line.direction(), strawSp);
}

template <CompositeSpacePoint Point_t>
int CompSpacePointAuxiliaries::strawSign(const Vector& pos, const Vector& dir,
                                         const Point_t& strawSp) {
  if (!strawSp.isStraw()) {
    return 0;
  }
  const double dist = Acts::detail::LineHelper::signedDistance(
      pos, dir, strawSp.localPosition(), strawSp.sensorDirection());
  return copySign(1, dist);
}
template <CompositeSpacePointContainer StrawCont_t>
std::vector<int> CompSpacePointAuxiliaries::strawSigns(
    const Line_t& line, const StrawCont_t& measurements) {
  return strawSigns(line.position(), line.direction(), measurements);
}
template <CompositeSpacePointContainer StrawCont_t>
std::vector<int> CompSpacePointAuxiliaries::strawSigns(
    const Vector& pos, const Vector& dir, const StrawCont_t& measurements) {
  std::vector<int> signs{};
  signs.reserve(measurements.size());
  for (const auto& strawSp : measurements) {
    signs.push_back(strawSign(pos, dir, *strawSp));
  }
  return signs;
}
template <CompositeSpacePoint Point_t>
void CompSpacePointAuxiliaries::updateSpatialResidual(
    const Line_t& line, const Point_t& spacePoint) {
  ACTS_DEBUG(__func__ << "() " << __LINE__ << " Update residual of "
                      << toString(line.position()) << " + "
                      << toString(line.direction()) << " w.r.t\n"
                      << toString(spacePoint));
  if (spacePoint.isStraw()) {
    /// Fetch the hit position & direction
    const auto& wireDir{spacePoint.sensorDirection()};
    /// Calculate the distance from the two reference points
    const Vector hitMinSeg = spacePoint.localPosition() - line.position();

    if (!updateStrawResidual(line, hitMinSeg, wireDir,
                             spacePoint.driftRadius())) {
      return;
    }

    if (m_cfg.calcAlongStraw && spacePoint.measuresLoc0()) {
      /// If the tube is a twin-tube, the hit position is no longer arbitrary
      /// along the wire. Calculate the distance along the wire towards the
      /// point of closest approach.
      updateAlongTheStraw(line, hitMinSeg, wireDir);
    }
  } else {
    updateStripResidual(line, spacePoint.planeNormal(),
                        spacePoint.toNextSensor(), spacePoint.sensorDirection(),
                        spacePoint.localPosition(), spacePoint.measuresLoc1(),
                        spacePoint.measuresLoc0());
  }
}

template <CompositeSpacePoint Point_t>
void CompSpacePointAuxiliaries::updateFullResidual(const Line_t& line,
                                                   const double timeOffset,
                                                   const Point_t& spacePoint,
                                                   const double driftV,
                                                   const double driftA) {
  /// Calculate first the spatial residual
  updateSpatialResidual(line, spacePoint);

  /// Calculate the time residual for strip-like measurements
  if (!spacePoint.isStraw()) {
    /// If the measurement does not provide time, then simply reset the time
    /// partial components
    if (!spacePoint.hasTime()) {
      resetTime();
      return;
    }
    updateTimeStripRes(spacePoint.toNextSensor(), spacePoint.sensorDirection(),
                       spacePoint.localPosition(), spacePoint.measuresLoc1(),
                       spacePoint.time(), timeOffset);
  } else {
    updateTimeStrawRes(line, spacePoint.localPosition() - line.position(),
                       spacePoint.sensorDirection(), spacePoint.driftRadius(),
                       driftV, driftA);
  }
}

}  // namespace Acts::Experimental::detail
