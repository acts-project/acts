// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/SpacePointUtility.hpp"

#include <iostream>
namespace Acts {

/// @brief Calculates (Delta theta)^2 + (Delta phi)^2 between two measurements
///
/// @param [in] pos1 position of the first measurement
/// @param [in] pos2 position the second measurement
/// @param [in] maxDistance Maximum distance between two measurements
/// @param [in] maxAngleTheta2 Maximum squared theta angle between two
/// measurements
/// @param [in] maxAnglePhi2 Maximum squared phi angle between two measurements
///
/// @return The squared sum within configuration parameters, otherwise -1
double SpacePointUtility::differenceOfMeasurementsChecked(
    const Vector3& pos1, const Vector3& pos2, const Vector3& posVertex,
    const double maxDistance, const double maxAngleTheta2,
    const double maxAnglePhi2) const {
  // Check if measurements are close enough to each other
  if ((pos1 - pos2).norm() > maxDistance) {
    return -1.;
  }

  // Calculate the angles of the vectors
  double phi1 = VectorHelpers::phi(pos1 - posVertex);
  double theta1 = VectorHelpers::theta(pos1 - posVertex);
  double phi2 = VectorHelpers::phi(pos2 - posVertex);
  double theta2 = VectorHelpers::theta(pos2 - posVertex);

  // Calculate the squared difference between the theta angles
  double diffTheta2 = (theta1 - theta2) * (theta1 - theta2);
  if (diffTheta2 > maxAngleTheta2) {
    return -1.;
  }
  // Calculate the squared difference between the phi angles
  double diffPhi2 = (phi1 - phi2) * (phi1 - phi2);
  if (diffPhi2 > maxAnglePhi2) {
    return -1.;
  }
  // Return the squared distance between both vector
  return diffTheta2 + diffPhi2;
}

std::pair<Vector3, Vector2> SpacePointUtility::globalCoords(
    const GeometryContext& gctx, const Measurement& meas) const {
  const auto& slink =
      std::visit([](const auto& x) { return &x.sourceLink(); }, meas);

  const auto geoId = slink->geometryId();

  const Surface* surface = m_config.trackingGeometry->findSurface(geoId);
  auto [localPos, localCov] = std::visit(
      [](const auto& measurement) {
        auto expander = measurement.expander();
        BoundVector par = expander * measurement.parameters();
        BoundSymMatrix cov =
            expander * measurement.covariance() * expander.transpose();
        // extract local position
        Vector2 lpar(par[eBoundLoc0], par[eBoundLoc1]);
        // extract local position covariance.
        SymMatrix2 lcov = cov.block<2, 2>(eBoundLoc0, eBoundLoc0);
        return std::make_pair(lpar, lcov);
      },
      meas);
  Vector3 globalPos = surface->localToGlobal(gctx, localPos, Vector3());
  RotationMatrix3 rotLocalToGlobal =
      surface->referenceFrame(gctx, globalPos, Vector3());

  // the space point requires only the variance of the transverse and
  // longitudinal position. reduce computations by transforming the
  // covariance directly from local to rho/z.
  //
  // compute Jacobian from global coordinates to rho/z
  //
  //         rho = sqrt(x² + y²)
  // drho/d{x,y} = (1 / sqrt(x² + y²)) * 2 * {x,y}
  //             = 2 * {x,y} / r
  //       dz/dz = 1 (duuh!)
  //
  auto x = globalPos[ePos0];
  auto y = globalPos[ePos1];
  auto scale = 2 / std::hypot(x, y);
  ActsMatrix<2, 3> jacXyzToRhoZ = ActsMatrix<2, 3>::Zero();
  jacXyzToRhoZ(0, ePos0) = scale * x;
  jacXyzToRhoZ(0, ePos1) = scale * y;
  jacXyzToRhoZ(1, ePos2) = 1;
  // compute Jacobian from local coordinates to rho/z
  ActsMatrix<2, 2> jac =
      jacXyzToRhoZ * rotLocalToGlobal.block<3, 2>(ePos0, ePos0);
  // compute rho/z variance
  ActsVector<2> var = (jac * localCov * jac.transpose()).diagonal();

  auto gcov = Vector2(var[0], var[1]);
  return std::make_pair(globalPos, gcov);
}

Vector2 SpacePointUtility::calcGlobalVars(const GeometryContext& gctx,
                                          const Measurement& measFront,
                                          const Measurement& measBack,
                                          const Vector3& globalPos,
                                          const double theta) const {
  const auto var1 = getLoc0Var(measFront);
  const auto var2 = getLoc0Var(measBack);
  // strip1 and strip2 are tilted at +/- theta/2

  double sigma_x = std::hypot(var1, var2) / (2 * sin(theta * 0.5));
  double sigma_y = std::hypot(var1, var2) / (2 * cos(theta * 0.5));

  // projection to the surface with strip1.
  double sig_x1 = sigma_x * cos(0.5 * theta) + sigma_y * sin(0.5 * theta);
  double sig_y1 = sigma_y * cos(0.5 * theta) + sigma_x * sin(0.5 * theta);
  SymMatrix2 lcov;
  lcov << sig_x1, 0, 0, sig_y1;

  const auto& slink_meas1 =
      std::visit([](const auto& x) { return &x.sourceLink(); }, measFront);

  const auto geoId = slink_meas1->geometryId();

  auto gcov = globalCov(gctx, geoId, globalPos, lcov);

  return gcov;
}

double SpacePointUtility::getLoc0Var(const Measurement& meas) const {
  auto cov = std::visit(
      [](const auto& x) {
        auto expander = x.expander();
        BoundSymMatrix bcov = expander * x.covariance() * expander.transpose();
        SymMatrix2 lcov = bcov.block<2, 2>(eBoundLoc0, eBoundLoc0);
        return lcov;
      },
      meas);
  return cov(0, 0);
}

Vector2 SpacePointUtility::globalCov(const GeometryContext& gctx,
                                     const GeometryIdentifier& geoId,
                                     const Vector3& globalPos,
                                     const SymMatrix2& localCov) const {
  Vector3 globalFakeMom(1, 1, 1);

  const Surface* surface = m_config.trackingGeometry->findSurface(geoId);

  RotationMatrix3 rotLocalToGlobal =
      surface->referenceFrame(gctx, globalPos, globalFakeMom);

  auto x = globalPos[ePos0];
  auto y = globalPos[ePos1];
  auto scale = 2 / std::hypot(x, y);
  ActsMatrix<2, 3> jacXyzToRhoZ = ActsMatrix<2, 3>::Zero();
  jacXyzToRhoZ(0, ePos0) = scale * x;
  jacXyzToRhoZ(0, ePos1) = scale * y;
  jacXyzToRhoZ(1, ePos2) = 1;
  // compute Jacobian from local coordinates to rho/z
  ActsMatrix<2, 2> jac =
      jacXyzToRhoZ * rotLocalToGlobal.block<3, 2>(ePos0, ePos0);
  // compute rho/z variance
  ActsVector<2> var = (jac * localCov * jac.transpose()).diagonal();

  auto gcov = Vector2(var[0], var[1]);

  return gcov;
}

/// @brief This function performs a straight forward calculation of a space
/// point and returns whether it was succesful or not.
///
/// @param [in] stripEnds1 Top and bottom end of the first strip detector
/// element
/// @param [in] stripEnds1 Top and bottom end of the second strip detector
/// element
/// @param [in] posVertex Position of the vertex
/// @param [in, out] spParams Data container of the calculations
/// @param [in] stripLengthTolerance Tolerance scaling factor on the strip
/// detector element length
///
/// @return Boolean statement whether the space point calculation was succesful
bool SpacePointUtility::calculateStripSPPosition(
    const std::pair<Vector3, Vector3>& stripEnds1,
    const std::pair<Vector3, Vector3>& stripEnds2, const Vector3& posVertex,
    SpacePointParameters& spParams, const double stripLengthTolerance) const {
  /// The following algorithm is meant for finding the position on the first
  /// strip if there is a corresponding Measurement on the second strip. The
  /// resulting point is a point x on the first surfaces. This point is
  /// along a line between the points a (top end of the strip)
  /// and b (bottom end of the strip). The location can be parametrized as
  /// 	2 * x = (1 + m) a + (1 - m) b
  /// as function of the scalar m. m is a parameter in the interval
  /// -1 < m < 1 since the hit was on the strip. Furthermore, the vector
  /// from the vertex to the Measurement on the second strip y is needed to be a
  /// multiple k of the vector from vertex to the hit on the first strip x.
  /// As a consequence of this demand y = k * x needs to be on the
  /// connecting line between the top (c) and bottom (d) end of
  /// the second strip. If both measurements correspond to each other, the
  /// condition
  /// 	y * (c X d) = k * x (c X d) = 0 ("X" represents a cross product)
  /// needs to be fulfilled. Inserting the first equation into this
  /// equation leads to the condition for m as given in the following
  /// algorithm and therefore to the calculation of x.
  /// The same calculation can be repeated for y. Its corresponding
  /// parameter will be named n.

  spParams.q = stripEnds1.first - stripEnds1.second;
  spParams.r = stripEnds2.first - stripEnds2.second;
  spParams.s = stripEnds1.first + stripEnds1.second - 2 * posVertex;
  spParams.t = stripEnds2.first + stripEnds2.second - 2 * posVertex;
  spParams.qs = spParams.q.cross(spParams.s);
  spParams.rt = spParams.r.cross(spParams.t);
  spParams.m = -spParams.s.dot(spParams.rt) / spParams.q.dot(spParams.rt);

  // Set the limit for the parameter
  if (spParams.limit == 1. && stripLengthTolerance != 0.) {
    spParams.limit = 1. + stripLengthTolerance;
  }

  // Check if m and n can be resolved in the interval (-1, 1)
  return (fabs(spParams.m) <= spParams.limit &&
          fabs(spParams.n = -spParams.t.dot(spParams.qs) /
                            spParams.r.dot(spParams.qs)) <= spParams.limit);
}

/// @brief This function tests if a space point can be estimated by a more
/// tolerant treatment of construction. In fact, this function indirectly
/// allows shifts of the vertex.
///
/// @param [in] spParams container that stores geometric parameters and rules of
/// the space point formation
/// @param [in] stripLengthGapTolerance Tolerance scaling factor of the gap
/// between strip detector elements
///
/// @return indicator if the test was successful
bool SpacePointUtility::recoverSpacePoint(
    SpacePointParameters& spParams, double stripLengthGapTolerance) const {
  /// Consider some cases that would allow an easy exit
  // Check if the limits are allowed to be increased
  if (stripLengthGapTolerance <= 0.) {
    return false;
  }
  spParams.qmag = spParams.q.norm();
  // Increase the limits. This allows a check if the point is just slightly
  // outside the SDE
  spParams.limitExtended =
      spParams.limit + stripLengthGapTolerance / spParams.qmag;
  // Check if m is just slightly outside
  if (fabs(spParams.m) > spParams.limitExtended) {
    return false;
  }
  // Calculate n if not performed previously
  if (spParams.n == 0.) {
    spParams.n = -spParams.t.dot(spParams.qs) / spParams.r.dot(spParams.qs);
  }
  // Check if n is just slightly outside
  if (fabs(spParams.n) > spParams.limitExtended) {
    return false;
  }

  /// The following code considers an overshoot of m and n in the same direction
  /// of their SDE. The term "overshoot" represents the amount of m or n outside
  /// its regular interval (-1, 1).
  /// It calculates which overshoot is worse. In order to compare both, the
  /// overshoot in n is projected onto the first surface by considering the
  /// normalized projection of r onto q.
  /// This allows a rescaling of the overshoot. The worse overshoot will be set
  /// to +/-1, the parameter with less overshoot will be moved towards 0 by the
  /// worse overshoot.
  /// In order to treat both SDEs equally, the rescaling eventually needs to be
  /// performed several times. If these shifts allows m and n to be in the
  /// limits, the space point can be stored.
  /// @note This shift can be understood as a shift of the particle's
  /// trajectory. This is leads to a shift of the vertex. Since these two points
  /// are treated independently from other measurement, it is also possible to
  /// consider this as a change in the slope of the particle's trajectory.
  ///  The would also move the vertex position.

  // Calculate the scaling factor to project lengths of the second SDE on the
  // first SDE
  double secOnFirstScale =
      spParams.q.dot(spParams.r) / (spParams.qmag * spParams.qmag);
  // Check if both overshoots are in the same direction
  if (spParams.m > 1. && spParams.n > 1.) {
    // Calculate the overshoots
    double mOvershoot = spParams.m - 1.;
    double nOvershoot =
        (spParams.n - 1.) * secOnFirstScale;  // Perform projection
    // Resolve worse overshoot
    double biggerOvershoot = std::max(mOvershoot, nOvershoot);
    // Move m and n towards 0
    spParams.m -= biggerOvershoot;
    spParams.n -= (biggerOvershoot / secOnFirstScale);
    // Check if this recovered the space point
    return fabs(spParams.m) < spParams.limit &&
           fabs(spParams.n) < spParams.limit;
  }
  // Check if both overshoots are in the same direction
  if (spParams.m < -1. && spParams.n < -1.) {
    // Calculate the overshoots
    double mOvershoot = -(spParams.m + 1.);
    double nOvershoot =
        -(spParams.n + 1.) * secOnFirstScale;  // Perform projection
    // Resolve worse overshoot
    double biggerOvershoot = std::max(mOvershoot, nOvershoot);
    // Move m and n towards 0
    spParams.m += biggerOvershoot;
    spParams.n += (biggerOvershoot / secOnFirstScale);
    // Check if this recovered the space point
    return fabs(spParams.m) < spParams.limit &&
           fabs(spParams.n) < spParams.limit;
  }
  // No solution could be found
  return false;
}

/// @brief Calculates a space point whithout using the vertex
/// @note This is mostly to resolve space points from cosmic data
/// @param stripEnds1 The ends of one strip
/// @param stripEnds2 The ends of another strip
/// @param spParams SpacePointParamaters for the SP
/// @return parameter that indicates the location of the space point; returns
/// 1. if it failed
/// @note The meaning of the parameter is explained in more detail in the
/// function body
double SpacePointUtility::calcPerpendicularProjection(
    const std::pair<Vector3, Vector3>& stripEnds1,
    const std::pair<Vector3, Vector3>& stripEnds2,
    SpacePointParameters& spParams) const {
  /// This approach assumes that no vertex is available. This option aims to
  /// approximate the space points from cosmic data.
  /// The underlying assumption is that the best point is given by the  closest
  /// distance between both lines describing the SDEs.
  /// The point x on the first SDE is parametrized as a + lambda0 * q with  the
  /// top end a of the strip and the vector q = a - b(ottom end of the  strip).
  /// An analogous parametrization is performed of the second SDE with y = c  +
  /// lambda1 * r.
  /// x get resolved by resolving lambda0 from the condition that |x-y| is  the
  /// shortest distance between two skew lines.

  spParams.q = stripEnds1.first - stripEnds1.second;
  spParams.r = stripEnds2.first - stripEnds2.second;

  Vector3 ac = stripEnds2.first - stripEnds1.first;
  double qr = (spParams.q).dot(spParams.r);
  double denom = spParams.q.dot(spParams.q) - qr * qr;

  // Check for numerical stability
  if (fabs(denom) > 10e-7) {
    // Return lambda0
    return (ac.dot(spParams.r) * qr -
            ac.dot(spParams.q) * (spParams.r).dot(spParams.r)) /
           denom;
  }
  // lambda0 is in the interval [-1,0]. This return serves as error check.
  return 1.;
}

}  // namespace Acts
