// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/CovarianceEngine.hpp"

#include "Acts/Surfaces/Surface.hpp"

namespace Acts {
namespace {

/// Projection Jacobian from bound parameters on the surface to free parameters.
BoundToFreeMatrix boundToFreeJacobian(const StepperState& state,
                                      const Surface& surface) {
  // convert free to local parameters position
  Vector2D loc = Vector2D::Zero();
  surface.globalToLocal(state.geoContext, state.pos, state.dir, loc);
  BoundVector pars;
  pars[eBoundLoc0] = loc[0];
  pars[eBoundLoc1] = loc[1];
  pars[eBoundTime] = state.t;
  pars[eBoundPhi] = VectorHelpers::phi(state.dir);
  pars[eBoundTheta] = VectorHelpers::theta(state.dir);
  pars[eBoundQOverP] = (state.q != 0) ? (state.q / state.p) : (1 / state.p);
  BoundToFreeMatrix jacToGlobal;
  surface.initJacobianToGlobal(state.geoContext, jacToGlobal, state.pos,
                               state.dir, pars);
  return jacToGlobal;
}

/// Projection Jacobian from curvilinear parameters to free parameters.
///
/// In principle, this should be equivalent to the projection Jacobian for
/// general bound parameters with an appropriate Curvilinear surface. Use an
/// explicit implementation for faster computation.
BoundToFreeMatrix curvilinearToFreeJacobian(const Vector3D& unitDirection) {
  // normalized direction vector components expressed as phi/theta
  const double x = unitDirection[0];  // == cos(phi) * sin(theta)
  const double y = unitDirection[1];  // == sin(phi) * sin(theta)
  const double z = unitDirection[2];  // == cos(theta)
  // extract phi/theta expressions
  const auto cosTheta = z;
  const auto sinTheta = std::hypot(x, y);
  const auto invSinTheta = 1 / sinTheta;
  const auto cosPhi = x * invSinTheta;
  const auto sinPhi = y * invSinTheta;
  BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();
  jacToGlobal(eFreePos0, eBoundLoc0) = -sinPhi;
  jacToGlobal(eFreePos0, eBoundLoc1) = -cosPhi * cosTheta;
  jacToGlobal(eFreePos1, eBoundLoc0) = cosPhi;
  jacToGlobal(eFreePos1, eBoundLoc1) = -sinPhi * cosTheta;
  jacToGlobal(eFreePos2, eBoundLoc1) = sinTheta;
  jacToGlobal(eFreeTime, eBoundTime) = 1;
  jacToGlobal(eFreeDir0, eBoundPhi) = -sinTheta * sinPhi;
  jacToGlobal(eFreeDir0, eBoundTheta) = cosTheta * cosPhi;
  jacToGlobal(eFreeDir1, eBoundPhi) = sinTheta * cosPhi;
  jacToGlobal(eFreeDir1, eBoundTheta) = cosTheta * sinPhi;
  jacToGlobal(eFreeDir2, eBoundTheta) = -sinTheta;
  jacToGlobal(eFreeQOverP, eBoundQOverP) = 1;
  return jacToGlobal;
}

/// Projection Jacobian from free to curvilinear parameters.
FreeToBoundMatrix freeToCurvilinearJacobian(const Vector3D& unitDirection) {
  // normalized direction vector components expressed as phi/theta
  const double x = unitDirection[0];  // == cos(phi) * sin(theta)
  const double y = unitDirection[1];  // == sin(phi) * sin(theta)
  const double z = unitDirection[2];  // == cos(theta)
  // extract phi/theta expressions
  const auto cosTheta = z;
  const auto sinTheta = std::hypot(x, y);
  const auto invSinTheta = 1 / sinTheta;
  const auto cosPhi = x * invSinTheta;
  const auto sinPhi = y * invSinTheta;
  FreeToBoundMatrix jacToCurv = FreeToBoundMatrix::Zero();
  if (std::abs(cosTheta) < s_curvilinearProjTolerance) {
    // We normally operate in curvilinear coordinates defined as above
    jacToCurv(eBoundLoc0, eFreePos0) = -sinPhi;
    jacToCurv(eBoundLoc0, eFreePos1) = cosPhi;
    jacToCurv(eBoundLoc1, eFreePos0) = -cosPhi * cosTheta;
    jacToCurv(eBoundLoc1, eFreePos1) = -sinPhi * cosTheta;
    jacToCurv(eBoundLoc1, eFreePos2) = sinTheta;
  } else {
    // Under grazing incidence to z, the above coordinate system definition
    // becomes numerically unstable, and we need to switch to another one
    const auto c = std::hypot(y, z);
    const auto invC = 1 / c;
    jacToCurv(eBoundLoc0, eFreePos1) = -z * invC;
    jacToCurv(eBoundLoc0, eFreePos2) = y * invC;
    jacToCurv(eBoundLoc1, eFreePos0) = c;
    jacToCurv(eBoundLoc1, eFreePos1) = -x * y * invC;
    jacToCurv(eBoundLoc1, eFreePos2) = -x * z * invC;
  }
  // Time parameter
  jacToCurv(eT, eFreeTime) = 1;
  // Directional and momentum parameters for curvilinear
  // TODO this becomes unstable for direction along z
  //      sin(theta) -> 0 -> 1/sin(theta) -> inf
  jacToCurv(eBoundPhi, eFreeDir0) = -sinPhi * invSinTheta;
  jacToCurv(eBoundPhi, eFreeDir1) = cosPhi * invSinTheta;
  jacToCurv(eBoundTheta, eFreeDir2) = -invSinTheta;
  jacToCurv(eBoundQOverP, eFreeQOverP) = 1;
  return jacToCurv;
}

/// Apply Jacobian modifications from projection onto the surface.
///
/// This function treats the modifications of the jacobian related to
/// the projection onto a surface. Since a variation of the start parameters
/// within a given uncertainty would lead to a variation of the end
/// parameters, these need to be propagated onto the target surface. This an
/// approximated approach to treat the (assumed) small change.
///
/// @param [in] position Particle position
/// @param [in] unitDirection Particle direction
/// @param [in] derivative Stepper propagation derivative
/// @param [in] surface Local surface to project onto
/// @param [in] referenceFrameT Transpose local reference frame on surface
/// @param [in] geoContext Geometry context
/// @param [in,out] jacToGlobal curvilinear-to-free Jacobian to modify
void applyDerivativeCorrectionBound(const Vector3D& position,
                                    const Vector3D& unitDirection,
                                    const FreeVector& derivative,
                                    const Surface& surface,
                                    const RotationMatrix3D& referenceFrameT,
                                    const GeometryContext& geoContext,
                                    BoundToFreeMatrix& jacToGlobal) {
  // Calculate the form factors for the derivatives
  const BoundRowVector sfactors = surface.derivativeFactors(
      geoContext, position, unitDirection, referenceFrameT, jacToGlobal);
  jacToGlobal -= derivative * sfactors;
}

/// Apply Jacobian modifications from projection onto the curvilinear frame.
///
/// @param [in] unitDirection Particle direction
/// @param [in] derivative Stepper propagation derivative
/// @param [in,out] jacToGlobal curvilinear-to-free Jacobian to modify
void applyDerivativeCorrectionCurvilinear(const Vector3D& unitDirection,
                                          const FreeVector& derivative,
                                          BoundToFreeMatrix& jacToGlobal) {
  // Transport the covariance
  const BoundRowVector sfactors =
      unitDirection.transpose() *
      jacToGlobal.template topLeftCorner<3, eBoundParametersSize>();
  jacToGlobal -= derivative * sfactors;
}

}  // namespace
namespace detail {

std::tuple<BoundParameters, BoundMatrix, double> boundState(
    StepperState& state, const Surface& surface) {
  // Transport the covariance to here
  std::optional<BoundSymMatrix> cov = std::nullopt;
  if (state.covTransport) {
    // Initialize the transport final frame jacobian
    transportCovarianceToBound(state, surface);
    cov = state.cov;
  }
  // Create the bound parameters
  BoundParameters parameters(state.geoContext, std::move(cov), state.pos,
                             state.p * state.dir, state.q, state.t,
                             surface.getSharedPtr());
  // Create the bound state
  return std::make_tuple(std::move(parameters), state.jacobian,
                         state.pathAccumulated);
}

std::tuple<CurvilinearParameters, BoundMatrix, double> curvilinearState(
    StepperState& state) {
  // Transport the covariance to here
  std::optional<BoundSymMatrix> cov = std::nullopt;
  if (state.covTransport) {
    transportCovarianceToCurvilinear(state);
    cov = state.cov;
  }
  // Create the curvilinear parameters
  CurvilinearParameters parameters(std::move(cov), state.pos,
                                   state.p * state.dir, state.q, state.t);
  // Create the bound state
  return std::make_tuple(std::move(parameters), state.jacobian,
                         state.pathAccumulated);
}

void transportCovarianceToBound(StepperState& state, const Surface& surface) {
  // Jacobian from current free parameters onto the surface
  FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
  const auto& referenceFrameTranspose = surface.initJacobianToLocal(
      state.geoContext, jacToLocal, state.pos, state.dir);
  // Add corrections due to propagation derivatives
  // TODO 20200319 msmk: why does this change the initial projection and not the
  //                     final one?
  applyDerivativeCorrectionBound(state.pos, state.dir, state.derivative,
                                 surface, referenceFrameTranspose,
                                 state.geoContext, state.jacToGlobal);
  // combine Jacobian for bound parameters from:
  // 1. projection from the initial bound frame to free parameters
  // 2. free transport Jacobian along trajectory as computed by stepper
  // 3. Jacobian onto the curvilinear frame at the current position
  // matrix multiplication uses right-to-left order, i.e. 3 * 2 * 1
  state.jacobian = jacToLocal * state.jacTransport * state.jacToGlobal;

  // transport the covariance matrix
  state.cov = state.jacobian * state.cov * state.jacobian.transpose();

  // reset Jacobians for a propagation starting on the surface
  state.jacToGlobal = boundToFreeJacobian(state, surface);
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
}

void transportCovarianceToCurvilinear(StepperState& state) {
  // Jacobian from current free parameters onto the curvilinear frame
  const FreeToBoundMatrix jacToLocal = freeToCurvilinearJacobian(state.dir);
  // Add corrections due to propagation derivatives
  // TODO 20200319 msmk: why does this change the initial projection and not the
  //                     final one?
  applyDerivativeCorrectionCurvilinear(state.dir, state.derivative,
                                       state.jacToGlobal);
  // see transportCovarianceToSurface for explanation
  state.jacobian = jacToLocal * state.jacTransport * state.jacToGlobal;

  // transport the covariance matrix
  state.cov = state.jacobian * state.cov * state.jacobian.transpose();

  // reset Jacobians for a propagation starting on the curvilinear frame
  state.jacToGlobal = curvilinearToFreeJacobian(state.dir);
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
}

}  // namespace detail
}  // namespace Acts
