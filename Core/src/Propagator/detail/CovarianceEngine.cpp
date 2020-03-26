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
namespace detail {

BoundParameters makeBoundParameters(const FreeVector& freeParams,
                                    const BoundSymMatrix& boundCov,
                                    bool covIsValid, const Surface& surface,
                                    const GeometryContext& geoCtx) {
  const auto& pos = freeParams.segment<3>(eFreePos0);
  const auto& dir = freeParams.segment<3>(eFreeDir0);

  // convert free to local position
  Vector2D loc = Vector2D::Zero();
  surface.globalToLocal(geoCtx, pos, dir, loc);
  // convert to bound parameters
  BoundVector boundParams;
  boundParams[eBoundLoc0] = loc[0];
  boundParams[eBoundLoc1] = loc[1];
  boundParams[eBoundTime] = freeParams[eFreeTime];
  boundParams[eBoundPhi] = VectorHelpers::phi(dir);
  boundParams[eBoundTheta] = VectorHelpers::theta(dir);
  boundParams[eBoundQOverP] = freeParams[eFreeQOverP];
  // construct the optional covariance
  auto cov = covIsValid ? std::make_optional(boundCov) : std::nullopt;
  // create the bound parameters
  return {geoCtx, std::move(cov), std::move(boundParams),
          surface.getSharedPtr()};
}

CurvilinearParameters makeCurvilinearParameters(
    const FreeVector& freeParams, const BoundSymMatrix& curvilinearCov,
    bool covIsValid) {
  const auto& pos = freeParams.segment<3>(eFreePos0);
  const auto& dir = freeParams.segment<3>(eFreeDir0);

  // extract momentum and charge
  // TODO use q/p directly without conversion. requires track parameter changes
  const auto p = std::abs(1 / freeParams[eFreeQOverP]);
  const auto q = std::copysign(1, freeParams[eFreeQOverP]);
  // extract time
  const auto t = freeParams[eFreeTime];

  // construct the optional covariance
  auto cov = covIsValid ? std::make_optional(curvilinearCov) : std::nullopt;
  // create the curvilinear parameters
  return {std::move(cov), pos, p * dir, q, t};
}

namespace {

/// Projection Jacobian from bound parameters on the surface to free parameters.
BoundToFreeMatrix boundToFreeJacobian(const FreeVector& freeParams,
                                      const Surface& surface,
                                      const GeometryContext& geoCtx) {
  const auto& pos = freeParams.segment<3>(eFreePos0);
  const auto& dir = freeParams.segment<3>(eFreeDir0);

  // convert free to local parameters position
  Vector2D loc = Vector2D::Zero();
  surface.globalToLocal(geoCtx, pos, dir, loc);
  BoundVector boundParams;
  boundParams[eBoundLoc0] = loc[0];
  boundParams[eBoundLoc1] = loc[1];
  boundParams[eBoundTime] = freeParams[eFreeTime];
  boundParams[eBoundPhi] = VectorHelpers::phi(dir);
  boundParams[eBoundTheta] = VectorHelpers::theta(dir);
  boundParams[eBoundQOverP] = freeParams[eFreeQOverP];
  BoundToFreeMatrix jacToGlobal;
  surface.initJacobianToGlobal(geoCtx, jacToGlobal, pos, dir, boundParams);
  return jacToGlobal;
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
void applyDerivativeCorrectionBound(const FreeVector& freeParams,
                                    const FreeVector& derivative,
                                    const Surface& surface,
                                    const RotationMatrix3D& referenceFrameT,
                                    const GeometryContext& geoContext,
                                    BoundToFreeMatrix& jacToGlobal) {
  // Calculate the form factors for the derivatives
  const BoundRowVector sfactors = surface.derivativeFactors(
      geoContext, freeParams.segment<3>(eFreePos0),
      freeParams.segment<3>(eFreeDir0), referenceFrameT, jacToGlobal);
  jacToGlobal -= derivative * sfactors;
}

}  // namespace

void updateJacobiansToBound(const FreeVector& freeParams,
                            const Surface& surface,
                            const GeometryContext& geoCtx,
                            BoundToFreeMatrix& jacToGlobal,
                            FreeMatrix& jacTransport, FreeVector& derivative,
                            BoundMatrix& jacobian) {
  // Jacobian from current free parameters onto the surface
  FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
  const auto& referenceFrameTranspose = surface.initJacobianToLocal(
      geoCtx, jacToLocal, freeParams.segment<3>(eFreePos0),
      freeParams.segment<3>(eFreeDir0));
  // Add corrections due to propagation derivatives
  // TODO 20200319 msmk: why does this change the initial projection and not the
  //                     final one?
  applyDerivativeCorrectionBound(freeParams, derivative, surface,
                                 referenceFrameTranspose, geoCtx, jacToGlobal);
  // combine Jacobian for bound parameters from:
  // 1. projection from the initial bound frame to free parameters
  // 2. free transport Jacobian along trajectory as computed by stepper
  // 3. Jacobian onto the bound frame on the surface at the current position
  // matrix multiplication uses right-to-left order, i.e. 3 * 2 * 1
  jacobian = jacToLocal * jacTransport * jacToGlobal;

  // reset Jacobians for a propagation starting on the surface
  jacToGlobal = boundToFreeJacobian(freeParams, surface, geoCtx);
  jacTransport = FreeMatrix::Identity();
  derivative = FreeVector::Zero();
}

namespace {

/// Projection Jacobian from curvilinear parameters to free parameters.
///
/// In principle, this should be equivalent to the projection Jacobian for
/// general bound parameters with an appropriate Curvilinear surface. Use an
/// explicit implementation for faster computation.
BoundToFreeMatrix curvilinearToFreeJacobian(const FreeVector& freeParams) {
  // normalized direction vector components expressed as phi/theta
  const double x = freeParams[eFreeDir0];  // == cos(phi) * sin(theta)
  const double y = freeParams[eFreeDir1];  // == sin(phi) * sin(theta)
  const double z = freeParams[eFreeDir2];  // == cos(theta)
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
FreeToBoundMatrix freeToCurvilinearJacobian(const FreeVector& freeParams) {
  // normalized direction vector components expressed as phi/theta
  const double x = freeParams[eFreeDir0];  // == cos(phi) * sin(theta)
  const double y = freeParams[eFreeDir1];  // == sin(phi) * sin(theta)
  const double z = freeParams[eFreeDir2];  // == cos(theta)
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

/// Apply Jacobian modifications from projection onto the curvilinear frame.
///
/// @param [in] unitDirection Particle direction
/// @param [in] derivative Stepper propagation derivative
/// @param [in,out] jacToGlobal curvilinear-to-free Jacobian to modify
void applyDerivativeCorrectionCurvilinear(const FreeVector& freeParams,
                                          const FreeVector& derivative,
                                          BoundToFreeMatrix& jacToGlobal) {
  // Transport the covariance
  const BoundRowVector sfactors =
      freeParams.segment<3>(eFreeDir0).transpose() *
      jacToGlobal.template topLeftCorner<3, eBoundParametersSize>();
  jacToGlobal -= derivative * sfactors;
}

}  // namespace

void updateJacobiansToCurvilinear(const FreeVector& freeParams,
                                  BoundToFreeMatrix& jacToGlobal,
                                  FreeMatrix& jacTransport,
                                  FreeVector& derivative,
                                  BoundMatrix& jacobian) {
  // Jacobian from current free parameters onto the curvilinear frame
  const FreeToBoundMatrix jacToLocal = freeToCurvilinearJacobian(freeParams);
  // Add corrections due to propagation derivatives
  // TODO 20200319 msmk: why does this change the initial projection and not the
  //                     final one?
  applyDerivativeCorrectionCurvilinear(freeParams, derivative, jacToGlobal);
  // combine Jacobian for bound parameters from:
  // 1. projection from the initial bound frame to free parameters
  // 2. free transport Jacobian along trajectory as computed by stepper
  // 3. Jacobian onto the curvilinear frame at the current position
  // matrix multiplication uses right-to-left order, i.e. 3 * 2 * 1
  jacobian = jacToLocal * jacTransport * jacToGlobal;

  // reset Jacobians for a propagation starting on the curvilinear frame
  jacToGlobal = curvilinearToFreeJacobian(freeParams);
  jacTransport = FreeMatrix::Identity();
  derivative = FreeVector::Zero();
}

}  // namespace detail
}  // namespace Acts
