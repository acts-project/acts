// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/JacobianEngine.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AlgebraHelpers.hpp"

namespace Acts {

FreeToBoundMatrix detail::freeToCurvilinearJacobian(const Vector3& direction) {
  auto [cosPhi, sinPhi, cosTheta, sinTheta] =
      VectorHelpers::evaluateTrigonomics(direction);
  ActsScalar invSinTheta = 1. / sinTheta;
  // Prepare the jacobian to curvilinear
  FreeToBoundMatrix freeToCurvJacobian = FreeToBoundMatrix::Zero();
  if (std::abs(cosTheta) < s_curvilinearProjTolerance) {
    // We normally operate in curvilinear coordinates defined as follows
    freeToCurvJacobian(eBoundLoc0, eFreePos0) = -sinPhi;
    freeToCurvJacobian(eBoundLoc0, eFreePos1) = cosPhi;
    freeToCurvJacobian(eBoundLoc1, eFreePos0) = -cosPhi * cosTheta;
    freeToCurvJacobian(eBoundLoc1, eFreePos1) = -sinPhi * cosTheta;
    freeToCurvJacobian(eBoundLoc1, eFreePos2) = sinTheta;
  } else {
    // Under grazing incidence to z, the above coordinate system definition
    // becomes numerically unstable, and we need to switch to another one
    const ActsScalar x = direction(0);  // == cos(phi) * sin(theta)
    const ActsScalar y = direction(1);  // == sin(phi) * sin(theta)
    const ActsScalar z = direction(2);  // == cos(theta)
    const ActsScalar c = std::hypot(y, z);
    const ActsScalar invC = 1. / c;
    freeToCurvJacobian(eBoundLoc0, eFreePos1) = -z * invC;
    freeToCurvJacobian(eBoundLoc0, eFreePos2) = y * invC;
    freeToCurvJacobian(eBoundLoc1, eFreePos0) = c;
    freeToCurvJacobian(eBoundLoc1, eFreePos1) = -x * y * invC;
    freeToCurvJacobian(eBoundLoc1, eFreePos2) = -x * z * invC;
  }
  // Time parameter
  freeToCurvJacobian(eBoundTime, eFreeTime) = 1.;
  // Directional and momentum parameters for curvilinear
  freeToCurvJacobian(eBoundPhi, eFreeDir0) = -sinPhi * invSinTheta;
  freeToCurvJacobian(eBoundPhi, eFreeDir1) = cosPhi * invSinTheta;
  freeToCurvJacobian(eBoundTheta, eFreeDir0) = cosPhi * cosTheta;
  freeToCurvJacobian(eBoundTheta, eFreeDir1) = sinPhi * cosTheta;
  freeToCurvJacobian(eBoundTheta, eFreeDir2) = -sinTheta;
  freeToCurvJacobian(eBoundQOverP, eFreeQOverP) = 1.;

  return freeToCurvJacobian;
}

BoundToFreeMatrix detail::curvilinearToFreeJacobian(const Vector3& direction) {
  auto [cosPhi, sinPhi, cosTheta, sinTheta] =
      VectorHelpers::evaluateTrigonomics(direction);

  // Prepare the jacobian to free
  BoundToFreeMatrix curvToFreeJacobian = BoundToFreeMatrix::Zero();

  curvToFreeJacobian(eFreePos0, eBoundLoc0) = -sinPhi;
  curvToFreeJacobian(eFreePos0, eBoundLoc1) = -cosPhi * cosTheta;
  curvToFreeJacobian(eFreePos1, eBoundLoc0) = cosPhi;
  curvToFreeJacobian(eFreePos1, eBoundLoc1) = -sinPhi * cosTheta;
  curvToFreeJacobian(eFreePos2, eBoundLoc1) = sinTheta;
  // Time parameter: stays as is
  curvToFreeJacobian(eFreeTime, eBoundTime) = 1;
  curvToFreeJacobian(eFreeDir0, eBoundPhi) = -sinTheta * sinPhi;
  curvToFreeJacobian(eFreeDir0, eBoundTheta) = cosTheta * cosPhi;
  curvToFreeJacobian(eFreeDir1, eBoundPhi) = sinTheta * cosPhi;
  curvToFreeJacobian(eFreeDir1, eBoundTheta) = cosTheta * sinPhi;
  curvToFreeJacobian(eFreeDir2, eBoundTheta) = -sinTheta;
  // Q/P parameter: stays as is
  curvToFreeJacobian(eFreeQOverP, eBoundQOverP) = 1;

  return curvToFreeJacobian;
}

void detail::boundToBoundTransportJacobian(
    const GeometryContext& geoContext, const Surface& surface,
    const FreeVector& freeParameters,
    const BoundToFreeMatrix& boundToFreeJacobian,
    const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives,
    BoundMatrix& fullTransportJacobian) {
  // Calculate the derivative of path length at the final surface or the
  // point-of-closest approach w.r.t. free parameters
  const FreeToPathMatrix freeToPath =
      surface.freeToPathDerivative(geoContext, freeParameters);
  // Calculate the jacobian from free to bound at the final surface
  FreeToBoundMatrix freeToBoundJacobian =
      surface.freeToBoundJacobian(geoContext, freeParameters);
  // https://acts.readthedocs.io/en/latest/white_papers/correction-for-transport-jacobian.html
  // Calculate the full jacobian from the local/bound parameters at the start
  // surface to local/bound parameters at the final surface
  // @note jac(locA->locB) = jac(gloB->locB)*(1+
  // pathCorrectionFactor(gloB))*jacTransport(gloA->gloB) *jac(locA->gloA)
  fullTransportJacobian =
      freeToBoundJacobian *
      (FreeMatrix::Identity() + freeToPathDerivatives * freeToPath) *
      freeTransportJacobian * boundToFreeJacobian;
}

void detail::boundToCurvilinearTransportJacobian(
    const Vector3& direction, const BoundToFreeMatrix& boundToFreeJacobian,
    const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives,
    BoundMatrix& fullTransportJacobian) {
  // Calculate the jacobian from global to local at the curvilinear surface
  FreeToBoundMatrix freeToBoundJacobian =
      CurvilinearSurface(direction).freeToBoundJacobian();

  // Update the jacobian to include the derivative of the path length at the
  // curvilinear surface w.r.t. the free parameters
  freeToBoundJacobian.topLeftCorner<6, 3>() +=
      (freeToBoundJacobian * freeToPathDerivatives) *
      (-1.0 * direction).transpose();

  // Calculate the full jocobian from the local parameters at the start surface
  // to curvilinear parameters
  // @note jac(locA->locB) = jac(gloB->locB)*(1+
  // pathCorrectionFactor(gloB))*jacTransport(gloA->gloB) *jac(locA->gloA)
  fullTransportJacobian =
      blockedMult(freeToBoundJacobian,
                  blockedMult(freeTransportJacobian, boundToFreeJacobian));
}

BoundToFreeMatrix detail::boundToFreeTransportJacobian(
    const BoundToFreeMatrix& boundToFreeJacobian,
    const FreeMatrix& freeTransportJacobian) {
  // Calculate the full jacobian, in this case simple a product of
  // jacobian(transport in free) * jacobian(bound to free)
  return freeTransportJacobian * boundToFreeJacobian;
}

void detail::freeToBoundTransportJacobian(
    const GeometryContext& geoContext, const Surface& surface,
    const FreeVector& freeParameters, const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives,
    FreeToBoundMatrix& fullTransportJacobian) {
  // Calculate the jacobian from free to bound at the final surface
  FreeToBoundMatrix freeToBoundJacobian =
      surface.freeToBoundJacobian(geoContext, freeParameters);
  FreeToPathMatrix sVec =
      surface.freeToPathDerivative(geoContext, freeParameters);
  // Return the jacobian to local
  fullTransportJacobian = freeToBoundJacobian * (freeTransportJacobian +
                                                 freeToPathDerivatives * sVec *
                                                     freeTransportJacobian);
}

FreeToBoundMatrix detail::freeToCurvilinearTransportJacobian(
    const Vector3& direction, const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives) {
  auto sfactors = direction.transpose() *
                  freeTransportJacobian.template topLeftCorner<3, 8>();

  // Since the jacobian to local needs to calculated for the bound parameters
  // here, it is convenient to do the same here
  return CurvilinearSurface(direction).freeToBoundJacobian() *
         (freeTransportJacobian - freeToPathDerivatives * sfactors);
}

Result<void> detail::reinitializeJacobians(
    const GeometryContext& geoContext, const Surface& surface,
    FreeMatrix& freeTransportJacobian, FreeVector& freeToPathDerivatives,
    BoundToFreeMatrix& boundToFreeJacobian, const FreeVector& freeParameters) {
  // Reset the jacobians
  freeTransportJacobian = FreeMatrix::Identity();
  freeToPathDerivatives = FreeVector::Zero();

  // Get the local position
  const Vector3 position = freeParameters.segment<3>(eFreePos0);
  const Vector3 direction = freeParameters.segment<3>(eFreeDir0);
  auto lpResult = surface.globalToLocal(geoContext, position, direction);
  if (!lpResult.ok()) {
    return lpResult.error();
  }
  // Transform from free to bound parameters
  Result<BoundVector> boundParameters = detail::transformFreeToBoundParameters(
      freeParameters, surface, geoContext);
  if (!boundParameters.ok()) {
    return boundParameters.error();
  }
  // Reset the jacobian from local to global
  boundToFreeJacobian =
      surface.boundToFreeJacobian(geoContext, *boundParameters);
  return Result<void>::success();
}

void detail::reinitializeJacobians(FreeMatrix& freeTransportJacobian,
                                   FreeVector& freeToPathDerivatives,
                                   BoundToFreeMatrix& boundToFreeJacobian,
                                   const Vector3& direction) {
  // Reset the jacobians
  freeTransportJacobian = FreeMatrix::Identity();
  freeToPathDerivatives = FreeVector::Zero();
  boundToFreeJacobian = CurvilinearSurface(direction).boundToFreeJacobian();
}

}  // namespace Acts
