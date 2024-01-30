// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/JacobianEngine.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <cmath>

namespace Acts {

FreeToBoundMatrix detail::freeToCurvilinearJacobian(const Vector3& direction) {
  // Optimized trigonometry on the propagation direction
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = std::hypot(x, y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;
  // prepare the jacobian to curvilinear
  FreeToBoundMatrix jacToCurv = FreeToBoundMatrix::Zero();
  if (std::abs(cosTheta) < s_curvilinearProjTolerance) {
    // We normally operate in curvilinear coordinates defined as follows
    jacToCurv(0, 0) = -sinPhi;
    jacToCurv(0, 1) = cosPhi;
    jacToCurv(1, 0) = -cosPhi * cosTheta;
    jacToCurv(1, 1) = -sinPhi * cosTheta;
    jacToCurv(1, 2) = sinTheta;
  } else {
    // Under grazing incidence to z, the above coordinate system definition
    // becomes numerically unstable, and we need to switch to another one
    const double c = std::hypot(y, z);
    const double invC = 1. / c;
    jacToCurv(0, 1) = -z * invC;
    jacToCurv(0, 2) = y * invC;
    jacToCurv(1, 0) = c;
    jacToCurv(1, 1) = -x * y * invC;
    jacToCurv(1, 2) = -x * z * invC;
  }
  // Time parameter
  jacToCurv(5, 3) = 1.;
  // Directional and momentum parameters for curvilinear
  jacToCurv(2, 4) = -sinPhi * invSinTheta;
  jacToCurv(2, 5) = cosPhi * invSinTheta;
  jacToCurv(3, 4) = cosPhi * cosTheta;
  jacToCurv(3, 5) = sinPhi * cosTheta;
  jacToCurv(3, 6) = -sinTheta;
  jacToCurv(4, 7) = 1.;

  return jacToCurv;
}

BoundToFreeMatrix detail::curvilinearToFreeJacobian(const Vector3& direction) {
  auto [cosPhi, sinPhi, cosTheta, sinTheta, invSinTheta] =
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

BoundMatrix detail::boundToBoundTransportJacobian(
    const GeometryContext& geoContext, const FreeVector& freeParameters,
    const BoundToFreeMatrix& boundToFreeJacobian,
    const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives, const Surface& surface) {
  BoundMatrix result;
  detail::boundToBoundTransportJacobian(
      geoContext, surface, freeParameters, boundToFreeJacobian,
      freeTransportJacobian, freeToPathDerivatives, result);
  return result;
}

void detail::boundToCurvilinearTransportJacobian(
    const Vector3& direction, const BoundToFreeMatrix& boundToFreeJacobian,
    const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives,
    BoundMatrix& fullTransportJacobian) {
  // Calculate the jacobian from global to local at the curvilinear surface
  FreeToBoundMatrix freeToBoundJacobian = freeToCurvilinearJacobian(direction);

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
  return freeToCurvilinearJacobian(direction) *
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
  boundToFreeJacobian = BoundToFreeMatrix::Zero();

  // Optimized trigonometry on the propagation direction
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = std::hypot(x, y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;

  boundToFreeJacobian(eFreePos0, eBoundLoc0) = -sinPhi;
  boundToFreeJacobian(eFreePos0, eBoundLoc1) = -cosPhi * cosTheta;
  boundToFreeJacobian(eFreePos1, eBoundLoc0) = cosPhi;
  boundToFreeJacobian(eFreePos1, eBoundLoc1) = -sinPhi * cosTheta;
  boundToFreeJacobian(eFreePos2, eBoundLoc1) = sinTheta;
  boundToFreeJacobian(eFreeTime, eBoundTime) = 1;
  boundToFreeJacobian(eFreeDir0, eBoundPhi) = -sinTheta * sinPhi;
  boundToFreeJacobian(eFreeDir0, eBoundTheta) = cosTheta * cosPhi;
  boundToFreeJacobian(eFreeDir1, eBoundPhi) = sinTheta * cosPhi;
  boundToFreeJacobian(eFreeDir1, eBoundTheta) = cosTheta * sinPhi;
  boundToFreeJacobian(eFreeDir2, eBoundTheta) = -sinTheta;
  boundToFreeJacobian(eFreeQOverP, eBoundQOverP) = 1;
}

}  // namespace Acts
