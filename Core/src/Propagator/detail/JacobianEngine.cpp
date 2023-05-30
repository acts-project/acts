// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

namespace Acts {

namespace detail {

FreeToBoundMatrix freeToCurvilinearJacobian(const Vector3& direction) {
  auto [cosPhi, sinPhi, cosTheta, sinTheta, invSinTheta] =
      VectorHelpers::evaluateTrigonomics(direction);
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

BoundToFreeMatrix curvilinearToFreeJacobian(const Vector3& direction) {
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

ActsMatrix<8, 7> anglesToDirectionJacobian(const Vector3& direction) {
  ActsMatrix<8, 7> jacobian = ActsMatrix<8, 7>::Zero();

  auto [cosPhi, sinPhi, cosTheta, sinTheta, invSinTheta] =
      VectorHelpers::evaluateTrigonomics(direction);

  jacobian(0, 0) = 1.;
  jacobian(1, 1) = 1.;
  jacobian(2, 2) = 1.;
  jacobian(3, 3) = 1.;
  jacobian(7, 6) = 1.;

  jacobian(4, 4) = -sinTheta * sinPhi;
  jacobian(4, 5) = cosTheta * cosPhi;
  jacobian(5, 4) = sinTheta * cosPhi;
  jacobian(5, 5) = cosTheta * sinPhi;
  jacobian(6, 5) = -sinTheta;
  return jacobian;
}

ActsMatrix<7, 8> directionToAnglesJacobian(const Vector3& direction) {
  ActsMatrix<7, 8> jacobian = ActsMatrix<7, 8>::Zero();

  auto [cosPhi, sinPhi, cosTheta, sinTheta, invSinTheta] =
      VectorHelpers::evaluateTrigonomics(direction);

  jacobian(0, 0) = 1.;
  jacobian(1, 1) = 1.;
  jacobian(2, 2) = 1.;
  jacobian(3, 3) = 1.;
  jacobian(6, 7) = 1.;

  jacobian(4, 4) = -sinPhi * invSinTheta;
  jacobian(4, 5) = cosPhi * invSinTheta;
  jacobian(5, 4) = cosPhi * cosTheta;
  jacobian(5, 5) = sinPhi * cosTheta;
  jacobian(5, 6) = -sinTheta;

  return jacobian;
}

BoundMatrix boundToBoundTransportJacobian(
    const GeometryContext& geoContext, const FreeVector& freeParameters,
    const BoundToFreeMatrix& boundToFreeJacobian,
    const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives, const Surface& surface) {
  // Calculate the derivative of path length at the final surface or the
  // point-of-closest approach w.r.t. free parameters
  const FreeToPathMatrix freeToPath =
      surface.freeToPathDerivative(geoContext, freeParameters);
  // Calculate the jacobian from free to bound at the final surface
  FreeToBoundMatrix freeToBoundJacobian =
      surface.freeToBoundJacobian(geoContext, freeParameters);
  // Calculate the full jacobian from the local/bound parameters at the start
  // surface to local/bound parameters at the final surface
  return freeToBoundJacobian *
         (FreeMatrix::Identity() + freeToPathDerivatives * freeToPath) *
         freeTransportJacobian * boundToFreeJacobian;
}

BoundMatrix boundToCurvilinearTransportJacobian(
    const Vector3& direction, const BoundToFreeMatrix& boundToFreeJacobian,
    const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives) {
  // Calculate the derivative of path length at the curvilinear surface
  // w.r.t. free parameters
  FreeToPathMatrix freeToPath = FreeToPathMatrix::Zero();
  freeToPath.segment<3>(eFreePos0) = -1.0 * direction;
  // Calculate the jacobian from global to local at the curvilinear surface
  FreeToBoundMatrix freeToBoundJacobian = freeToCurvilinearJacobian(direction);
  // Calculate the full jacobian from the local parameters at the start surface
  // to curvilinear parameters
  return freeToBoundJacobian *
         (FreeMatrix::Identity() + freeToPathDerivatives * freeToPath) *
         freeTransportJacobian * boundToFreeJacobian;
}

BoundToFreeMatrix boundToFreeTransportJacobian(
    const BoundToFreeMatrix& boundToFreeJacobian,
    const FreeMatrix& freeTransportJacobian) {
  // Calculate the full jacobian, in this case simple a product of
  // jacobian(transport in free) * jacobian(bound to free)
  return (freeTransportJacobian * boundToFreeJacobian);
}

FreeToBoundMatrix freeToBoundTransportJacobian(
    const GeometryContext& geoContext, const FreeVector& freeParameters,
    const ActsMatrix<7, 8>& directionToAnglesJacobian,
    const ActsMatrix<8, 7>& anglesToDirectionJacobian,
    const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives, const Surface& surface) {
  // Calculate the jacobian from free to bound at the final surface
  FreeToBoundMatrix freeToBoundJacobian =
      surface.freeToBoundJacobian(geoContext, freeParameters);
  // Calculate the form factors for the derivatives
  auto transport = freeTransportJacobian * anglesToDirectionJacobian;
  FreeToPathMatrix sVec =
      surface.freeToPathDerivative(geoContext, freeParameters);
  // Return the jacobian to local
  return freeToBoundJacobian *
         (transport + freeToPathDerivatives * sVec * transport) *
         directionToAnglesJacobian;
}

FreeToBoundMatrix freeToCurvilinearTransportJacobian(
    const Vector3& direction, const ActsMatrix<7, 8>& directionToAnglesJacobian,
    const ActsMatrix<8, 7>& anglesToDirectionJacobian,
    const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives) {
  const ActsMatrix<8, 7> transport =
      freeTransportJacobian * anglesToDirectionJacobian;
  auto sfactors =
      direction.transpose() * transport.template topLeftCorner<3, 7>();

  // Since the jacobian to local needs to calculated for the bound parameters
  // here, it is convenient to do the same here
  return freeToCurvilinearJacobian(direction) *
         (transport - freeToPathDerivatives * sfactors) *
         directionToAnglesJacobian;
}

FreeMatrix freeToFreeTransportJacobian(
    const ActsMatrix<7, 8>& directionToAnglesJacobian,
    const ActsMatrix<8, 7>& anglesToDirectionJacobian,
    const FreeMatrix& freeTransportJacobian) {
  return freeTransportJacobian * anglesToDirectionJacobian *
         directionToAnglesJacobian;
}

}  // namespace detail
}  // namespace Acts
