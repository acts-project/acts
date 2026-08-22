// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TransformationHelpers.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

Acts::FreeVector Acts::transformBoundToFreeParameters(
    const Acts::Surface& surface, const GeometryContext& geoCtx,
    const Acts::BoundVector& boundParams) {
  // convert angles to global unit direction vector
  Vector3 direction = makeDirectionFromPhiTheta(boundParams[eBoundPhi],
                                                boundParams[eBoundTheta]);

  // convert local position to global position vector
  Vector2 local(boundParams[eBoundLoc0], boundParams[eBoundLoc1]);
  Vector3 position = surface.localToGlobal(geoCtx, local, direction);

  // construct full free-vector. time and q/p stay as-is.
  FreeVector freeParams = FreeVector::Zero();
  freeParams[eFreePos0] = position[ePos0];
  freeParams[eFreePos1] = position[ePos1];
  freeParams[eFreePos2] = position[ePos2];
  freeParams[eFreeTime] = boundParams[eBoundTime];
  freeParams[eFreeDir0] = direction[eMom0];
  freeParams[eFreeDir1] = direction[eMom1];
  freeParams[eFreeDir2] = direction[eMom2];
  freeParams[eFreeQOverP] = boundParams[eBoundQOverP];
  return freeParams;
}

Acts::Result<Acts::BoundVector> Acts::transformFreeToBoundParameters(
    const FreeVector& freeParams, const Surface& surface,
    const GeometryContext& geoCtx, double tolerance) {
  // initialize the bound vector
  BoundVector bp = BoundVector::Zero();
  // convert global to local position on the surface
  auto position = freeParams.segment<3>(eFreePos0);
  auto direction = freeParams.segment<3>(eFreeDir0);
  auto result = surface.globalToLocal(geoCtx, position, direction, tolerance);
  if (!result.ok()) {
    return Result<Acts::BoundVector>::failure(result.error());
  }

  auto localPosition = result.value();
  bp[eBoundLoc0] = localPosition[ePos0];
  bp[eBoundLoc1] = localPosition[ePos1];

  bp[eBoundTime] = freeParams[eFreeTime];
  bp[eBoundPhi] = VectorHelpers::phi(direction);
  bp[eBoundTheta] = VectorHelpers::theta(direction);
  bp[eBoundQOverP] = freeParams[eFreeQOverP];
  return Result<Acts::BoundVector>::success(bp);
}

Acts::Result<Acts::BoundVector> Acts::transformFreeToBoundParameters(
    const Acts::Vector3& position, double time, const Acts::Vector3& direction,
    double qOverP, const Acts::Surface& surface,
    const Acts::GeometryContext& geoCtx, double tolerance) {
  // initialize the bound vector
  BoundVector bp = BoundVector::Zero();
  // convert global to local position on the surface
  auto result = surface.globalToLocal(geoCtx, position, direction, tolerance);
  if (!result.ok()) {
    return Result<Acts::BoundVector>::failure(result.error());
  }

  auto localPosition = result.value();
  bp[eBoundLoc0] = localPosition[ePos0];
  bp[eBoundLoc1] = localPosition[ePos1];

  bp[eBoundTime] = time;
  bp[eBoundPhi] = VectorHelpers::phi(direction);
  bp[eBoundTheta] = VectorHelpers::theta(direction);
  bp[eBoundQOverP] = qOverP;
  return Result<Acts::BoundVector>::success(bp);
}

Acts::BoundVector Acts::transformFreeToCurvilinearParameters(double time,
                                                             double phi,
                                                             double theta,
                                                             double qOverP) {
  BoundVector bp = BoundVector::Zero();
  // local coordinates are zero by construction
  bp[eBoundTime] = time;
  bp[eBoundPhi] = phi;
  bp[eBoundTheta] = theta;
  bp[eBoundQOverP] = qOverP;
  return bp;
}

Acts::BoundVector Acts::transformFreeToCurvilinearParameters(
    double time, const Vector3& direction, double qOverP) {
  BoundVector bp = BoundVector::Zero();
  // local coordinates are zero by construction
  bp[eBoundTime] = time;
  bp[eBoundPhi] = VectorHelpers::phi(direction);
  bp[eBoundTheta] = VectorHelpers::theta(direction);
  bp[eBoundQOverP] = qOverP;
  return bp;
}

std::pair<Acts::Vector4, Acts::SquareMatrix<7>>
Acts::transformBoundToCartesianFourPositionMomentum(
    const GeometryContext& gctx, const BoundTrackParameters& params,
    Vector3& momentum) {
  const Vector4 pos4 = params.fourPosition(gctx);
  const Vector3 dir = params.direction();
  momentum = params.momentum();

  // d(FreeVector)/d(BoundVector)
  // rows ordered [pos0,pos1,pos2, time, dir0,dir1,dir2, qOverP].
  const BoundToFreeMatrix jacToFree =
      params.referenceSurface().boundToFreeJacobian(
          gctx, pos4.segment<3>(ePos0), dir);

  const double p = params.absoluteMomentum();
  const double qOverP = params.qOverP();
  // since p = charge/qOverP
  const double dpdqop = -p / qOverP;

  // d(px,py,pz)/d(freeDir0,freeDir1,freeDir2,qOverP)
  Matrix<3, 8> momentumToFree = Matrix<3, 8>::Zero();
  momentumToFree.block<3, 3>(0, eFreeDir0) = p * SquareMatrix3::Identity();
  momentumToFree.block<3, 1>(0, eFreeQOverP) = dpdqop * dir;

  // d(x,y,z,t,px,py,pz)/d(FreeVector)
  Matrix<7, 8> freeToCartesian = Matrix<7, 8>::Zero();
  freeToCartesian.block<4, 4>(0, 0) = SquareMatrix<4>::Identity();
  freeToCartesian.block<3, 8>(4, 0) = momentumToFree;

  // d(x,y,z,t,px,py,pz)/d(BoundVector)
  const Matrix<7, 6> jacToCartesian = freeToCartesian * jacToFree;

  const SquareMatrix<7> cov7 =
      jacToCartesian * params.covariance().value() * jacToCartesian.transpose();

  return {pos4, cov7};
}
