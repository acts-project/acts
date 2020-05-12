// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Alignment/detail/AlignmentEngine.hpp"

namespace Acts {
namespace detail {

AlignmentToBoundMatrix alignmentToLocalDerivative(
    const GeometryContext& gctx, const BoundParameters& boundParams,
    const FreeVector& derivatives) const {
  // Initialize the derivative of bound parameters w.r.t. alignment parameters
  AlignmentToBoundMatrix alignToLocal = AlignmentToBoundMatrix::Zero();

  // The reference surface
  auto surface = &boundParams.referenceSurface();
  // Transform of the reference surface
  const auto& transform = surface.transform(gctx);
  // Get the translation and rotation
  const auto& translation = transform->translation();
  const auto& rotation = transform->matrix();

  // The position and momentum direction
  const auto pos = boundParams.position();
  const auto dir = boundParams.momentum().normalized();
  // Initialize the jacobian from free parameters to bound parameters
  FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
  // Initalize the jacobian to local, returns the transposed ref frame (i.e.
  // measurement frame)
  auto rframeT = surface->initJacobianToLocal(geoContext, jacToLocal, pos, dir);
  auto rframe = rframeT.transpose();
  // The transformation from alignment frame to measurement frame (need for
  // derivative of bound parameters w.r.t. rotation)
  auto alignToLocTransform = rframe * rotation.inverse();
  // The axis of measurement frame
  auto localXAxis = rframe.col(0);
  auto localYAxis = rframe.col(1);
  auto localZAxis = rframe.col(2);

  // The cosine of angle between momentum direction and z axis of measurement
  // frame
  const double cosThetaDir = localZAxis * dir;
  // Derivative of propagation path w.r.t. geometry object position
  const auto posToPathDerVec = localZAxis / cosThetaDir;
  // Derivative of eLOC_0
  const auto jacToLoc0 = jacToLocal.block<1, 8>(0, 0);
  alignToLocal(eLOC_0, ePos_X) =
      -1.0 * localXAxis(0) + jacToLoc0 * derivatives * posToPathDerVec(0);
  alignToLocal(eLOC_0, ePos_Y) =
      -1.0 * localXAxis(1) + jacToLoc0 * derivatives * posToPathDerVec(1);
  alignToLocal(eLOC_0, ePos_Z) =
      -1.0 * localXAxis(2) + jacToLoc0 * derivatives * posToPathDerVec(2);
  // Derivative of eLOC_1
  const auto jacToLoc1 = jacToLocal.block<1, 8>(1, 0);
  alignToLocal(eLOC_1, ePos_X) =
      -1.0 * localYAxis(0) + jacToLoc1 * derivatives * posToPathDerVec(0);
  alignToLocal(eLOC_1, ePos_Y) =
      -1.0 * localYAxis(1) + jacToLoc1 * derivatives * posToPathDerVec(1);
  alignToLocal(eLOC_1, ePos_Z) =
      -1.0 * localYAxis(2) + jacToLoc1 * derivatives * posToPathDerVec(2);
  // Derivative of ePHI
  const auto jacToPhi = jacToLocal.block<1, 8>(2, 0);
  alignToLocal(ePHI, ePos_X) = jacToPhi * derivatives * posToPathDerVec(0);
  alignToLocal(ePHI, ePos_Y) = jacToPhi * derivatives * posToPathDerVec(1);
  alignToLocal(ePHI, ePos_Z) = jacToPhi * derivatives * posToPathDerVec(2);
  // Derivative of eTHETA
  const auto jacToTheta = jacToLocal.block<1, 8>(3, 0);
  alignToLocal(eTHETA, ePos_X) = jacToTheta * derivatives * posToPathDerVec(0);
  alignToLocal(eTHETA, ePos_Y) = jacToTheta * derivatives * posToPathDerVec(1);
  alignToLocal(eTHETA, ePos_Z) = jacToTheta * derivatives * posToPathDerVec(2);
  // Derivative of eQOP
  const auto jacToQoP = jacToLocal.block<1, 8>(4, 0);
  alignToLocal(eQOP, ePos_X) = jacToQoP * derivatives * posToPathDerVec(0);
  alignToLocal(eQOP, ePos_Y) = jacToQoP * derivatives * posToPathDerVec(1);
  alignToLocal(eQOP, ePos_Z) = jacToQoP * derivatives * posToPathDerVec(2);
  // Derivative of eT
  const auto jacToT = jacToLocal.block<1, 8>(5, 0);
  alignToLocal(eT, ePos_X) = jacToT * derivatives * posToPathDerVec(0);
  alignToLocal(eT, ePos_Y) = jacToT * derivatives * posToPathDerVec(1);
  alignToLocal(eT, ePos_Z) = jacToT * derivatives * posToPathDerVec(2);

  //@Todo: calculate derivative w.r.t. geometry object rotation

  return alignToLocal;
}

}  // namespace detail
}  // namespace Acts
