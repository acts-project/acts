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
    const FreeVector& derivatives, const Vector3D& rframeOrigin,
    const Cartesian3DToLocal2DMatrix& cartesianToLocal) const {
  // The position and momentum direction
  const auto pos = boundParams.position();
  const auto dir = boundParams.momentum().normalized();
  // The reference surface
  auto surface = &boundParams.referenceSurface();

  // Initialize the jacobian from free parameters to bound parameters
  FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
  // Set the jacobian to local, returns the transposed ref frame (i.e.
  // measurement frame)
  auto rframeT = surface->initJacobianToLocal(geoContext, jacToLocal, pos, dir);
  auto rframe = rframeT.transpose();
  // The cartesian axis of local measurement frame
  auto localXAxis = rframe.col(0);
  auto localYAxis = rframe.col(1);
  auto localZAxis = rframe.col(2);

  // Get Euler angles for rotation representated by rotZ * rotY * rotX, i.e.
  // first rotation around x axis, then y axis, last z axis
  //@Todo: add parameter for rotation order
  const Vector3D rotAngles = rframe.eulerAngles(0, 1, 2);
  double alpha = rotAngles(0);  // rot_X
  double beta = rotAngles(1);   // rot_Y
  double gamma = rotAngles(2);  // rot_Z
  double sinAlpha = std::sin(alpha);
  double cosAlpha = std::cos(alpha);
  double sinBeta = std::sin(beta);
  double cosBeta = std::cos(beta);
  double sinGam = std::sin(gamma);
  double cosGam = std::cos(gamma);
  //@Todo: add function for the following calculation
  RotationMatrix3D rotToLocalXAxisDer = RotationMatrix3D::Zero();
  rotToLocalXAxisDer.col(0) =
      (0, 0, 0);  // Derivative of local x axis w.r.t. rot_X
  rotToLocalXAxisDer.col(1) =
      (-sinBeta * cosGam, -sinBeta * sinGam,
       -cosBeta);  // Derivative of local x axis w.r.t. rot_Y
  rotToLocalXAxisDer.col(2) = (-cosBeta * sinGam, cosBeta * cosGam,
                               0);  // Derivative of local x axis w.r.t. rot_Z
  RotationMatrix3D rotToLocalYAxisDer = RotationMatrix3D::Zero();
  rotToLocalYAxisDer.col(0) =
      (cosAlpha * sinBeta * cosGam + sinAlpha * sinGam,
       cosAlpha * sinBeta * sinGam - sinAlpha * cosGam,
       cosAlpha * cosBeta);  // Derivative of local y axis w.r.t. rot_X
  rotToLocalYAxisDer.col(1) =
      (sinAlpha * cosBeta * cosGam, sinAlpha * cosBeta * sinGam,
       -sinAlpha * sinBeta);  // Derivative of local y axis w.r.t. rot_Y
  rotToLocalYAxisDer.col(2) = (-sinAlpha * sinBeta * sinGam - cosAlpha * cosGam,
                               sinAlpha * sinBeta * cosGam - cosAlpha * sinGam,
                               0);  // Derivative of local y axis w.r.t. rot_z
  RotationMatrix3D rotToLocalZAxisDer = RotationMatrix3D::Zero();
  rotToLocalZAxisDer.col(0) =
      (-sinAlpha * sinBeta * cosGam + cosAlpha * sinGam,
       -sinAlpha * sinBeta * sinGam - cosAlpha * cosGam,
       -sinAlpha * cosBeta);  // Derivative of local z axis w.r.t. rot_X
  rotToLocalZAxisDer.col(1) =
      (cosAlpha * cosBeta * cosGam, cosAlpha * cosBeta * sinGam,
       -cosAlpha * sinBeta);  // Derivative of local z axis w.r.t. rot_Y
  rotToLocalZAxisDer.col(2) = (-cosAlpha * sinBeta * sinGam + sinAlpha * cosGam,
                               cosAlpha * sinBeta * cosGam + sinAlpha * sinGam,
                               0);  // Derivative of local z axis w.r.t. rot_Z

  // The vector of track position in reference frame
  const auto localPosRowVec = (pos - rframeOrigin).transpose();
  // Initialize the derivative of local 3D cartesian coordinates w.r.t.
  // alignment parameters (without path correction)
  AlignmentToCartesian3DMatrix alignToCartesian =
      AlignmentToCartesian3DMatrix::Zero();
  //@Todo: add enum for cartesian coordinates
  alignToCartesian.block<1, 3>(ePos_X, ePos_X) = -localXAxis.transpose();
  alignToCartesian.block<1, 3>(ePos_X, eRot_X) =
      localPosRowVec * rotToLocalXAxisDer;
  alignToCartesian.block<1, 3>(ePos_Y, ePos_X) = -localYAxis.transpose();
  alignToCartesian.block<1, 3>(ePos_Y, eRot_X) =
      localPosRowVec * rotToLocalYAxisDer;
  alignToCartesian.block<1, 3>(ePos_Z, ePos_X) = -localZAxis.transpose();
  alignToCartesian.block<1, 3>(ePos_Z, eRot_X) =
      localPosRowVec * rotToLocalZAxisDer;

  // Cosine of angle between momentum direction and z axis of measurement frame
  const double cosThetaDir = localZAxis * dir;
  // Derivative of propagation path w.r.t. reference frame origin and rotation
  const auto posToPathDerRowVec = localZAxis.transpose() / cosThetaDir;
  const auto rotToPathDerRowVec =
      -localPosRowVec * rotToLocalZAxisDer / cosThetaDir;

  // Initialize the derivative of bound parameters w.r.t. alignment parameters
  AlignmentToBoundMatrix alignToLocal = AlignmentToBoundMatrix::Zero();
  // Derivative of eLOC_0
  const auto jacToLoc0 = jacToLocal.block<1, eFreeParametersSize>(0, 0);
  alignToLocal.block<1, 3>(eLOC_0, ePos_X) =
      cartesianToLocal.block<1, 3>(eLOC_0, ePos_X) *
          alignToCartesian.block<3, 3>(ePos_X, ePos_X) +
      jacToLoc0 * derivatives * posToPathDerRowVec;
  alignToLocal.block<1, 3>(eLOC_0, eRot_X) =
      cartesianToLocal.block<1, 3>(eLOC_0, ePos_X) *
          alignToCartesian.block<3, 3>(ePos_X, eRot_X) +
      jacToLoc0 * derivatives * rotToPathDerRowVec;
  // Derivative of eLOC_1
  const auto jacToLoc1 = jacToLocal.block<1, eFreeParametersSize>(1, 0);
  alignToLocal.block<1, 3>(eLOC_1, ePos_X) =
      cartesianToLocal.block<1, 3>(eLOC_1, ePos_X) *
          alignToCartesian.block<3, 3>(ePos_X, ePos_X) +
      jacToLoc1 * derivatives * posToPathDerRowVec;
  alignToLocal.block<1, 3>(eLOC_1, eRot_X) =
      cartesianToLocal.block<1, 3>(eLOC_1, ePos_X) *
          alignToCartesian.block<3, 3>(ePos_X, eRot_X) +
      jacToLoc1 * derivatives * rotToPathDerRowVec;
  // Derivative of ePHI
  const auto jacToPhi = jacToLocal.block<1, eFreeParametersSize>(2, 0);
  alignToLocal.block<1, 3>(ePHI, ePos_X) =
      jacToPhi * derivatives * posToPathDerRowVec;
  alignToLocal.block<1, 3>(ePHI, eRot_X) =
      jacToPhi * derivatives * rotToPathDerRowVec;
  // Derivative of eTHETA
  const auto jacToTheta = jacToLocal.block<1, eFreeParametersSize>(3, 0);
  alignToLocal.block<1, 3>(eTHETA, ePos_X) =
      jacToTheta * derivatives * posToPathDerRowVec;
  alignToLocal.block<1, 3>(eTHETA, eRot_X) =
      jacToTheta * derivatives * rotToPathDerRowVec;
  // Derivative of eQOP
  const auto jacToQoP = jacToLocal.block<1, eFreeParametersSize>(4, 0);
  alignToLocal.block<1, 3>(eQOP, ePos_X) =
      jacToQoP * derivatives * posToPathDerRowVec;
  alignToLocal.block<1, 3>(eQOP, eRot_X) =
      jacToQoP * derivatives * rotToPathDerRowVec;
  // Derivative of eT
  const auto jacToT = jacToLocal.block<1, eFreeParametersSize>(5, 0);
  alignToLocal.block<1, 3>(eT, ePos_X) =
      jacToT * derivatives * posToPathDerRowVec;
  alignToLocal.block<1, 3>(eT, eRot_X) =
      jacToT * derivatives * rotToPathDerRowVec;

  return alignToLocal;
}

}  // namespace detail
}  // namespace Acts
