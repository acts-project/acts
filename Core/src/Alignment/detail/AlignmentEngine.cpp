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
AlignmentToBoundMatrix alignmentToBoundDerivative(
    const GeometryContext& gctx, const BoundParameters& boundParams,
    const FreeVector& derivatives, const Vector3D& rframeOrigin,
    CartesianToBoundLocalMatrix& locCartesianToLocBound) {
  // The position and momentum direction
  const auto pos = boundParams.position();
  const auto dir = boundParams.momentum().normalized();
  // The reference surface
  const auto surface = &boundParams.referenceSurface();

  // Initialize the jacobian from free parameters to bound parameters
  FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
  // Set the jacobian to local, returns the transposed ref frame (i.e.
  // measurement frame)
  const auto rframeT =
      surface->initJacobianToLocal(geoContext, jacToLocal, pos, dir);
  const auto rframe = rframeT.transpose();
  // The axes of local measurement frame in Cartesian coordinates
  const auto localXAxis = rframe.col(0);
  const auto localYAxis = rframe.col(1);
  const auto localZAxis = rframe.col(2);
  // Get the derivative of local axes w.r.t. local axes rotation
  const auto [rotToLocalXAxis, rotToLocalYAxis, rotToLocalZAxis] =
      rotationToLocalAxesDerivative(rframe);
  // The vector of track position in reference frame
  const auto localPosRowVec = (pos - rframeOrigin).transpose();

  // Initialize the derivative of local 3D Cartesian coordinates w.r.t.
  // alignment parameters (without path correction)
  AlignmentToCartesianMatrix alignToLocCartesian =
      AlignmentToCartesianMatrix::Zero();
  // Derivative of local Cartesian 3D coordinates w.r.t. reference frame origin
  alignToLocCartesian.block<1, 3>(eX, eOrigin_X) = -localXAxis.transpose();
  alignToLocCartesian.block<1, 3>(eY, eOrigin_X) = -localYAxis.transpose();
  alignToLocCartesian.block<1, 3>(eZ, eOrigin_X) = -localZAxis.transpose();
  // Derivative of local Cartesian 3D coordinates w.r.t. reference frame
  // rotation
  alignToLocCartesian.block<1, 3>(eX, eRotation_X) =
      localPosRowVec * rotToLocalXAxis;
  alignToLocCartesian.block<1, 3>(eY, eRotation_X) =
      localPosRowVec * rotToLocalYAxis;
  alignToLocCartesian.block<1, 3>(eZ, eRotation_X) =
      localPosRowVec * rotToLocalZAxis;

  // Cosine of angle between momentum direction and z axis of measurement frame
  const double cosThetaDir = localZAxis.transpose() * dir;
  // Derivative of propagation path w.r.t. reference frame origin and rotation
  AlignmentRowVector alignmentToPath = AlignmentRowVector::Zero();
  alignmentToPath.block<1, 3>(0, eOrigin_X) =
      localZAxis.transpose() / cosThetaDir;
  alignmentToPath.block<1, 3>(0, eRotation_X) =
      -localPosRowVec * rotToLocalZAxis / cosThetaDir;

  // Initialize the derivative of bound parameters w.r.t. alignment parameters
  AlignmentToBoundMatrix alignToBound = AlignmentToBoundMatrix::Zero();
  const auto jacToLoc0 = jacToLocal.block<1, eFreeParametersSize>(0, 0);
  const auto jacToLoc1 = jacToLocal.block<1, eFreeParametersSize>(1, 0);
  const auto jacToPhi = jacToLocal.block<1, eFreeParametersSize>(2, 0);
  const auto jacToTheta = jacToLocal.block<1, eFreeParametersSize>(3, 0);
  const auto jacToQoP = jacToLocal.block<1, eFreeParametersSize>(4, 0);
  const auto jacToT = jacToLocal.block<1, eFreeParametersSize>(5, 0);
  alignToBound.block<1, 6>(eLOC_0, eOrigin_X) =
      locCartesianToLocBound.block<1, 3>(eLOC_0, eX) * alignToLocCartesian +
      jacToLoc0 * derivatives * alignmentToPath;
  alignToBound.block<1, 6>(eLOC_1, eOrigin_X) =
      locCartesianToLocBound.block<1, 3>(eLOC_1, eX) * alignToLocCartesian +
      jacToLoc1 * derivatives * alignmentToPath;
  alignToBound.block<1, 6>(ePHI, eOrigin_X) =
      jacToPhi * derivatives * alignmentToPath;
  alignToBound.block<1, 6>(eTHETA, eOrigin_X) =
      jacToTheta * derivatives * alignmentToPath;
  alignToBound.block<1, 6>(eQOP, eOrigin_X) =
      jacToQoP * derivatives * alignmentToPath;
  alignToBound.block<1, 6>(eT, eOrigin_X) =
      jacToT * derivatives * alignmentToPath;

  return alignToBound;
}

std::tuple<RotationMatrix3D, RotationMatrix3D, RotationMatrix3D>
rotationToLocalAxesDerivative(const RotationMatrix3D& rframe) {
  // Get Euler angles for rotation representated by Z1 * Y2 * X3, i.e.
  // first rotation around z axis, then y axis, last x axis
  const Vector3D rotAngles = rframe.eulerAngles(2, 1, 0);
  double s1 = std::sin(rotAngles(0));
  double c1 = std::cos(rotAngles(0));
  double s2 = std::sin(rotAngles(1));
  double c2 = std::cos(rotAngles(1));
  double s3 = std::sin(rotAngles(2));
  double c3 = std::cos(rotAngles(2));
  // Z1 * Y2 * X3 =
  // [ c1c2  c1s2s3-c3s1  s1s3+c1c3s2 ]
  // [ c2s1  c1c3+s1s2s3  c3s1s2-c1s3 ]
  // [ -s2   c2s3         c2c3        ]

  // Derivative of local x axis w.r.t. rotation around global x/y/z
  RotationMatrix3D rotToLocalXAxis = RotationMatrix3D::Zero();
  rotToLocalXAxis.col(0) = Vector3D(-s1 * c2, c1 * c2, 0);
  rotToLocalXAxis.col(1) = Vector3D(-c1 * s2, -s1 * s2, -c2);
  rotToLocalXAxis.col(2) = Vector3D(0, 0, 0);
  // Derivative of local y axis w.r.t. rotation around global x/y/z
  RotationMatrix3D rotToLocalYAxis = RotationMatrix3D::Zero();
  rotToLocalYAxis.col(0) =
      Vector3D(-s1 * s2 * s3 - c1 * c3, c1 * s2 * s3 - s1 * c3, 0);
  rotToLocalYAxis.col(1) = Vector3D(c1 * c2 * s3, s1 * c2 * s3, -s2 * s3);
  rotToLocalYAxis.col(2) =
      Vector3D(c1 * s2 * c3 + s1 * s3, s1 * s2 * c3 - c1 * s3, c2 * c3);
  // Derivative of local z axis w.r.t. rotation around global x/y/z
  RotationMatrix3D rotToLocalZAxis = RotationMatrix3D::Zero();
  rotToLocalZAxis.col(0) =
      Vector3D(c1 * s3 - s1 * s2 * c3, c1 * s2 * c3 + s1 * s3, 0);
  rotToLocalZAxis.col(1) = Vector3D(c1 * c2 * c3, s1 * c2 * c3, -s2 * c3);
  rotToLocalZAxis.col(2) =
      Vector3D(s1 * c3 - c1 * s2 * s3, -s1 * s2 * s3 - c1 * c3, -c2 * s3);

  return std::make_tuple(std::move(rotToLocalXAxis), std::move(rotToLocalYAxis),
                         std::move(rotToLocalZAxis));
}

}  // namespace detail
}  // namespace Acts
