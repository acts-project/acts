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
AlignmentToBoundMatrix surfaceAlignmentToBoundDerivative(
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
  const auto& [rotToLocalXAxis, rotToLocalYAxis, rotToLocalZAxis] =
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
  alignmentToPath.segment<3>(eOrigin_X) = localZAxis.transpose() / cosThetaDir;
  alignmentToPath.segment<3>(eRotation_X) =
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

AlignmentToBoundMatrix layerAlignmentToBoundDerivative(
    const GeometryContext& gctx, const BoundParameters& boundParams,
    const FreeVector& derivatives, const Vector3D& rframeOrigin,
    const AlignmentMatrix& layerAlignToSurfaceAlign,
    const CartesianToBoundLocalMatrix& locCartesianToLocBound =
        CartesianToBoundLocalMatrix::Identity(2,
                                              eCartesianCoordinatesDimension)) {
  const auto surfaceAlignToBound = surfaceAlignmentToBoundDerivative(
      gctx, boundParams, derivatives, rframeOrigin, locCartesianToLocBound);
  return surfaceAlignToBound * layerAlignToSurfaceAlign;
}

AlignmentToBoundMatrix volumeAlignmentToBoundDerivative(
    const GeometryContext& gctx, const BoundParameters& boundParams,
    const FreeVector& derivatives, const Vector3D& rframeOrigin,
    const AlignmentMatrix& volumeAlignToSurfaceAlign,
    const CartesianToBoundLocalMatrix& locCartesianToLocBound =
        CartesianToBoundLocalMatrix::Identity(2,
                                              eCartesianCoordinatesDimension)) {
  const auto surfaceAlignToBound = surfaceAlignmentToBoundDerivative(
      gctx, boundParams, derivatives, rframeOrigin, locCartesianToLocBound);
  return surfaceAlignToBound * volumeAlignToSurfaceAlign;
}

std::tuple<RotationMatrix3D, RotationMatrix3D, RotationMatrix3D>
rotationToLocalAxesDerivative(const RotationMatrix3D& rframe) {
  // Get Euler angles for rotation representated by rotZ * rotY * rotX, i.e.
  // first rotation around x axis, then y axis, last z axis
  // The elements stored in rotAngles is (rotZ, rotY, rotX)
  const Vector3D rotAngles = rframe.eulerAngles(2, 1, 0);
  double sx = std::sin(rotAngles(2));
  double cx = std::cos(rotAngles(2));
  double sy = std::sin(rotAngles(1));
  double cy = std::cos(rotAngles(1));
  double sz = std::sin(rotAngles(0));
  double cz = std::cos(rotAngles(0));
  // rotZ * rotY * rotX =
  // [ cz*cy  cz*sy*sx-cx*sz  sz*sx+cz*cx*sy ]
  // [ cy*sz  cz*cx+sz*sy*sx  cx*sz*sy-cz*sx ]
  // [ -sy   cy*sx         cy*cx        ]

  // Derivative of local x axis w.r.t. (rotX, rotY, rotZ)
  RotationMatrix3D rotToLocalXAxis = RotationMatrix3D::Zero();
  rotToLocalXAxis.col(0) = Vector3D(0, 0, 0);
  rotToLocalXAxis.col(1) = Vector3D(-cz * sy, -sz * sy, -cy);
  rotToLocalXAxis.col(2) = Vector3D(-sz * cy, cz * cy, 0);
  // Derivative of local y axis w.r.t. (rotX, rotY, rotZ)
  RotationMatrix3D rotToLocalYAxis = RotationMatrix3D::Zero();
  rotToLocalYAxis.col(0) =
      Vector3D(cz * sy * cx + sz * sx, sz * sy * cx - cz * sx, cy * cx);
  rotToLocalYAxis.col(1) = Vector3D(cz * cy * sx, sz * cy * sx, -sy * sx);
  rotToLocalYAxis.col(2) =
      Vector3D(-sz * sy * sx - cz * cx, cz * sy * sx - sz * cx, 0);
  // Derivative of local z axis w.r.t. (rotX, rotY, rotZ)
  RotationMatrix3D rotToLocalZAxis = RotationMatrix3D::Zero();
  rotToLocalZAxis.col(0) =
      Vector3D(sz * cx - cz * sy * sx, -sz * sy * sx - cz * cx, -cy * sx);
  rotToLocalZAxis.col(1) = Vector3D(cz * cy * cx, sz * cy * cx, -sy * cx);
  rotToLocalZAxis.col(2) =
      Vector3D(cz * sx - sz * sy * cx, cz * sy * cx + sz * sx, 0);

  return std::make_tuple(std::move(rotToLocalXAxis), std::move(rotToLocalYAxis),
                         std::move(rotToLocalZAxis));
}

}  // namespace detail
}  // namespace Acts
