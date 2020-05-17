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
    const FreeVector& derivatives,
    CartesianToBoundLocalMatrix& locCartesianToLocBound) {
  // The position and momentum direction
  const auto pos = boundParams.position();
  const auto dir = boundParams.momentum().normalized();
  // The reference surface
  const auto surface = &boundParams.referenceSurface();
  // The local frame transform
  const auto& sTransform = surface->transform(gctx);
  const auto& rotation = sTransform.rotation();

  // Calcuate the derivative of local frame axes w.r.t its rotation
  const auto& rotToAxes = rotationToLocalAxesDerivative(rotation);
  // Calculate the derivative of local 3D Cartesian coordinates w.r.t.
  // alignment parameters (without path correction)
  const auto& alignToLocCartesian =
      alignmentToLocalCartesianDerivative(sTransform, rotToAxes, pos);
  // Calculate the derivative of propagation path w.r.t. alignment parameters
  // @Todo: need re-implementation for LineSurface, CylinderSurface and
  // ConeSurface?
  const auto& alignToPath =
      alignmentToPathDerivative(sTransform, rotToAxes, pos, dir);

  // Initialize the jacobian from free parameters to bound parameters
  FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
  // Set the jacobian to local
  surface->initJacobianToLocal(geoContext, jacToLocal, pos, dir);

  // Initialize the derivative of bound parameters w.r.t. alignment parameters
  AlignmentToBoundMatrix alignToBound = AlignmentToBoundMatrix::Zero();
  const auto jacToLoc0 =
      jacToLocal.block<1, eFreeParametersSize>(eLOC_0, eFreePos0);
  const auto jacToLoc1 =
      jacToLocal.block<1, eFreeParametersSize>(eLOC_1, eFreePos0);
  alignToBound.block<1, eAlignmentParametersSize>(eLOC_0, eCenter_X) =
      locCartesianToLocBound.block<1, 3>(eLOC_0, eX) * alignToLocCartesian +
      jacToLoc0 * derivatives * alignToPath;
  alignToBound.block<1, eAlignmentParametersSize>(eLOC_1, eCenter_X) =
      locCartesianToLocBound.block<1, 3>(eLOC_1, eX) * alignToLocCartesian +
      jacToLoc1 * derivatives * alignToPath;
  alignToBound.block<4, eAlignmentParametersSize>(ePHI, eCenter_X) =
      jacToLocal.block<4, eFreeParametersSize>(ePHI, eFreePos0) * derivatives *
      alignToPath;

  return alignToBound;
}

AlignmentToBoundMatrix layerAlignmentToBoundDerivative(
    const GeometryContext& gctx, const BoundParameters& boundParams,
    const FreeVector& derivatives,
    const AlignmentMatrix& layerAlignToSurfaceAlign,
    const CartesianToBoundLocalMatrix& locCartesianToLocBound) {
  const auto surfaceAlignToBound = surfaceAlignmentToBoundDerivative(
      gctx, boundParams, derivatives, locCartesianToLocBound);
  return surfaceAlignToBound * layerAlignToSurfaceAlign;
}

AlignmentToBoundMatrix volumeAlignmentToBoundDerivative(
    const GeometryContext& gctx, const BoundParameters& boundParams,
    const FreeVector& derivatives,
    const AlignmentMatrix& volumeAlignToSurfaceAlign,
    const CartesianToBoundLocalMatrix& locCartesianToLocBound) {
  const auto surfaceAlignToBound = surfaceAlignmentToBoundDerivative(
      gctx, boundParams, derivatives, locCartesianToLocBound);
  return surfaceAlignToBound * volumeAlignToSurfaceAlign;
}

AlignmentToCartesianMatrix alignmentToLocalCartesianDerivative(
    const Transform3D& sTransform, const RotationToAxes& rotToAxes,
    const Vector3D& position) {
  const auto& center = sTransform.translation();
  const auto& rotation = sTransform.rotation();
  // The vector between position and local frame origin
  const auto localPosRowVec = (position - center).transpose();
  // The axes of local frame
  const auto localXAxis = rotation.col(0);
  const auto localYAxis = rotation.col(1);
  const auto localZAxis = rotation.col(2);

  // Get the derivative of local axes w.r.t. local axes rotation
  const auto& [rotToLocalXAxis, rotToLocalYAxis, rotToLocalZAxis] = rotToAxes;

  // Initialize the derivative of local 3D Cartesian coordinates w.r.t.
  // alignment parameters (without path correction)
  AlignmentToCartesianMatrix alignToLocCartesian =
      AlignmentToCartesianMatrix::Zero();
  // Derivative of local Cartesian 3D coordinates w.r.t. local frame translation
  alignToLocCartesian.block<1, 3>(eX, eCenter_X) = -localXAxis.transpose();
  alignToLocCartesian.block<1, 3>(eY, eCenter_X) = -localYAxis.transpose();
  alignToLocCartesian.block<1, 3>(eZ, eCenter_X) = -localZAxis.transpose();
  // Derivative of local Cartesian 3D coordinates w.r.t. local frame
  // rotation
  alignToLocCartesian.block<1, 3>(eX, eRotation_X) =
      localPosRowVec * rotToLocalXAxis;
  alignToLocCartesian.block<1, 3>(eY, eRotation_X) =
      localPosRowVec * rotToLocalYAxis;
  alignToLocCartesian.block<1, 3>(eZ, eRotation_X) =
      localPosRowVec * rotToLocalZAxis;

  return alignToLocCartesian;
}

AlignmentRowVector alignToPathDerivative(const Transform3D& sTransform,
                                         const RotationToAxes& rotToAxes,
                                         const Vector3D& position,
                                         const Vector3D& direction) {
  const auto& center = sTransform.translation();
  const auto& rotation = sTransform.rotation();
  // The vector between position and local frame origin
  const auto localPosRowVec = (position - center).transpose();
  // The axes of local frame
  const auto localXAxis = rotation.col(0);
  const auto localYAxis = rotation.col(1);
  const auto localZAxis = rotation.col(2);

  // Get the derivative of local axes w.r.t. local axes rotation
  const auto& [rotToLocalXAxis, rotToLocalYAxis, rotToLocalZAxis] = rotToAxes;

  // Cosine of angle between momentum direction and local frame z axis
  const double cosThetaDir = localZAxis.transpose() * direction;
  // Initialize the derivative of propagation path w.r.t. local frame
  // translation (origin) and rotation
  AlignmentRowVector alignToPath = AlignmentRowVector::Zero();
  alignToPath.segment<3>(eCenter_X) = localZAxis.transpose() / cosThetaDir;
  alignToPath.segment<3>(eRotation_X) =
      -localPosRowVec * rotToLocalZAxis / cosThetaDir;

  return alignToPath;
}

RotationToAxes rotationToLocalAxesDerivative(const RotationMatrix3D& rotation) {
  // Get Euler angles for rotation representated by rotZ * rotY * rotX, i.e.
  // first rotation around x axis, then y axis, last z axis
  // The elements stored in rotAngles is (rotZ, rotY, rotX)
  const Vector3D rotAngles = rotation.eulerAngles(2, 1, 0);
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
