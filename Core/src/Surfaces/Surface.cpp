// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/detail/AlignmentHelper.hpp"

#include <iomanip>
#include <iostream>
#include <utility>

Acts::Surface::Surface(std::shared_ptr<const Transform3D> tform)
    : GeometryObject(), m_transform(std::move(tform)) {}

Acts::Surface::Surface(const DetectorElementBase& detelement)
    : GeometryObject(),
      m_transform(nullptr),
      m_associatedDetElement(&detelement) {}

Acts::Surface::Surface(const Surface& other)
    : GeometryObject(other),
      std::enable_shared_from_this<Surface>(),
      m_transform(other.m_transform),
      m_surfaceMaterial(other.m_surfaceMaterial) {}

Acts::Surface::Surface(const GeometryContext& gctx, const Surface& other,
                       const Transform3D& shift)
    : GeometryObject(),
      m_transform(std::make_shared<const Transform3D>(
          Transform3D(shift * other.transform(gctx)))),
      m_associatedLayer(nullptr),
      m_surfaceMaterial(other.m_surfaceMaterial) {}

Acts::Surface::~Surface() = default;

bool Acts::Surface::isOnSurface(const GeometryContext& gctx,
                                const Vector3D& position,
                                const Vector3D& momentum,
                                const BoundaryCheck& bcheck) const {
  // create the local position
  Vector2D lposition{0., 0.};
  // global to local transformation
  bool gtlSuccess = globalToLocal(gctx, position, momentum, lposition);
  if (gtlSuccess) {
    return bcheck ? bounds().inside(lposition, bcheck) : true;
  }
  // did not succeed
  return false;
}

const Acts::AlignmentToBoundMatrix Acts::Surface::alignmentToBoundDerivative(
    const GeometryContext& gctx, const FreeVector& derivatives,
    const Vector3D& position, const Vector3D& direction) const {
  // The vector between position and center
  const ActsRowVector<double, 3> pcRowVec =
      (position - center(gctx)).transpose();
  // The local frame rotation
  const auto& rotation = transform(gctx).rotation();
  // The axes of local frame
  const Vector3D localXAxis = rotation.col(0);
  const Vector3D localYAxis = rotation.col(1);
  const Vector3D localZAxis = rotation.col(2);

  // 1) Calcuate the derivative of local frame axes w.r.t its rotation
  const auto& [rotToLocalXAxis, rotToLocalYAxis, rotToLocalZAxis] =
      detail::rotationToLocalAxesDerivative(rotation);
  // 2) Calculate the derivative of local 3D Cartesian coordinates w.r.t.
  // alignment parameters (without path correction)
  AlignmentToLocalCartesianMatrix alignToLoc3D =
      AlignmentToLocalCartesianMatrix::Zero();
  alignToLoc3D.block<1, 3>(eX, eAlignmentCenter0) = -localXAxis.transpose();
  alignToLoc3D.block<1, 3>(eY, eAlignmentCenter0) = -localYAxis.transpose();
  alignToLoc3D.block<1, 3>(eZ, eAlignmentCenter0) = -localZAxis.transpose();
  alignToLoc3D.block<1, 3>(eX, eAlignmentRotation0) =
      pcRowVec * rotToLocalXAxis;
  alignToLoc3D.block<1, 3>(eY, eAlignmentRotation0) =
      pcRowVec * rotToLocalYAxis;
  alignToLoc3D.block<1, 3>(eZ, eAlignmentRotation0) =
      pcRowVec * rotToLocalZAxis;
  // 3) Calculate the derivative of track position represented in
  // (local) bound track parameters (could be in non-Cartesian coordinates)
  // w.r.t. track position represented in local 3D Cartesian coordinates.
  const auto& loc3DToLocBound =
      localCartesianToBoundLocalDerivative(gctx, position);
  // 4) Calculate the derivative of path length w.r.t. alignment parameters
  const auto& alignToPath =
      alignmentToPathDerivative(gctx, rotToLocalZAxis, position, direction);
  // 5) Calculate the jacobian from free parameters to bound parameters
  FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
  initJacobianToLocal(gctx, jacToLocal, position, direction);
  // 6) Initialize the derivative of bound parameters w.r.t. alignment
  // parameters
  AlignmentToBoundMatrix alignToBound = AlignmentToBoundMatrix::Zero();
  // -> For bound track parameters eLOC_0, eLOC_1, it's
  // loc3DToLocBound*alignToLoc3D +
  // jacToLocal*derivatives*alignToPath
  alignToBound.block<2, eAlignmentParametersSize>(eLOC_0, eAlignmentCenter0) =
      loc3DToLocBound * alignToLoc3D +
      jacToLocal.block<2, eFreeParametersSize>(eLOC_0, eFreePos0) *
          derivatives * alignToPath;
  // -> For bound track parameters ePHI, eTHETA, eQOP, eT, it's
  // jacToLocal*derivatives*alignToPath
  alignToBound.block<4, eAlignmentParametersSize>(ePHI, eAlignmentCenter0) =
      jacToLocal.block<4, eFreeParametersSize>(ePHI, eFreePos0) * derivatives *
      alignToPath;

  return alignToBound;
}

const Acts::AlignmentRowVector Acts::Surface::alignmentToPathDerivative(
    const GeometryContext& gctx, const RotationMatrix3D& rotToLocalZAxis,
    const Vector3D& position, const Vector3D& direction) const {
  // The vector between position and center
  const ActsRowVector<double, 3> pcRowVec =
      (position - center(gctx)).transpose();
  // The local frame rotation
  const auto& rotation = transform(gctx).rotation();
  // The local frame z axis
  const Vector3D localZAxis = rotation.col(2);

  // Cosine of angle between momentum direction and local frame z axis
  const double dirZ = localZAxis.dot(direction);
  // Initialize the derivative of propagation path w.r.t. local frame
  // translation (origin) and rotation
  AlignmentRowVector alignToPath = AlignmentRowVector::Zero();
  alignToPath.segment<3>(eAlignmentCenter0) = localZAxis.transpose() / dirZ;
  alignToPath.segment<3>(eAlignmentRotation0) =
      -pcRowVec * rotToLocalZAxis / dirZ;

  return alignToPath;
}

std::shared_ptr<Acts::Surface> Acts::Surface::getSharedPtr() {
  return shared_from_this();
}

std::shared_ptr<const Acts::Surface> Acts::Surface::getSharedPtr() const {
  return shared_from_this();
}

Acts::Surface& Acts::Surface::operator=(const Surface& other) {
  if (&other != this) {
    GeometryObject::operator=(other);
    // detector element, identifier & layer association are unique
    m_transform = other.m_transform;
    m_associatedLayer = other.m_associatedLayer;
    m_surfaceMaterial = other.m_surfaceMaterial;
    m_associatedDetElement = other.m_associatedDetElement;
  }
  return *this;
}

bool Acts::Surface::operator==(const Surface& other) const {
  // (a) fast exit for pointer comparison
  if (&other == this) {
    return true;
  }
  // (b) fast exit for type
  if (other.type() != type()) {
    return false;
  }
  // (c) fast exit for bounds
  if (other.bounds() != bounds()) {
    return false;
  }
  // (d) compare  detector elements
  if (m_associatedDetElement != other.m_associatedDetElement) {
    return false;
  }
  // (e) compare transform values
  if (m_transform != nullptr && other.m_transform != nullptr) {
    if (!m_transform->isApprox(*other.m_transform, 1e-9)) {
      return false;
    }
  }
  // (f) compare material
  if (m_surfaceMaterial != other.m_surfaceMaterial) {
    return false;
  }

  // we should be good
  return true;
}

// overload dump for stream operator
std::ostream& Acts::Surface::toStream(const GeometryContext& gctx,
                                      std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(4);
  sl << name() << std::endl;
  const Vector3D& sfcenter = center(gctx);
  sl << "     Center position  (x, y, z) = (" << sfcenter.x() << ", "
     << sfcenter.y() << ", " << sfcenter.z() << ")" << std::endl;
  Acts::RotationMatrix3D rot(transform(gctx).matrix().block<3, 3>(0, 0));
  Acts::Vector3D rotX(rot.col(0));
  Acts::Vector3D rotY(rot.col(1));
  Acts::Vector3D rotZ(rot.col(2));
  sl << std::setprecision(6);
  sl << "     Rotation:             colX = (" << rotX(0) << ", " << rotX(1)
     << ", " << rotX(2) << ")" << std::endl;
  sl << "                           colY = (" << rotY(0) << ", " << rotY(1)
     << ", " << rotY(2) << ")" << std::endl;
  sl << "                           colZ = (" << rotZ(0) << ", " << rotZ(1)
     << ", " << rotZ(2) << ")" << std::endl;
  sl << "     Bounds  : " << bounds();
  sl << std::setprecision(-1);
  return sl;
}

bool Acts::Surface::operator!=(const Acts::Surface& sf) const {
  return !(operator==(sf));
}
