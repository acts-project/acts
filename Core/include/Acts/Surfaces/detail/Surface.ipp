// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/detail/TransformationBoundToFree.hpp"

inline Acts::Vector3D Acts::Surface::center(const GeometryContext& gctx) const {
  // fast access via tranform matrix (and not translation())
  auto tMatrix = transform(gctx).matrix();
  return Vector3D(tMatrix(0, 3), tMatrix(1, 3), tMatrix(2, 3));
}

inline Acts::Vector3D Acts::Surface::normal(const GeometryContext& gctx,
                                            const Vector3D& /*unused*/) const {
  return normal(gctx, s_origin2D);
}

inline const Acts::Transform3D& Acts::Surface::transform(
    const GeometryContext& gctx) const {
  if (m_associatedDetElement != nullptr) {
    return m_associatedDetElement->transform(gctx);
  }
  return m_transform;
}

inline bool Acts::Surface::insideBounds(const Vector2D& lposition,
                                        const BoundaryCheck& bcheck) const {
  return bounds().inside(lposition, bcheck);
}

inline Acts::RotationMatrix3D Acts::Surface::referenceFrame(
    const GeometryContext& gctx, const Vector3D& /*unused*/,
    const Vector3D& /*unused*/) const {
  return transform(gctx).matrix().block<3, 3>(0, 0);
}

inline Acts::BoundToFreeMatrix Acts::Surface::jacobianLocalToGlobal(
    const GeometryContext& gctx, const BoundVector& boundParams) const {
  // Transform from bound to free parameters
  FreeVector freeParams =
      detail::transformBoundToFreeParameters(*this, gctx, boundParams);
  // The global position
  const Vector3D position = freeParams.segment<3>(eFreePos0);
  // The direction
  const Vector3D direction = freeParams.segment<3>(eFreeDir0);
  // Get the sines and cosines directly
  const double cos_theta = std::cos(boundParams[eBoundTheta]);
  const double sin_theta = std::sin(boundParams[eBoundTheta]);
  const double cos_phi = std::cos(boundParams[eBoundPhi]);
  const double sin_phi = std::sin(boundParams[eBoundPhi]);
  // retrieve the reference frame
  const auto rframe = referenceFrame(gctx, position, direction);
  // Initialize the jacobian from local to global
  BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();
  // the local error components - given by reference frame
  jacToGlobal.topLeftCorner<3, 2>() = rframe.topLeftCorner<3, 2>();
  // the time component
  jacToGlobal(eFreeTime, eBoundTime) = 1;
  // the momentum components
  jacToGlobal(eFreeDir0, eBoundPhi) = (-sin_theta) * sin_phi;
  jacToGlobal(eFreeDir0, eBoundTheta) = cos_theta * cos_phi;
  jacToGlobal(eFreeDir1, eBoundPhi) = sin_theta * cos_phi;
  jacToGlobal(eFreeDir1, eBoundTheta) = cos_theta * sin_phi;
  jacToGlobal(eFreeDir2, eBoundTheta) = (-sin_theta);
  jacToGlobal(eFreeQOverP, eBoundQOverP) = 1;
  return jacToGlobal;
}

inline Acts::FreeToBoundMatrix Acts::Surface::jacobianGlobalToLocal(
    const GeometryContext& gctx, const FreeVector& parameters) const {
  // The global position
  const auto position = parameters.segment<3>(eFreePos0);
  // The direction
  const auto direction = parameters.segment<3>(eFreeDir0);
  // Optimized trigonometry on the propagation direction
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;
  // The measurement frame of the surface
  RotationMatrix3D rframeT =
      referenceFrame(gctx, position, direction).transpose();
  // Initalize the jacobian from global to local
  FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
  // Local position component given by the refernece frame
  jacToLocal.block<2, 3>(eBoundLoc0, eFreePos0) = rframeT.block<2, 3>(0, 0);
  // Time component
  jacToLocal(eBoundTime, eFreeTime) = 1;
  // Directional and momentum elements for reference frame surface
  jacToLocal(eBoundPhi, eFreeDir0) = -sinPhi * invSinTheta;
  jacToLocal(eBoundPhi, eFreeDir1) = cosPhi * invSinTheta;
  jacToLocal(eBoundTheta, eFreeDir0) = cosPhi * cosTheta;
  jacToLocal(eBoundTheta, eFreeDir1) = sinPhi * cosTheta;
  jacToLocal(eBoundTheta, eFreeDir2) = -sinTheta;
  jacToLocal(eBoundQOverP, eFreeQOverP) = 1;
  return jacToLocal;
}

inline Acts::FreeRowVector Acts::Surface::freeToPathDerivative(
    const GeometryContext& gctx, const FreeVector& parameters) const {
  // The global position
  const auto position = parameters.segment<3>(eFreePos0);
  // The direction
  const auto direction = parameters.segment<3>(eFreeDir0);
  // The measurement frame of the surface
  const RotationMatrix3D rframe = referenceFrame(gctx, position, direction);
  // The measurement frame z axis
  const Vector3D refZAxis = rframe.col(2);
  // Cosine of angle between momentum direction and measurement frame z axis
  const double dz = refZAxis.dot(direction);
  // Initialize the derivative
  FreeRowVector freeToPath = FreeRowVector::Zero();
  freeToPath.segment<3>(eFreePos0) = -1.0 * refZAxis.transpose() / dz;
  return freeToPath;
}

inline const Acts::DetectorElementBase*
Acts::Surface::associatedDetectorElement() const {
  return m_associatedDetElement;
}

inline const Acts::Layer* Acts::Surface::associatedLayer() const {
  return (m_associatedLayer);
}

inline const Acts::ISurfaceMaterial* Acts::Surface::surfaceMaterial() const {
  return m_surfaceMaterial.get();
}

inline const std::shared_ptr<const Acts::ISurfaceMaterial>&
Acts::Surface::surfaceMaterialSharedPtr() const {
  return m_surfaceMaterial;
}

inline void Acts::Surface::assignSurfaceMaterial(
    std::shared_ptr<const Acts::ISurfaceMaterial> material) {
  m_surfaceMaterial = std::move(material);
}

inline void Acts::Surface::associateLayer(const Acts::Layer& lay) {
  m_associatedLayer = (&lay);
}
