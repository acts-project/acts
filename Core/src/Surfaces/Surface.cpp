// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/Surface.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/detail/AlignmentHelper.hpp"
#include "Acts/Utilities/JacobianHelpers.hpp"
#include "Acts/Visualization/ViewConfig.hpp"

#include <iomanip>
#include <utility>

namespace Acts {

Surface::Surface(const Transform3& transform)
    : GeometryObject(), m_transform(std::make_unique<Transform3>(transform)) {}

Surface::Surface(const SurfacePlacementBase& placement) noexcept
    : GeometryObject(), m_placement(&placement) {}

Surface::Surface(const Surface& other) noexcept
    : GeometryObject(other),
      std::enable_shared_from_this<Surface>(),
      m_placement(other.m_placement),
      m_surfaceMaterial(other.m_surfaceMaterial) {
  if (other.m_transform) {
    m_transform = std::make_unique<Transform3>(*other.m_transform);
  }
}

Surface::Surface(const GeometryContext& gctx, const Surface& other,
                 const Transform3& shift) noexcept
    : GeometryObject(),
      m_transform(std::make_unique<Transform3>(
          shift * other.localToGlobalTransform(gctx))),
      m_surfaceMaterial(other.m_surfaceMaterial) {}

Surface::~Surface() noexcept = default;

std::ostream& operator<<(std::ostream& os, Surface::SurfaceType type) {
  return os << Surface::s_surfaceTypeNames[static_cast<std::size_t>(type)];
}

bool Surface::isOnSurface(const GeometryContext& gctx, const Vector3& position,
                          const Vector3& direction,
                          const BoundaryTolerance& boundaryTolerance,
                          double tolerance) const {
  // global to local transformation
  auto lpResult = globalToLocal(gctx, position, direction, tolerance);
  if (!lpResult.ok()) {
    return false;
  }
  return bounds().inside(lpResult.value(), boundaryTolerance);
}

AlignmentToBoundMatrix Surface::alignmentToBoundDerivative(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const FreeVector& pathDerivative) const {
  assert(isOnSurface(gctx, position, direction, BoundaryTolerance::Infinite()));

  // 1) Calculate the derivative of bound parameter local position w.r.t.
  // alignment parameters without path length correction
  const auto alignToBoundWithoutCorrection =
      alignmentToBoundDerivativeWithoutCorrection(gctx, position, direction);
  // 2) Calculate the derivative of path length w.r.t. alignment parameters
  const auto alignToPath = alignmentToPathDerivative(gctx, position, direction);
  // 3) Calculate the jacobian from free parameters to bound parameters
  FreeToBoundMatrix jacToLocal = freeToBoundJacobian(gctx, position, direction);
  // 4) The derivative of bound parameters w.r.t. alignment
  // parameters is alignToBoundWithoutCorrection +
  // jacToLocal*pathDerivative*alignToPath
  AlignmentToBoundMatrix alignToBound =
      alignToBoundWithoutCorrection + jacToLocal * pathDerivative * alignToPath;

  return alignToBound;
}

AlignmentToBoundMatrix Surface::alignmentToBoundDerivativeWithoutCorrection(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  static_cast<void>(direction);
  assert(isOnSurface(gctx, position, direction, BoundaryTolerance::Infinite()));

  // The vector between position and center
  const auto pcRowVec = (position - center(gctx)).transpose().eval();
  // The local frame rotation
  const auto& rotation = localToGlobalTransform(gctx).rotation();
  // The axes of local frame
  const auto& localXAxis = rotation.col(0);
  const auto& localYAxis = rotation.col(1);
  const auto& localZAxis = rotation.col(2);
  // Calculate the derivative of local frame axes w.r.t its rotation
  const auto [rotToLocalXAxis, rotToLocalYAxis, rotToLocalZAxis] =
      detail::rotationToLocalAxesDerivative(rotation);
  // Calculate the derivative of local 3D Cartesian coordinates w.r.t.
  // alignment parameters (without path correction)
  AlignmentToPositionMatrix alignToLoc3D = AlignmentToPositionMatrix::Zero();
  alignToLoc3D.block<1, 3>(eX, eAlignmentCenter0) = -localXAxis.transpose();
  alignToLoc3D.block<1, 3>(eY, eAlignmentCenter0) = -localYAxis.transpose();
  alignToLoc3D.block<1, 3>(eZ, eAlignmentCenter0) = -localZAxis.transpose();
  alignToLoc3D.block<1, 3>(eX, eAlignmentRotation0) =
      pcRowVec * rotToLocalXAxis;
  alignToLoc3D.block<1, 3>(eY, eAlignmentRotation0) =
      pcRowVec * rotToLocalYAxis;
  alignToLoc3D.block<1, 3>(eZ, eAlignmentRotation0) =
      pcRowVec * rotToLocalZAxis;
  // The derivative of bound local w.r.t. local 3D Cartesian coordinates
  Matrix<2, 3> loc3DToBoundLoc =
      localCartesianToBoundLocalDerivative(gctx, position);
  // Initialize the derivative of bound parameters w.r.t. alignment
  // parameters without path correction
  AlignmentToBoundMatrix alignToBound = AlignmentToBoundMatrix::Zero();
  // It's only relevant with the bound local position without path correction
  alignToBound.block<2, eAlignmentSize>(eBoundLoc0, eAlignmentCenter0) =
      loc3DToBoundLoc * alignToLoc3D;
  return alignToBound;
}

AlignmentToPathMatrix Surface::alignmentToPathDerivative(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  assert(isOnSurface(gctx, position, direction, BoundaryTolerance::Infinite()));

  // The vector between position and center
  const auto pcRowVec = (position - center(gctx)).transpose().eval();
  // The local frame rotation
  const auto& rotation = localToGlobalTransform(gctx).rotation();
  // The local frame z axis
  const auto& localZAxis = rotation.col(2);
  // Cosine of angle between momentum direction and local frame z axis
  const auto dz = localZAxis.dot(direction);
  // Calculate the derivative of local frame axes w.r.t its rotation
  const auto [rotToLocalXAxis, rotToLocalYAxis, rotToLocalZAxis] =
      detail::rotationToLocalAxesDerivative(rotation);
  // Initialize the derivative of propagation path w.r.t. local frame
  // translation (origin) and rotation
  AlignmentToPathMatrix alignToPath = AlignmentToPathMatrix::Zero();
  alignToPath.segment<3>(eAlignmentCenter0) = localZAxis.transpose() / dz;
  alignToPath.segment<3>(eAlignmentRotation0) =
      -pcRowVec * rotToLocalZAxis / dz;

  return alignToPath;
}

std::shared_ptr<Surface> Surface::getSharedPtr() {
  return shared_from_this();
}

std::shared_ptr<const Surface> Surface::getSharedPtr() const {
  return shared_from_this();
}

Surface& Surface::operator=(const Surface& other) {
  if (&other != this) {
    GeometryObject::operator=(other);
    // detector element, identifier & layer association are unique
    if (other.m_transform) {
      m_transform = std::make_unique<Transform3>(*other.m_transform);
    } else {
      m_transform.reset();
    }
    m_associatedLayer = other.m_associatedLayer;
    m_surfaceMaterial = other.m_surfaceMaterial;
    m_placement = other.m_placement;
    m_isSensitive = other.m_isSensitive;
  }
  return *this;
}

bool Surface::operator==(const Surface& other) const {
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
  if (m_placement != other.m_placement) {
    return false;
  }
  // (e) compare transform values
  if (m_transform && other.m_transform &&
      !m_transform->isApprox((*other.m_transform), 1e-9)) {
    return false;
  }
  // (f) compare material
  if (m_surfaceMaterial != other.m_surfaceMaterial) {
    return false;
  }
  // (g) compare sensitivity
  if (m_isSensitive != other.m_isSensitive) {
    return false;
  }

  // we should be good
  return true;
}

std::ostream& Surface::toStreamImpl(const GeometryContext& gctx,
                                    std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(4);
  sl << name() << std::endl;
  const Vector3& sfcenter = center(gctx);
  sl << "     Center position  (x, y, z) = (" << sfcenter.x() << ", "
     << sfcenter.y() << ", " << sfcenter.z() << ")" << std::endl;
  RotationMatrix3 rot(localToGlobalTransform(gctx).matrix().block<3, 3>(0, 0));
  Vector3 rotX(rot.col(0));
  Vector3 rotY(rot.col(1));
  Vector3 rotZ(rot.col(2));
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

std::string Surface::toString(const GeometryContext& gctx) const {
  std::stringstream ss;
  ss << toStream(gctx);
  return ss.str();
}

Vector3 Surface::center(const GeometryContext& gctx) const {
  // fast access via transform matrix (and not translation())
  auto tMatrix = localToGlobalTransform(gctx).matrix();
  return Vector3(tMatrix(0, 3), tMatrix(1, 3), tMatrix(2, 3));
}

const Transform3& Surface::transform(const GeometryContext& gctx) const {
  return localToGlobalTransform(gctx);
}

const Transform3& Surface::localToGlobalTransform(
    const GeometryContext& gctx) const {
  if (m_placement != nullptr) {
    return m_placement->localToGlobalTransform(gctx);
  }
  return *m_transform;
}

Vector2 Surface::closestPointOnBoundary(const Vector2& lposition,
                                        const SquareMatrix2& metric) const {
  return bounds().closestPoint(lposition, metric);
}

double Surface::distanceToBoundary(const Vector2& lposition) const {
  return bounds().distance(lposition);
}

bool Surface::insideBounds(const Vector2& lposition,
                           const BoundaryTolerance& boundaryTolerance) const {
  return bounds().inside(lposition, boundaryTolerance);
}

RotationMatrix3 Surface::referenceFrame(const GeometryContext& gctx,
                                        const Vector3& /*position*/,
                                        const Vector3& /*direction*/) const {
  return localToGlobalTransform(gctx).matrix().block<3, 3>(0, 0);
}

BoundToFreeMatrix Surface::boundToFreeJacobian(const GeometryContext& gctx,
                                               const Vector3& position,
                                               const Vector3& direction) const {
  assert(isOnSurface(gctx, position, direction, BoundaryTolerance::Infinite()));

  // retrieve the reference frame
  const auto rframe = referenceFrame(gctx, position, direction);

  // Initialize the jacobian from local to global
  BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();
  // the local error components - given by reference frame
  jacToGlobal.topLeftCorner<3, 2>() = rframe.topLeftCorner<3, 2>();
  // the time component
  jacToGlobal(eFreeTime, eBoundTime) = 1;
  // the momentum components
  jacToGlobal.block<3, 2>(eFreeDir0, eBoundPhi) =
      sphericalToFreeDirectionJacobian(direction);
  jacToGlobal(eFreeQOverP, eBoundQOverP) = 1;
  return jacToGlobal;
}

FreeToBoundMatrix Surface::freeToBoundJacobian(const GeometryContext& gctx,
                                               const Vector3& position,
                                               const Vector3& direction) const {
  assert(isOnSurface(gctx, position, direction, BoundaryTolerance::Infinite()));

  // The measurement frame of the surface
  RotationMatrix3 rframeT =
      referenceFrame(gctx, position, direction).transpose();

  // Initialize the jacobian from global to local
  FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
  // Local position component given by the reference frame
  jacToLocal.block<2, 3>(eBoundLoc0, eFreePos0) = rframeT.block<2, 3>(0, 0);
  // Time component
  jacToLocal(eBoundTime, eFreeTime) = 1;
  // Directional and momentum elements for reference frame surface
  jacToLocal.block<2, 3>(eBoundPhi, eFreeDir0) =
      freeToSphericalDirectionJacobian(direction);
  jacToLocal(eBoundQOverP, eFreeQOverP) = 1;
  return jacToLocal;
}

FreeToPathMatrix Surface::freeToPathDerivative(const GeometryContext& gctx,
                                               const Vector3& position,
                                               const Vector3& direction) const {
  assert(isOnSurface(gctx, position, direction, BoundaryTolerance::Infinite()));

  // The measurement frame of the surface
  const RotationMatrix3 rframe = referenceFrame(gctx, position, direction);
  // The measurement frame z axis
  const Vector3 refZAxis = rframe.col(2);
  // Cosine of angle between momentum direction and measurement frame z axis
  const double dz = refZAxis.dot(direction);
  // Initialize the derivative
  FreeToPathMatrix freeToPath = FreeToPathMatrix::Zero();
  freeToPath.segment<3>(eFreePos0) = -1.0 * refZAxis.transpose() / dz;
  return freeToPath;
}

const SurfacePlacementBase* Surface::surfacePlacement() const {
  return m_placement;
}

const DetectorElementBase* Surface::associatedDetectorElement() const {
  return dynamic_cast<const DetectorElementBase*>(m_placement);
}

double Surface::thickness() const {
  return m_thickness;
}

void Surface::assignThickness(double thick) {
  assert(thick >= 0.);
  m_thickness = thick;
}

const Layer* Surface::associatedLayer() const {
  return m_associatedLayer;
}

const ISurfaceMaterial* Surface::surfaceMaterial() const {
  return m_surfaceMaterial.get();
}

const std::shared_ptr<const ISurfaceMaterial>&
Surface::surfaceMaterialSharedPtr() const {
  return m_surfaceMaterial;
}

void Surface::assignDetectorElement(const SurfacePlacementBase& detelement) {
  assignSurfacePlacement(detelement);
}

void Surface::assignSurfacePlacement(const SurfacePlacementBase& placement) {
  m_placement = &placement;
  // resetting the transform as it will be handled through the detector element
  // now
  m_transform.reset();
  // reset sensitivity flag
  m_isSensitive = false;
}

void Surface::assignSurfaceMaterial(
    std::shared_ptr<const ISurfaceMaterial> material) {
  m_surfaceMaterial = std::move(material);
}

void Surface::associateLayer(const Layer& lay) {
  m_associatedLayer = (&lay);
}

void Surface::visualize(IVisualization3D& helper, const GeometryContext& gctx,
                        const ViewConfig& viewConfig) const {
  Polyhedron polyhedron =
      polyhedronRepresentation(gctx, viewConfig.quarterSegments);
  polyhedron.visualize(helper, viewConfig);
}

void Surface::assignIsSensitive(bool isSensitive) {
  if (m_placement != nullptr) {
    throw std::logic_error(
        "Cannot assign sensitivity to a surface associated to a detector "
        "element.");
  }
  m_isSensitive = isSensitive;
}

bool Surface::isSensitive() const {
  if (m_placement != nullptr) {
    return m_placement->isSensitive();
  }
  return m_isSensitive;
}
bool Surface::isAlignable() const {
  return m_placement != nullptr;
}

}  // namespace Acts
