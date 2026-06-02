// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfacePlacementBase.hpp"
#include "ActsTests/CommonHelpers/LineSurfaceStub.hpp"

namespace Acts {
class PlanarBounds;
class CylinderBounds;
class DiscBounds;
class LineBounds;
class ISurfaceMaterial;
}  // namespace Acts

namespace ActsTests {

/// @class DetectorElementStub
///
/// This is a lightweight type of detector element,
/// it simply implements the base class.
class DetectorElementStub : public Acts::SurfacePlacementBase {
 public:
  /// Default constructor
  DetectorElementStub() = default;

  /// Constructor for single sided detector element
  /// - no surface bound to it
  /// @param transform places the element in global frame
  explicit DetectorElementStub(const Acts::Transform3& transform);

  /// Constructor for single sided detector element
  /// - bound to a Line Surface
  ///
  /// @param transform places the element in global frame
  /// @param cBounds is the cylindrical bounds
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  DetectorElementStub(
      const Acts::Transform3& transform,
      std::shared_ptr<const Acts::CylinderBounds> cBounds, double thickness,
      std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr);

  /// Constructor for single sided detector element
  /// - bound to a Plane Surface
  ///
  /// @param transform places the element in global frame
  /// @param pBounds is the planar bounds for the planar detector element
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  DetectorElementStub(
      const Acts::Transform3& transform,
      std::shared_ptr<const Acts::PlanarBounds> pBounds, double thickness,
      std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr);

  /// Constructor for single sided detector element
  /// - bound to a Line Surface
  ///
  /// @param transform places the element in global frame
  /// @param dBounds is the line bounds for the line like detector element
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  DetectorElementStub(
      const Acts::Transform3& transform,
      std::shared_ptr<const Acts::LineBounds> lBounds, double thickness,
      std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr);

  ///  Destructor
  ~DetectorElementStub() override = default;

  /// Return local to global transform associated with this identifier
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  ///
  /// @note this is called from the surface().transform() in the PROXY mode
  const Acts::Transform3& localToGlobalTransform(
      const Acts::GeometryContext& gctx) const override;

  /// Return surface associated with this detector element
  const Acts::Surface& surface() const override;

  /// Non-const access to surface associated with this detector element
  Acts::Surface& surface() override;

  /// The maximal thickness of the detector element wrt normal axis
  double thickness() const;

  /// Is the detector element a sensitive element
  bool isSensitive() const override { return true; }

 private:
  /// the transform for positioning in 3D space
  Acts::Transform3 m_elementTransform;
  /// the surface represented by it
  std::shared_ptr<Acts::Surface> m_elementSurface{nullptr};
  /// the element thickness
  double m_elementThickness{0.};
};

inline const Acts::Transform3& DetectorElementStub::localToGlobalTransform(
    const Acts::GeometryContext& /*gctx*/) const {
  return m_elementTransform;
}

inline const Acts::Surface& DetectorElementStub::surface() const {
  return *m_elementSurface;
}

inline Acts::Surface& DetectorElementStub::surface() {
  return *m_elementSurface;
}

inline double DetectorElementStub::thickness() const {
  return m_elementThickness;
}

}  // namespace ActsTests
