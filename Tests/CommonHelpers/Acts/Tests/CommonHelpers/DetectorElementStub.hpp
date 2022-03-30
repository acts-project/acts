// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DetectorElementStub.h, Acts project, Generic Detector plugin
///////////////////////////////////////////////////////////////////

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/Identification/Identifier.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Tests/CommonHelpers/LineSurfaceStub.hpp"

namespace Acts {

class PlanarBounds;
class DiscBounds;
class ISurfaceMaterial;
class LineBounds;

namespace Test {

/// @class DetectorElementStub
///
/// This is a lightweight type of detector element,
/// it simply implements the base class.
class DetectorElementStub : public IdentifiedDetectorElement {
 public:
  DetectorElementStub() : IdentifiedDetectorElement() {}

  DetectorElementStub(const Transform3& transform)
      : IdentifiedDetectorElement(), m_elementTransform(std::move(transform)) {}

  /// Constructor for single sided detector element
  /// - bound to a Plane Surface
  ///
  /// @param transform is the transform that element the layer in 3D frame
  /// @param pBounds is the planar bounds for the planar detector element
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  DetectorElementStub(
      const Transform3& transform, std::shared_ptr<const PlanarBounds> pBounds,
      double thickness,
      std::shared_ptr<const ISurfaceMaterial> material = nullptr,
      std::shared_ptr<const Acts::DigitizationModule> digitizationModule =
          nullptr)
      : IdentifiedDetectorElement(),
        m_elementTransform(transform),
        m_elementThickness(thickness),
        m_digitizationModule(digitizationModule) {
    auto mutableSurface = Surface::makeShared<PlaneSurface>(pBounds, *this);
    mutableSurface->assignSurfaceMaterial(material);
    m_elementSurface = mutableSurface;
  }

  /// Constructor for single sided detector element
  /// - bound to a Line Surface
  ///
  /// @param transform is the transform that element the layer in 3D frame
  /// @param dBounds is the line bounds for the line like detector element
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  DetectorElementStub(
      const Transform3& transform, std::shared_ptr<const LineBounds> lBounds,
      double thickness,
      std::shared_ptr<const ISurfaceMaterial> material = nullptr,
      std::shared_ptr<const Acts::DigitizationModule> digitizationModule =
          nullptr)
      : IdentifiedDetectorElement(),
        m_elementTransform(transform),
        m_elementThickness(thickness),
        m_digitizationModule(digitizationModule) {
    auto mutableSurface = Surface::makeShared<LineSurfaceStub>(lBounds, *this);
    mutableSurface->assignSurfaceMaterial(material);
    m_elementSurface = mutableSurface;
  }

  ///  Destructor
  ~DetectorElementStub() override { /*nop */
  }

  /// Return local to global transform associated with this identifier
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  ///
  /// @note this is called from the surface().transform() in the PROXY mode
  const Transform3& transform(const GeometryContext& gctx) const override;

  /// Return surface associated with this detector element
  const Surface& surface() const override;

  /// The maximal thickness of the detector element wrt normal axis
  double thickness() const override;

  /// Identifier
  Identifier identifier() const override final;

  /// Retrieve the DigitizationModule
  const std::shared_ptr<const Acts::DigitizationModule> digitizationModule()
      const final override;

 private:
  /// the transform for positioning in 3D space
  Transform3 m_elementTransform;
  /// the surface represented by it
  std::shared_ptr<const Surface> m_elementSurface{nullptr};
  /// the element thickness
  double m_elementThickness{0.};
  /// identifier
  Identifier m_elementIdentifier;
  /// The Digitization module
  std::shared_ptr<const Acts::DigitizationModule> m_digitizationModule =
      nullptr;
};

inline const Transform3& DetectorElementStub::transform(
    const GeometryContext& /*gctx*/) const {
  return m_elementTransform;
}

inline const Surface& DetectorElementStub::surface() const {
  return *m_elementSurface;
}

inline double DetectorElementStub::thickness() const {
  return m_elementThickness;
}

inline Identifier DetectorElementStub::identifier() const {
  return m_elementIdentifier;
}

inline const std::shared_ptr<const Acts::DigitizationModule>
DetectorElementStub::digitizationModule() const {
  return m_digitizationModule;
}

}  // namespace Test
}  // namespace Acts
