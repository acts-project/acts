// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/Plugins/Identification/Identifier.hpp"

namespace Acts {
class Surface;
class PlanarBounds;
class DiscBounds;
class ISurfaceMaterial;
class DigitizationModule;
}  // namespace Acts

namespace ActsExamples {

namespace Generic {

/// @class GenericDetectorElement
///
/// This is a lightweight type of detector element,
/// it simply implements the base class.
///
class GenericDetectorElement : public Acts::IdentifiedDetectorElement {
 public:
  /// Broadcast the ContextType
  using ContextType = Acts::GeometryContext;

  /// Constructor for single sided detector element
  /// - bound to a Plane Surface
  ///
  /// @param identifier is the module identifier
  /// @param transform is the transform that element the layer in 3D frame
  /// @param pBounds is the planar bounds for the planar detector element
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  GenericDetectorElement(
      const Identifier identifier,
      std::shared_ptr<const Acts::Transform3> transform,
      std::shared_ptr<const Acts::PlanarBounds> pBounds, double thickness,
      std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr,
      std::shared_ptr<const Acts::DigitizationModule> digitizationModule =
          nullptr);

  /// Constructor for single sided detector element
  /// - bound to a Disc Surface
  ///
  /// @param identifier is the module identifier
  /// @param transform is the transform that element the layer in 3D frame
  /// @param dBounds is the planar bounds for the disc like detector element
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  GenericDetectorElement(
      const Identifier identifier,
      std::shared_ptr<const Acts::Transform3> transform,
      std::shared_ptr<const Acts::DiscBounds> dBounds, double thickness,
      std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr,
      std::shared_ptr<const Acts::DigitizationModule> digitizationModule =
          nullptr);

  /// Identifier
  Identifier identifier() const final;

  /// Return local to global transform associated with this detector element
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  ///
  /// @note this is called from the surface().transform(gctx) in the PROXY
  /// mode
  const Acts::Transform3& transform(
      const Acts::GeometryContext& gctx) const override;

  /// Return surface associated with this detector element
  const Acts::Surface& surface() const override;

  /// Set the identifier after construction (sometimes needed)
  void assignIdentifier(const Identifier& identifier);

  /// The maximal thickness of the detector element wrt normal axis
  double thickness() const override;

  /// Retrieve the DigitizationModule
  const std::shared_ptr<const Acts::DigitizationModule> digitizationModule()
      const override;

 private:
  /// the element representation
  /// identifier
  Identifier m_elementIdentifier;
  /// the transform for positioning in 3D space
  std::shared_ptr<const Acts::Transform3> m_elementTransform;
  /// the surface represented by it
  std::shared_ptr<const Acts::Surface> m_elementSurface;
  /// the element thickness
  double m_elementThickness;
  /// store either
  std::shared_ptr<const Acts::PlanarBounds> m_elementPlanarBounds = nullptr;
  std::shared_ptr<const Acts::DiscBounds> m_elementDiscBounds = nullptr;
  /// The Digitization module
  std::shared_ptr<const Acts::DigitizationModule> m_digitizationModule =
      nullptr;
};

inline void ActsExamples::Generic::GenericDetectorElement::assignIdentifier(
    const Identifier& identifier) {
  m_elementIdentifier = identifier;
}

inline Identifier ActsExamples::Generic::GenericDetectorElement::identifier()
    const {
  return m_elementIdentifier;
}

inline const Acts::Transform3&
ActsExamples::Generic::GenericDetectorElement::transform(
    const Acts::GeometryContext& /*gctx*/) const {
  return *m_elementTransform;
}

inline const Acts::Surface&
ActsExamples::Generic::GenericDetectorElement::surface() const {
  return *m_elementSurface;
}

inline double ActsExamples::Generic::GenericDetectorElement::thickness() const {
  return m_elementThickness;
}

inline const std::shared_ptr<const Acts::DigitizationModule>
ActsExamples::Generic::GenericDetectorElement::digitizationModule() const {
  return m_digitizationModule;
}

}  // end of namespace Generic

}  // end of namespace ActsExamples
