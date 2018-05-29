// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DetectorElementStub.h, Acts project, Generic Detector plugin
///////////////////////////////////////////////////////////////////

#pragma once
// Algebra and Identifier
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Identifier.hpp"
// Geometry module
#include "Acts/Detector/DetectorElementBase.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "LineSurfaceStub.hpp"
//#include "Acts/Surfaces/DiscSurface.hpp"

namespace Acts {

class PlanarBounds;
class DiscBounds;
class SurfaceMaterial;
class LineBounds;

/// @class DetectorElementStub
///
/// This is a lightweight type of detector element,
/// it simply implements the base class.
///
class DetectorElementStub : public DetectorElementBase
{
public:
  DetectorElementStub() : DetectorElementBase() {}

  /// Constructor for single sided detector element
  /// - bound to a Plane Surface
  ///
  /// @param identifier is the module identifier
  /// @param transform is the transform that element the layer in 3D frame
  /// @param pBounds is the planar bounds for the planar detector element
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  DetectorElementStub(const Identifier                       identifier,
                      std::shared_ptr<const Transform3D>     transform,
                      std::shared_ptr<const PlanarBounds>    pBounds,
                      double                                 thickness,
                      std::shared_ptr<const SurfaceMaterial> material = nullptr)
    : DetectorElementBase()
    , m_elementIdentifier(std::move(identifier))
    , m_elementTransform(std::move(transform))
    , m_elementSurface(new PlaneSurface(pBounds, *this))
    , m_elementThickness(thickness)
    , m_elementSurfaces({m_elementSurface})
    , m_elementPlanarBounds(std::move(pBounds))
    , m_elementDiscBounds(nullptr)
    , m_elementLineBounds(nullptr)
  {
    auto mutableSurface = std::const_pointer_cast<Surface>(m_elementSurface);
    mutableSurface->setAssociatedMaterial(material);
  }

  /// Constructor for single sided detector element
  /// - bound to a Line Surface
  ///
  /// @param identifier is the module identifier
  /// @param transform is the transform that element the layer in 3D frame
  /// @param dBounds is the line bounds for the line like detector element
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  DetectorElementStub(const Identifier                       identifier,
                      std::shared_ptr<const Transform3D>     transform,
                      std::shared_ptr<const LineBounds>      lBounds,
                      double                                 thickness,
                      std::shared_ptr<const SurfaceMaterial> material = nullptr)
    : DetectorElementBase()
    , m_elementIdentifier(std::move(identifier))
    , m_elementTransform(std::move(transform))
    , m_elementSurface(new LineSurfaceStub(lBounds, *this))
    , m_elementThickness(thickness)
    , m_elementSurfaces({m_elementSurface})
    , m_elementPlanarBounds(nullptr)
    , m_elementDiscBounds(nullptr)
    , m_elementLineBounds(std::move(lBounds))
  {
    auto mutableSurface = std::const_pointer_cast<Surface>(m_elementSurface);
    mutableSurface->setAssociatedMaterial(material);
  }

  void
  assignIdentifier(const Identifier& identifier) override
  {
    m_elementIdentifier = identifier;
  }

  ///  Destructor
  ~DetectorElementStub() { /*nop */}

  /// Identifier
  Identifier
  identify() const override;

  /// Return local to global transform associated with this identifier
  ///
  /// @note this is called from the surface().transform() in the PROXY mode
  ///
  /// @param identifier is ignored for this simple detector element
  const Transform3D&
  transform(const Identifier& identifier = Identifier()) const override;

  /// Return the DigitizationModule
  /// @return optionally the DigitizationModule
  virtual std::shared_ptr<const DigitizationModule>
  digitizationModule() const final override;

  /// Return surface associated with this identifier,
  ///
  /// @param identifier is ignored in this case
  ///
  /// @param identifier is ignored for this simple detector element
  const Surface&
  surface(const Identifier& identifier = Identifier()) const override;

  /// Returns the full list of all detection surfaces associated
  /// to this detector element
  const std::vector<std::shared_ptr<const Surface>>&
  surfaces() const override;

  /// The maximal thickness of the detector element wrt normal axis
  double
  thickness() const override;

private:
  /// the element representation
  /// identifier
  Identifier m_elementIdentifier;
  /// the transform for positioning in 3D space
  std::shared_ptr<const Transform3D> m_elementTransform;
  /// the surface represented by it
  std::shared_ptr<const Surface> m_elementSurface;
  /// the element thickness
  double m_elementThickness;

  /// the cache
  std::vector<std::shared_ptr<const Surface>> m_elementSurfaces;
  /// store either
  std::shared_ptr<const PlanarBounds> m_elementPlanarBounds;
  std::shared_ptr<const DiscBounds>   m_elementDiscBounds;
  std::shared_ptr<const LineBounds>   m_elementLineBounds;
};

inline Identifier
DetectorElementStub::identify() const
{
  return m_elementIdentifier;
}

inline const Transform3D&
DetectorElementStub::transform(const Identifier&) const
{
  return *m_elementTransform;
}

inline std::shared_ptr<const DigitizationModule>
DetectorElementStub::digitizationModule() const
{
  return DetectorElementBase::digitizationModule();
}

inline const Surface&
DetectorElementStub::surface(const Identifier&) const
{
  return *m_elementSurface;
}

inline const std::vector<std::shared_ptr<const Surface>>&
DetectorElementStub::surfaces() const
{
  return m_elementSurfaces;
}

inline double
DetectorElementStub::thickness() const
{
  return m_elementThickness;
}

}  // end of ns