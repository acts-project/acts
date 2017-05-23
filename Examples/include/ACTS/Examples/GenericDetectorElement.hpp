// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// GenericDetectorElement.h, ACTS project, Generic Detector plugin
///////////////////////////////////////////////////////////////////

#ifndef AGD_GENERICDETECTORELEMENT_GENERICDETECTORELEMENT
#define AGD_GENERICDETECTORELEMENT_GENERICDETECTORELEMENT 1

#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Identifier.hpp"
#include "ACTS/Detector/DetectorElementBase.hpp"

namespace Acts {

class Surface;
class PlanarBounds;
class DiscBounds;
class SurfaceMaterial;
class DigitizationModule;

/// @class GenericDetectorElement
///
/// This is a lightweight type of detector element,
/// it simply implements the base class.
///
class GenericDetectorElement : public DetectorElementBase
{
public:
  /// Constructor for single sided detector element
  /// - bound to a Plane Surface
  ///
  /// @param identifier is the module identifier
  /// @param transform is the transform that element the layer in 3D frame
  /// @param pBounds is the planar bounds for the planar detector element
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  GenericDetectorElement(const Identifier                       identifier,
                         std::shared_ptr<const Transform3D>     transform,
                         std::shared_ptr<const PlanarBounds>    pBounds,
                         double                                 thickness,
                         std::shared_ptr<const SurfaceMaterial> material
                         = nullptr,
                         std::shared_ptr<const DigitizationModule> dModule
                         = nullptr);

  /// Constructor for single sided detector element
  /// - bound to a Disc Surface
  ///
  /// @param identifier is the module identifier
  /// @param transform is the transform that element the layer in 3D frame
  /// @param dBounds is the planar bounds for the disc like detector element
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  GenericDetectorElement(const Identifier                       identifier,
                         std::shared_ptr<const Transform3D>     transform,
                         std::shared_ptr<const DiscBounds>      dBounds,
                         double                                 thickness,
                         std::shared_ptr<const SurfaceMaterial> material
                         = nullptr);

  ///  Destructor
  ~GenericDetectorElement();

  /// Identifier
  Identifier
  identify() const final override;

  /// Return local to global transform associated with this identifier
  ///
  /// @note this is called from the surface().transform() in the PROXY mode
  ///
  /// @param identifier is ignored for this simple detector element
  const Transform3D&
  transform(const Identifier& identifier = Identifier()) const final override;

  /// Return surface associated with this identifier,
  ///
  /// @param identifier is ignored in this case
  ///
  /// @param identifier is ignored for this simple detector element
  const Surface&
  surface(const Identifier& identifier = Identifier()) const final override;

  /// Returns the full list of all detection surfaces associated
  /// to this detector element
  const std::vector<std::shared_ptr<const Surface>>&
  surfaces() const final override;

  /// Return the DigitizationModule
  /// @return optionally the DigitizationModule
  std::shared_ptr<const DigitizationModule>
  digitizationModule() const final override;

  /// Set the identifier after construction (sometimes needed)
  virtual void
  assignIdentifier(const Identifier& identifier) final override;

  /// The maximal thickness of the detector element wrt normal axis
  double
  thickness() const final override;

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

  // the digitization module
  std::shared_ptr<const DigitizationModule> m_digitizationModule;
};

inline std::shared_ptr<const DigitizationModule>
GenericDetectorElement::digitizationModule() const
{
  return m_digitizationModule;
}

inline void
GenericDetectorElement::assignIdentifier(const Identifier& identifier)
{
  m_elementIdentifier = identifier;
}

inline Identifier
GenericDetectorElement::identify() const
{
  return m_elementIdentifier;
}

inline const Transform3D&
GenericDetectorElement::transform(const Identifier&) const
{
  return *m_elementTransform;
}

inline const Surface&
GenericDetectorElement::surface(const Identifier&) const
{
  return *m_elementSurface;
}

inline const std::vector<std::shared_ptr<const Surface>>&
GenericDetectorElement::surfaces() const
{
  return m_elementSurfaces;
}

inline double
GenericDetectorElement::thickness() const
{
  return m_elementThickness;
}

}  // end of ns

#endif
