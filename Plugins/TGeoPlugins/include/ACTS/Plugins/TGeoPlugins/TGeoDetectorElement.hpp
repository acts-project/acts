// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TGeoDetectorElement.h, ACTS project, TGeoDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ACTS_TGEODETECTORELEMENT_TGEODETECTORELEMENT
#define ACTS_TGEODETECTORELEMENT_TGEODETECTORELEMENT 1

#include "ACTS/Detector/DetectorElementBase.hpp"
#include <iostream>
#include "TGeoManager.h"

namespace Acts {

/// @class TGeoDetectorElement
///
/// DetectorElement plugin for ROOT TGeo shapes. Added possibility to hand over
/// transformation matrix.
///
/// @TODO what if shape conversion failes? add implementation of more than one
/// surface per module, implementing also for other shapes->Cone,ConeSeg,Tube? what
/// if not used with DD4hep?
///
  class TGeoDetectorElement : public DetectorElementBase
{
public:
  /// Constructor  
  /// @param identifier is the detector identifier
  /// @param tGeoDetElement is the TGeoNode which should be represented
  /// @Param motherTransform is the (optional) transform applied to the TGeoNode
  TGeoDetectorElement(const Identifier&                  identifier,
                      TGeoNode*                          tGeoDetElement,
                      std::shared_ptr<Acts::Transform3D> motherTransform = nullptr,
                      double scalor = 0.);

  ///  Destructor 
  virtual ~TGeoDetectorElement();

  /// Identifier 
  virtual Identifier
  identify() const override;

  /// Return local to global transform associated with this identifier
  virtual const Transform3D&
  transform(const Identifier& identifier = Identifier()) const final;

  /// Return surface associated with this identifier, which should come from the 
  virtual const Surface&
  surface(const Identifier& identifier = Identifier()) const final;

  /// Returns the full list of all detection surfaces associated
  /// to this detector element 
  virtual const std::vector<std::shared_ptr<const Surface>>&
  surfaces() const override;

  /// Returns the thickness of the module 
  virtual double
  thickness() const override;

private:
  /// DD4hep detector element, just linked - not owned
  TGeoNode* m_detElement;
  /// Transformation of the detector element
  std::shared_ptr<const Acts::Transform3D> m_transform;
  /// Center position of the detector element
  mutable std::shared_ptr<const Vector3D> m_center;
  /// Normal vector to the detector element
  mutable std::shared_ptr<const Vector3D> m_normal;
  /// Identifier of the detector element
  const Identifier m_identifier;
  /// Boundaries of the detector element
  std::shared_ptr<const SurfaceBounds> m_bounds;
  ///  Thickness of this detector element
  double m_thickness;  //@TODO implement thickness from TGeoMode
  /// Corresponding Surface
  std::shared_ptr<const Surface> m_surface;
  /// possible contained surfaces
  std::vector<std::shared_ptr<const Surface>> m_surfaces;
};

inline Identifier
TGeoDetectorElement::identify() const
{
  return m_identifier;
}

inline const Transform3D&
TGeoDetectorElement::transform(const Identifier&) const
{
  return (*m_transform);
}

inline const Surface&
TGeoDetectorElement::surface(const Identifier&) const
{
  return (*m_surface);
}

inline const std::vector<std::shared_ptr<const Surface>>&
TGeoDetectorElement::surfaces() const
{
  return (m_surfaces);
}

inline double
TGeoDetectorElement::thickness() const
{
  return m_thickness;
}
}

#endif
