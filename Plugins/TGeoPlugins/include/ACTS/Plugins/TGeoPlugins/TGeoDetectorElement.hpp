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

#include <iostream>
#include "ACTS/Detector/DetectorElementBase.hpp"
#include "TGeoManager.h"

namespace Acts {

class DigitizationModule;
class SurfaceMaterial;

/// @class TGeoDetectorElement
///
/// DetectorElement plugin for ROOT TGeo shapes. Added possibility to hand over
/// transformation matrix.
///
/// @todo what if shape conversion failes? add implementation of more than one
/// surface per module, implementing also for other shapes->Cone,ConeSeg,Tube?
/// what
/// if not used with DD4hep?
///
class TGeoDetectorElement : public DetectorElementBase
{
public:
  /// Constructor
  /// @param identifier is the detector identifier
  /// @param tGeoDetElement is the TGeoNode which should be represented
  /// @param mtoglobal to global is the (optional) transform applied to the
  /// TGeoNode
  /// @param axes is the axis orientation with respect to the tracking frame
  ///        it is a string of the three characters x, y and z (standing for the
  ///        three axes)
  ///        there is a distinction between capital and lower case characters :
  ///        - capital      -> positive orientation of the axis
  ///        - lower case   -> negative orientation of the axis
  ///        example options are "XYZ" -> identical frame definition (default
  ///        value)
  ///                            "YZX" -> node y axis is tracking x axis, etc.
  ///                            "XzY" -> negative node z axis is tracking y
  ///                            axis, etc.
  /// @param scalor is the scale factor for unit conversion if needed
  /// @param material Possible material of detector element
  TGeoDetectorElement(const Identifier&  identifier,
                      TGeoNode*          tGeoDetElement,
                      const TGeoMatrix*  mtoglobal = nullptr,
                      const std::string& axes      = "XYZ",
                      double             scalor    = 1.,
                      std::shared_ptr<const Acts::SurfaceMaterial> material
                      = nullptr);

  /// Alternative Constructor
  /// when the localToGlobal transform is already known for the detector element
  /// (e.g. usage in DD4hepPlugin)
  /// @param identifier is the detector identifier
  /// @param transform the full transformation matrix of the detector element
  /// (local to global)
  /// @param tGeoDetElement is the TGeoNode which should be represented
  /// @param axes is the axis orientation with respect to the tracking frame
  ///        it is a string of the three characters x, y and z (standing for the
  ///        three axes)
  ///        There is a distinction between
  /// capital and lower case
  /// characters :
  /// 	- capital      -> positive orientation of the axis
  ///		- lower case   -> negative oriantation of the axis
  ///
  ///
  /// Example options are:
  /// 	- "XYZ" -> identical frame definition (default value)
  /// 	- "YZX" -> node y axis is tracking x axis, etc.
  ///		- "XzY" -> negative node z axis is tracking y axis, etc.
  /// @param scalor is the scale factor for unit conversion if needed
  /// @param material Possible material of detector element
  TGeoDetectorElement(const Identifier&  identifier,
                      const TGeoMatrix&  transform,
                      TGeoNode*          tGeoDetElement,
                      const std::string& axes   = "XYZ",
                      double             scalor = 1.,
                      std::shared_ptr<const Acts::SurfaceMaterial> material
                      = nullptr);

  ///  Destructor
  virtual ~TGeoDetectorElement();

  /// Identifier
  virtual Identifier
  identify() const final override;

  /// Return local to global transform associated with this identifier
  virtual const Transform3D&
  transform(const Identifier& identifier = Identifier()) const final override;

  /// Set the identifier after construction (sometimes needed)
  virtual void
  assignIdentifier(const Identifier& identifier) final override;

  /// Return surface associated with this identifier, which should come from the
  virtual const Surface&
  surface(const Identifier& identifier = Identifier()) const final override;

  /// Returns the full list of all detection surfaces associated
  /// to this detector element
  virtual const std::vector<std::shared_ptr<const Surface>>&
  surfaces() const final override;

  /// Return the DigitizationModule
  /// @return optionally the DigitizationModule
  virtual std::shared_ptr<const DigitizationModule>
  digitizationModule() const override;

  /// Returns the thickness of the module
  virtual double
  thickness() const final override;

private:
  /// DD4hep detector element, just linked - not owned
  TGeoNode* m_detElement;
  /// Transformation of the detector element
  std::shared_ptr<const Acts::Transform3D> m_transform;
  /// Center position of the detector element
  std::shared_ptr<const Vector3D> m_center;
  /// Normal vector to the detector element
  std::shared_ptr<const Vector3D> m_normal;
  /// Identifier of the detector element
  Identifier m_identifier;
  /// Boundaries of the detector element
  std::shared_ptr<const SurfaceBounds> m_bounds;
  ///  Thickness of this detector element
  double m_thickness;  //@todo implement thickness from TGeoMode
  /// Corresponding Surface
  std::shared_ptr<Surface> m_surface;
  /// possible contained surfaces
  std::vector<std::shared_ptr<const Surface>> m_surfaces;
};

inline Identifier
TGeoDetectorElement::identify() const
{
  return m_identifier;
}

inline void
TGeoDetectorElement::assignIdentifier(const Identifier& identifier)
{
  m_identifier = identifier;
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
