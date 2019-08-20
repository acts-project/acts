// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TGeoDetectorElement.hpp, Acts project, TGeoDetector plugin
///////////////////////////////////////////////////////////////////

#pragma once
#include <iostream>
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/Plugins/Identification/Identifier.hpp"
#include "TGeoManager.h"

namespace Acts {

class ISurfaceMaterial;
class SurfaceBounds;
class DigitizationModule;

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
class TGeoDetectorElement : public IdentifiedDetectorElement {
 public:
  /// Broadcast the context type
  using ContextType = GeometryContext;

  /// Constructor
  /// @param identifier is the detector identifier
  /// @param tGeoDetElement is the TGeoNode which should be represented
  /// @param mGlobal to global is the (optional) transform applied to the
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
  /// @note This parameter only needs to be set for plane modules
  /// @param scalor is the scale factor for unit conversion if needed
  /// @param isDisc in case the sensitive detector module should be translated
  ///        as disc (e.g. for endcaps) this flag should be set to true
  /// @note In the translation from a 3D geometry (TGeo) which only knows tubes
  ///       to a 2D geometry (Tracking geometry) a distinction if the module
  ///       should be described as a cylinder or a disc surface needs to be
  ///       done. Since this information can not be taken just from the geometry
  ///       description (both can be described as TGeoTubeSeg), one needs to
  ///       set the flag 'isDisc' in case a volume with shape \c TGeoTubeSeg
  ///       should be translated to a disc surface. Per default it will be
  ///       translated into a cylindrical surface.
  /// @param material Possible material of detector element
  /// @param digitizationModule Shared pointer to the geometric digitization
  /// description
  TGeoDetectorElement(
      const Identifier& identifier, TGeoNode* tGeoDetElement,
      const TGeoMatrix* mGlobal = nullptr, const std::string& axes = "XYZ",
      double scalor = 1., bool isDisc = false,
      std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr,
      std::shared_ptr<const Acts::DigitizationModule> digitizationModule =
          nullptr);

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
  ///   - "XzY" -> negative node z axis is tracking y axis, etc.
  /// @note This parameter only needs to be set for plane modules
  /// @param scalor is the scale factor for unit conversion if needed
  /// @param isDisc in case the sensitive detector module should be translated
  ///        as disc (e.g. for endcaps) this flag should be set to true
  /// @note In the translation from a 3D geometry (TGeo) which only knows tubes
  ///       to a 2D geometry (Tracking geometry) a distinction if the module
  ///       should be described as a cylinder or a disc surface needs to be
  ///       done. Since this information can not be taken just from the geometry
  ///       description (both can be described as TGeoTubeSeg), one needs to
  ///       set the flag 'isDisc' in case a volume with shape \c TGeoTubeSeg
  ///       should be translated to a disc surface. Per default it will be
  ///       translated into a cylindrical surface.
  /// @param material Possible material of detector element
  /// @param digitizationModule Shared pointer to the geometric digitization
  /// description
  ///
  /// Note - it calles the construct() method directly
  TGeoDetectorElement(
      const Identifier& identifier, const TGeoMatrix& transform,
      TGeoNode* tGeoDetElement, const std::string& axes = "XYZ",
      double scalor = 1., bool isDisc = false,
      std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr,
      std::shared_ptr<const Acts::DigitizationModule> digitizationModule =
          nullptr);

  ///  Destructor
  ~TGeoDetectorElement() override;

  /// Identifier
  Identifier identifier() const final;

  /// Return local to global transform associated with this identifier
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  const Transform3D& transform(const GeometryContext& gctx) const final;

  /// Return surface associated with this identifier, which should come from the
  const Surface& surface() const final;

  /// Returns the thickness of the module
  double thickness() const final;

  /// Retrieve the DigitizationModule
  const std::shared_ptr<const Acts::DigitizationModule> digitizationModule()
      const final;

 protected:
  /// Construct methdo called from the constructor directly
  /// @note parameter definitions are given in the constructor
  void construct(
      const Double_t* rotation = nullptr, const Double_t* translation = nullptr,
      const std::string& axes = "XYZ", double scalor = 1., bool isDisc = false,
      std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr);

 private:
  /// DD4hep detector element, just linked - not owned
  TGeoNode* m_detElement{nullptr};
  /// Transformation of the detector element
  std::shared_ptr<const Acts::Transform3D> m_transform{nullptr};
  /// Center position of the detector element
  std::shared_ptr<const Vector3D> m_center{nullptr};
  /// Normal vector to the detector element
  std::shared_ptr<const Vector3D> m_normal{nullptr};
  /// Identifier of the detector element
  Identifier m_identifier;
  /// Boundaries of the detector element
  std::shared_ptr<const SurfaceBounds> m_bounds{nullptr};
  ///  Thickness of this detector element
  double m_thickness{0.};  //@todo implement thickness from TGeoMode
  /// Corresponding Surface
  std::shared_ptr<const Surface> m_surface{nullptr};
  /// The Digitization module
  std::shared_ptr<const Acts::DigitizationModule> m_digitizationModule{nullptr};
};

inline Identifier TGeoDetectorElement::identifier() const {
  return m_identifier;
}

inline const Transform3D& TGeoDetectorElement::transform(
    const GeometryContext& /*gctx*/) const {
  return (*m_transform);
}

inline const Surface& TGeoDetectorElement::surface() const {
  return (*m_surface);
}

inline double TGeoDetectorElement::thickness() const {
  return m_thickness;
}

inline const std::shared_ptr<const Acts::DigitizationModule>
TGeoDetectorElement::digitizationModule() const {
  return m_digitizationModule;
}
}  // namespace Acts
