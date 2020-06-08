// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
/// what if not used with DD4hep?
///
class TGeoDetectorElement : public IdentifiedDetectorElement {
 public:
  /// Broadcast the context type
  using ContextType = GeometryContext;

  /// Constructor
  /// @param identifier is the detector identifier
  /// @param tGeoNode is the TGeoNode which should be represented
  /// @param tGeoMatrix The Matrix to global (i.e. ACTS transform)
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
      const Identifier& identifier, const TGeoNode& tGeoNode,
      const TGeoMatrix& tGeoMatrix = TGeoIdentity(),
      const std::string& axes = "XYZ", double scalor = 10.,
      std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr);

  ~TGeoDetectorElement() override;

  Identifier identifier() const final;

  /// Return local to global transform associated with this identifier
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  const Transform3D& transform(const GeometryContext& gctx) const final;

  /// Return surface associated with this identifier, which should come from the
  const Surface& surface() const final;

  /// Retrieve the DigitizationModule
  const std::shared_ptr<const DigitizationModule> digitizationModule()
      const final {
    return nullptr;
  };

  /// Returns the thickness of the module
  double thickness() const final;

 private:
  /// Pointer to TGeoNode (not owned)
  const TGeoNode* m_detElement{nullptr};
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
  double m_thickness{0.};
  /// Corresponding Surface
  std::shared_ptr<Surface> m_surface{nullptr};
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

}  // namespace Acts
