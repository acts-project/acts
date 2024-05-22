// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <memory>

class GeoFullPhysVol;

namespace Acts {

class ISurfaceMaterial;
class Surface;
class PlanarBounds;

/// @class GeoModelDetectorElement
///
/// Detector element representative for GeoModel based
/// sensitive elements.
class GeoModelDetectorElement : public DetectorElementBase {
 public:
  /// Broadcast the context type
  using ContextType = GeometryContext;

  // Deleted default constructor
  GeoModelDetectorElement() = delete;

  /// @brief Factory to create a planar detector element with connected surfcace
  ///
  /// @param geoPhysVol reprsenting the physical volume
  /// @param pBounds the planar bounds
  /// @param sfTransform the surface transform
  /// @param thickness the thickness of the detector element
  /// @return
  static std::shared_ptr<GeoModelDetectorElement> createPlanarElement(
      const GeoFullPhysVol& geoPhysVol,
      const std::shared_ptr<PlanarBounds> pBounds,
      const Transform3& sfTransform, ActsScalar thickness);

  /// Constructor
  GeoModelDetectorElement(const GeoFullPhysVol& geoPhysVol,
                          const std::shared_ptr<PlanarBounds> pBounds,
                          const Transform3& sfTransform, ActsScalar thickness);

  /// Return local to global transform associated with this detector element
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  const Transform3& transform(const GeometryContext& gctx) const override;

  /// Return surface associated with this detector element
  const Surface& surface() const override;

  /// Non-const access to surface associated with this detector element
  Surface& surface() override;

  /// Return the thickness of this detector element
  ActsScalar thickness() const override;

  /// @return to the Geant4 physical volume
  const GeoFullPhysVol& physicalVolume() const;

 private:
  /// The GeoModel full physical volume
  const GeoFullPhysVol* m_geoPhysVol{nullptr};
  /// The surface
  std::shared_ptr<Surface> m_surface;
  /// The global transformation before the volume
  Transform3 m_surfaceTransform;
  ///  Thickness of this detector element
  ActsScalar m_thickness{0.};
};

}  // namespace Acts
