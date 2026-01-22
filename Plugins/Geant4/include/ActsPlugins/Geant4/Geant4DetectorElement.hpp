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

#include <memory>

class G4VPhysicalVolume;

namespace Acts {

class ISurfaceMaterial;
class Surface;

}  // namespace Acts

namespace ActsPlugins {
/// @addtogroup geant4_plugin
/// @{

/// @class Geant4DetectorElement
///
/// Detector element representative for Geant4 sensitive
/// elements.
class Geant4DetectorElement : public Acts::SurfacePlacementBase {
 public:
  /// Broadcast the context type
  using ContextType = Acts::GeometryContext;

  /// @brief  Constructor with arguments
  /// @param surface the surface representing this detector element
  /// @param g4physVol the physical volume representing this detector element
  /// @param toGlobal the global transformation before the volume
  /// @param thickness the thickness of this detector element
  Geant4DetectorElement(std::shared_ptr<Acts::Surface> surface,
                        const G4VPhysicalVolume& g4physVol,
                        const Acts::Transform3& toGlobal, double thickness);

  /// Return local to global transform associated with this detector element
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  const Acts::Transform3& localToGlobalTransform(
      const Acts::GeometryContext& gctx) const override;
  /// @return Reference to the local-to-global transformation matrix

  /// Return surface associated with this detector element
  const Acts::Surface& surface() const override;
  /// @return Const reference to the associated surface

  /// Non-const access to surface associated with this detector element
  Acts::Surface& surface() override;
  /// @return Mutable reference to the associated surface

  /// Return the thickness of this detector element
  /// @return The thickness value in length units
  double thickness() const override;

  /// @return to the Geant4 physical volume
  const G4VPhysicalVolume& g4PhysicalVolume() const;
  /// Is the detector element a sensitive element
  bool isSensitive() const final { return true; }

 private:
  /// Corresponding Surface
  std::shared_ptr<Acts::Surface> m_surface;
  /// The GEant4 physical volume
  const G4VPhysicalVolume* m_g4physVol{nullptr};
  /// The global transformation before the volume
  Acts::Transform3 m_toGlobal;
  ///  Thickness of this detector element
  double m_thickness{0.};
};

/// @}
}  // namespace ActsPlugins
