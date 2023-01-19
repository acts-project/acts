// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"

class G4VPhysicalVolume;

namespace Acts {

class ISurfaceMaterial;

/// @class Geant4DetectorElement
///
/// Detector element representative for Geant4 sensitive
/// elements.
class Geant4DetectorElement : public DetectorElementBase {
 public:
  /// Broadcast the context type
  using ContextType = GeometryContext;

  /// @brief  Constructor with arguments
  /// @param surface the surface representing this detector element
  /// @param thickness the thickness of this detector element
  /// @param g4phys the physical volume representing this detector element
  Geant4DetectorElement(std::shared_ptr<Surface> surface, ActsScalar thickness,
                        const G4VPhysicalVolume& g4phys);

  /// Return local to global transform associated with this identifier
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  const Transform3& transform(const GeometryContext& gctx) const override;

  /// Return surface associated with this identifier, which should come from the
  const Surface& surface() const override;

  /// Return the thickness of this detector element
  ActsScalar thickness() const override;

  /// @return to the Geant4 physical volume
  const G4VPhysicalVolume& g4PhysicalVolume() const;

 private:
  /// Corresponding Surface
  std::shared_ptr<Surface> m_surface{nullptr};
  ///  Thickness of this detector element
  ActsScalar m_thickness{0.};
  /// The GEant4 physical volume
  const G4VPhysicalVolume* m_g4physVol{nullptr};
};

}  // namespace Acts
