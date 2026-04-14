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

#include <GeoModelKernel/GeoFullPhysVol.h>

class GeoVPhysVol;

namespace Acts {

class CylinderBounds;
class LineBounds;
class PlanarBounds;
}  // namespace Acts

namespace ActsPlugins {

/// @addtogroup geomodel_plugin
/// @{

/// @class GeoModelDetectorElement
///
/// Detector element representative for GeoModel based
/// sensitive elements.
class GeoModelDetectorElement : public Acts::SurfacePlacementBase {
 public:
  /// Broadcast the context type
  using ContextType = Acts::GeometryContext;

  // Deleted default constructor
  GeoModelDetectorElement() = delete;

  /// @brief Factory to create a detector element wired to an existing surface
  ///
  /// @param geoPhysVol representing the physical volume
  /// @param sfTransform the surface transform
  /// @param thickness the thickness of the detector element
  /// @param surface the surface to associate with this element
  ///
  /// @return a shared pointer to an instance of the detector element
  static std::shared_ptr<GeoModelDetectorElement> createDetectorElement(
      const PVConstLink& geoPhysVol, const Acts::Transform3& sfTransform,
      double thickness, const std::shared_ptr<Acts::Surface>& surface) {
    auto detElement = std::make_shared<GeoModelDetectorElement>(
        geoPhysVol, sfTransform, thickness);
    detElement->assignSurface(surface);
    return detElement;
  }

  /// Constructor with arguments
  ///
  /// @param geoPhysVol representing the physical volume
  /// @param surface the representing surface
  /// @param sfTransform the surface transform
  /// @param thickness the thickness of the detector element
  GeoModelDetectorElement(PVConstLink geoPhysVol,
                          const Acts::Transform3& sfTransform,
                          double thickness);

  /// Return local to global transform associated with this detector element
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @return The local to global transform
  const Acts::Transform3& localToGlobalTransform(
      const Acts::GeometryContext& gctx) const override;

  /// Return the nominal - non-contextual transform
  /// @return The nominal transform
  const Acts::Transform3& nominalTransform() const;

  /// Return the thickness of this detector element
  /// @return The thickness
  double thickness() const;

  /// @return to the Geant4 physical volume
  PVConstLink physicalVolume() const;

  /// Get the name of the logical volume
  /// @return The logical volume name
  const std::string& logVolName() const;

  /// Get the string identifier of the corresponding database entry
  /// Note: This is not by defnitition a unique identifier, there can be
  /// several detector elements created from a single database entry.
  /// @return The database entry name
  const std::string& databaseEntryName() const { return m_entryName; };

  /// Set the corresponding database entry string
  /// @param n The database entry name
  void setDatabaseEntryName(const std::string& n) { m_entryName = n; };
  /// Is the detector element a sensitive element
  /// @return True if sensitive
  bool isSensitive() const final { return true; }

 protected:
 private:
  std::string m_entryName;
  /// The GeoModel full physical volume
  PVConstLink m_geoPhysVol{nullptr};
  /// The global transformation before the volume
  Acts::Transform3 m_surfaceTransform{Acts::Transform3::Identity()};
  ///  Thickness of this detector element
  double m_thickness{0.};
};

/// Collect the sensitive surface & detector element
using GeoModelSensitiveSurface =
    std::tuple<std::shared_ptr<GeoModelDetectorElement>,
               std::shared_ptr<Acts::Surface>>;

/// @}

}  // namespace ActsPlugins
