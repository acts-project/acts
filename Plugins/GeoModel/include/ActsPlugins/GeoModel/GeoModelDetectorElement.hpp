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

  /// @brief Factory to create a planar detector element with connected surface
  ///
  /// @tparam SurfaceType the surface type
  /// @tparam BoundsType the bounds type
  ///
  /// @param geoPhysVol representing the physical volume
  /// @param bounds the bounds class
  /// @param sfTransform the surface transform
  /// @param thickness the thickness of the detector element
  ///
  /// @return a shared pointer to an instance of the detector element
  template <typename SurfaceType, typename BoundsType>
  static std::shared_ptr<GeoModelDetectorElement> createDetectorElement(
      const PVConstLink& geoPhysVol, const std::shared_ptr<BoundsType> bounds,
      const Acts::Transform3& sfTransform, double thickness) {
    // First create the detector element with a nullptr
    auto detElement = std::make_shared<GeoModelDetectorElement>(
        geoPhysVol, nullptr, sfTransform, thickness);
    auto surface = Acts::Surface::makeShared<SurfaceType>(bounds, *detElement);
    detElement->attachSurface(surface);
    return detElement;
  }

  /// Constructor with arguments
  ///
  /// @param geoPhysVol representing the physical volume
  /// @param surface the representing surface
  /// @param sfTransform the surface transform
  /// @param thickness the thickness of the detector element
  GeoModelDetectorElement(PVConstLink geoPhysVol,
                          std::shared_ptr<Acts::Surface> surface,
                          const Acts::Transform3& sfTransform,
                          double thickness);

  /// Return local to global transform associated with this detector element
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  const Acts::Transform3& localToGlobalTransform(
      const Acts::GeometryContext& gctx) const override;

  /// Return the nominal - non-contextual transform
  const Acts::Transform3& nominalTransform() const;

  /// Return surface associated with this detector element
  const Acts::Surface& surface() const override;

  /// Non-const access to surface associated with this detector element
  Acts::Surface& surface() override;

  /// Return the thickness of this detector element
  double thickness() const;

  /// @return to the Geant4 physical volume
  PVConstLink physicalVolume() const;

  /// Get the name of the logical volume
  const std::string& logVolName() const;

  /// Get the string identifier of the corresponding database entry
  /// Note: This is not by defnitition a unique identifier, there can be
  /// several detector elements created from a single database entry.
  const std::string& databaseEntryName() const { return m_entryName; };

  /// Set the corresponding database entry string
  void setDatabaseEntryName(const std::string& n) { m_entryName = n; };
  /// Is the detector element a sensitive element
  bool isSensitive() const final { return true; }

 protected:
  /// Attach a surface
  ///
  /// @param surface The surface to attach
  void attachSurface(std::shared_ptr<Acts::Surface> surface) {
    m_surface = std::move(surface);
  }

 private:
  std::string m_entryName;
  /// The GeoModel full physical volume
  PVConstLink m_geoPhysVol{nullptr};
  /// The surface
  std::shared_ptr<Acts::Surface> m_surface;
  /// The global transformation before the volume
  Acts::Transform3 m_surfaceTransform;
  ///  Thickness of this detector element
  double m_thickness{0.};
};

/// Collect the sensitive surface & detector element
using GeoModelSensitiveSurface =
    std::tuple<std::shared_ptr<GeoModelDetectorElement>,
               std::shared_ptr<Acts::Surface>>;

/// @}

}  // namespace ActsPlugins
