// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/MultiIndex.hpp"
#include "ActsPlugins/GeoModel/GeoModelDetectorElement.hpp"

namespace ActsPlugins {

/// @addtogroup geomodel_plugin
/// @{

/// Identifier class for ITk detector elements
class ITkIdentifier {
  Acts::MultiIndex<std::size_t, 1, 2, 20, 1, 19, 20, 1> m_identifier{};

 public:
  /// Constructor
  /// @param hardware Hardware type (pixel=0, strip=1)
  /// @param barrelEndcap Barrel/endcap specifier (-2,0,2)
  /// @param layerWheel Layer or wheel number
  /// @param etaModule Eta module index
  /// @param phiModule Phi module index
  /// @param side Side for double sided strip modules
  ITkIdentifier(int hardware, int barrelEndcap, int layerWheel, int etaModule,
                int phiModule, int side);

  /// Access the hardware specifier (pixel=0, strip=1)
  /// @return Hardware type
  int hardware() const;

  /// Access the barrel-endcap specifier (-2,0,2)
  /// @return Barrel/endcap identifier
  int barrelEndcap() const;

  /// Access the layer specifier
  /// @return Layer or wheel number
  int layerWheel() const;

  /// Access the phi module specifier
  /// @return Phi module index
  int phiModule() const;

  /// Access the eta module specifier
  /// @return Eta module index
  int etaModule() const;

  /// Access the side (for double sided strip modules)
  /// @return Side identifier
  int side() const;

  /// A unique identifier that represents the combination of specifiers
  /// @return Unique identifier value
  std::size_t value() const;
};

/// Output stream operator for ITkIdentifier
/// @param os The output stream
/// @param id The identifier to output
/// @return The output stream
std::ostream& operator<<(std::ostream& os, const ITkIdentifier& id);

/// Specialization of the GeoModelDetectorElement for the ITk. This allows
/// mapping of Acts::GeometryIdentifiers to ITk modules in a straight-forward
/// way.
class GeoModelDetectorElementITk : public GeoModelDetectorElement {
 public:
  /// Constructor
  /// @param geoPhysVol GeoModel physical volume
  /// @param surface Surface object
  /// @param sfTransform Surface transform
  /// @param thickness Detector element thickness
  /// @param hardware Hardware identifier
  /// @param barrelEndcap Barrel/endcap identifier
  /// @param layerWheel Layer/wheel identifier
  /// @param etaModule Eta module identifier
  /// @param phiModule Phi module identifier
  /// @param side Side identifier
  GeoModelDetectorElementITk(const PVConstLink& geoPhysVol,
                             std::shared_ptr<Acts::Surface> surface,
                             const Acts::Transform3& sfTransform,
                             double thickness, int hardware, int barrelEndcap,
                             int layerWheel, int etaModule, int phiModule,
                             int side)
      : GeoModelDetectorElement(geoPhysVol, std::move(surface), sfTransform,
                                thickness),
        m_identifier(hardware, barrelEndcap, layerWheel, etaModule, phiModule,
                     side) {}

  /// Get the ITk identifier
  /// @return The ITk identifier
  ITkIdentifier identifier() const { return m_identifier; }

  /// Convert a GeoModelDetectorElement to a GeoModelDetectorElementITk
  /// A new surface is constructed.
  /// @param detEl Detector element to convert
  /// @param srf Surface object
  /// @param gctx Geometry context
  /// @param hardware Hardware identifier
  /// @param barrelEndcap Barrel/endcap identifier
  /// @param layerWheel Layer/wheel identifier
  /// @param etaModule Eta module identifier
  /// @param phiModule Phi module identifier
  /// @param side Side identifier
  /// @return Converted detector element and surface
  /// @todo Remove redundancy in signature once plugin is refactored
  static std::tuple<std::shared_ptr<GeoModelDetectorElementITk>,
                    std::shared_ptr<Acts::Surface>>
  convertFromGeomodel(std::shared_ptr<GeoModelDetectorElement> detEl,
                      std::shared_ptr<Acts::Surface> srf,
                      const Acts::GeometryContext& gctx, int hardware,
                      int barrelEndcap, int layerWheel, int etaModule,
                      int phiModule, int side);

 private:
  ITkIdentifier m_identifier;
};

/// @}

}  // namespace ActsPlugins
