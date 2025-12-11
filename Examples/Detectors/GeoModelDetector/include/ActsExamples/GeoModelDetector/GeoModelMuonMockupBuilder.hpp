// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ITrackingGeometryBuilder.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TransformRange.hpp"
#include "ActsPlugins/GeoModel/GeoModelDetectorObjectFactory.hpp"

namespace ActsExamples {

/// @brief Tracking Geometry Builder implementation to connect the GeoModel detector description with the
///        translation to a gen-3 tracking geometry for a mockup muon
///        spectrometer.

class GeoModelMuonMockupBuilder : public Acts::ITrackingGeometryBuilder {
 public:
  /** @brief Recycle the tuple of Volume, DetectorVolume, PVConstLink */
  using SensitiveSurfaces = std::vector<ActsPlugins::GeoModelSensitiveSurface>;
  using ConvertedVolList_t =
      ActsPlugins::GeoModelDetectorObjectFactory::ConvertedVolList_t;

  struct Config {
    /// The converted GeoModel volume objects
    ConvertedVolList_t volumeBoxFPVs{};

    /// The station names to be built (e.g for barrel: BIL, BML etc)
    std::vector<std::string> stationNames{};

    /// Pointer to the volume bound factory to share the bounds across several
    /// volumes
    std::shared_ptr<Acts::VolumeBoundFactory> volumeBoundFactory =
        std::make_shared<Acts::VolumeBoundFactory>();
  };

  explicit GeoModelMuonMockupBuilder(
      const Config& cfg,
      std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
          "GeoModelMuonMockupBuilder", Acts::Logging::DEBUG));

  ~GeoModelMuonMockupBuilder() override = default;

  std::unique_ptr<const Acts::TrackingGeometry> trackingGeometry(
      const Acts::GeometryContext& gctx) const override;

 private:
  Config m_cfg;

  std::shared_ptr<Acts::Experimental::StaticBlueprintNode> buildBarrelNode(
      const ConvertedVolList_t& boundingBoxes, const std::string& name,
      Acts::VolumeBoundFactory& boundFactory,
      const Acts::GeometryIdentifier& geoId) const;

  /// Private access method to the logger
  const Acts::Logger& logger() const { return *m_logger; }

  /// logging instance
  std::unique_ptr<const Acts::Logger> m_logger;
};
}  // namespace ActsExamples
