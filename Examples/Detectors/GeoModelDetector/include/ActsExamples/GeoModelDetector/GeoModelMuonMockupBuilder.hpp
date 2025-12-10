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

  enum class StationIdx : uint8_t {
    BI, BM, BO,
    EAI, EAM, EAO,
    ECI, ECM, ECO,
    nStations
  };

  struct cylBounds {
    double rMin{0.0};
    double rMax{0.0};
    double zMin{0.0};
    double zMax{0.0};

    void update(const double rMinIn, const double rMaxIn, const double zMinIn, const double zMaxIn) {
      rMin = std::min(rMin, rMinIn);
      rMax = std::max(rMax, rMaxIn);
      zMin = std::min(zMin, zMinIn);
      zMax = std::max(zMax, zMaxIn);
    }
  };

  using Node_t =
      std::shared_ptr<Acts::Experimental::StaticBlueprintNode>;

  explicit GeoModelMuonMockupBuilder(
      const Config& cfg,
      std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
          "GeoModelMuonMockupBuilder", Acts::Logging::DEBUG));

  ~GeoModelMuonMockupBuilder() override = default;

  std::unique_ptr<const Acts::TrackingGeometry> trackingGeometry(
      const Acts::GeometryContext& gctx) const override;

 private:
  Config m_cfg;

  Node_t processStation(
      const std::span<ConvertedVolList_t::value_type> boundingBoxes,
      const std::string& station, const bool isBarrel,
      Acts::VolumeBoundFactory& boundFactory,
      const Acts::GeometryIdentifier& geoId) const;
  
  std::unique_ptr<Acts::TrackingVolume> buildChildChamber(
      const ConvertedVolList_t::value_type& box,
      Acts::VolumeBoundFactory& boundFactory) const;

  // Helper function returning the station idx from a box volume
  StationIdx getStationIdx(const ConvertedVolList_t::value_type& box) const;
  // Helper function converting the station idx to string
  std::string stationIdxToString(StationIdx idx) const;

  /// Private access method to the logger
  const Acts::Logger& logger() const { return *m_logger; }

  /// logging instance
  std::unique_ptr<const Acts::Logger> m_logger;
};
}  // namespace ActsExamples
