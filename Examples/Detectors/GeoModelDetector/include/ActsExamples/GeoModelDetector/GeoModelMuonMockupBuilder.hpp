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
  // Enum class for the station indices
  enum class StationIdx : std::uint8_t {
    BI,   // Inner Barrel
    BM,   // Middle Barrel
    BO,   // Outer Barrel
    EAI,  // Endcap Inner A-side
    EAM,  // Endcap Middle A-side
    EAO,  // Endcap Outer A-side
    ECI,  // Endcap Inner C-side
    ECM,  // Endcap Middle C-side
    ECO,  // Endcap Outer C-side
    nStations
  };
  // Enum class for the first level of container indices
  enum class FirstContainerIdx : std::uint8_t {
    Central,  // Central container for barrel & NSWs
    BW_A,     // Big Wheel A-side
    BW_C,     // Big Wheel C-side
    nFirstContainers
  };
  // Enum class for the second level of container indices when the first is
  // Central
  enum class SecondContainerIdx : std::uint8_t {
    Barrel,  // Barrel container
    NSWs,    // NSW container
    nSecondContainers
  };

  Config m_cfg;

  using Box_t = ConvertedVolList_t::value_type;
  using Node_t = Acts::Experimental::StaticBlueprintNode;
  using NodePtr_t = std::shared_ptr<Node_t>;

  /// @brief Produce a station node from the provided converted volume boxes
  NodePtr_t processStation(const std::span<Box_t> boundingBoxes,
                           const std::string& station, const bool isBarrel,
                           Acts::VolumeBoundFactory& boundFactory,
                           const Acts::GeometryIdentifier& geoId) const;

  /// @brief Build a child chamber volume from the provided converted volume box
  std::unique_ptr<Acts::TrackingVolume> buildChildChamber(
      const Box_t& box, Acts::VolumeBoundFactory& boundFactory) const;

  /// @brief Helper struct to store cylinder bounds, used to compute the overall bounds
  ///        of a station tracking volume from its component volumes
  struct cylBounds {
    /// @brief Lowest radius
    double rMin{std::numeric_limits<double>::max()};
    /// @brief Highest radius
    double rMax{std::numeric_limits<double>::lowest()};
    /// @brief Lowest longitudinal coordinate
    double zMin{std::numeric_limits<double>::max()};
    /// @brief Highest longitudinal coordinate
    double zMax{std::numeric_limits<double>::lowest()};
  };
  /// @brief Helper function to update the cylinder bounds of each station across its component volumes (chambers)
  /// @tparam VolBounds_t: The bounds type of the chamber volume
  /// @param volume: The chamber volume to extract the bounds from
  /// @param bounds: The cylBounds object to be updated
  template <Acts::VolumeBounds::BoundsType VolBounds_t>
  void updateBounds(const Acts::TrackingVolume& volume,
                    cylBounds& bounds) const;

  // Helper function returning the station idx from a box volume
  StationIdx getStationIdx(const Box_t& box) const;
  // Helper function returning the first-level container idx from a station idx
  FirstContainerIdx getFirstContainerIdx(const StationIdx& stationIdx) const;
  // Helper function converting the station idx to string
  static std::string stationIdxToString(const StationIdx idx);
  // Helper function converting the first-level container idx to string
  static std::string firstContainerIdxToString(const FirstContainerIdx idx);

  /// Private access method to the logger
  const Acts::Logger& logger() const { return *m_logger; }

  /// logging instance
  std::unique_ptr<const Acts::Logger> m_logger;
};
}  // namespace ActsExamples
