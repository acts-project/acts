// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/ITrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DetectorCommons/Detector.hpp"

#include <memory>

namespace ActsExamples {

/// @brief Detector implementation to connect the GeoModel detector description with the
///        transcription to the G4 geometry and with the tracking geometry.
struct GeoModelDetector : public Detector {
  struct Config {
    /// @brief Configured instance to the GeoModel loaded geoModel tree
    Acts::GeoModelTree geoModelTree{};
    /// @brief Path to the GeoModel file. Used if the GeoModelTree remains unconfigured
    std::string path{};
    /// @brief Logging level of the child tools
    Acts::Logging::Level logLevel = Acts::Logging::INFO;
  };
  /// @brief Standard constructor to configure the Detector building
  /// @param cfg: Options to steer the building process
  explicit GeoModelDetector(const Config& cfg);

  /// @brief Build the Geant4 detector construction
  /// @note This throws an exception if Geant4 is not enabled
  /// @param options The Geant4 construction options
  /// @return The Geant4 detector construction
  std::unique_ptr<G4VUserDetectorConstruction> buildGeant4DetectorConstruction(
      const Geant4ConstructionOptions& options) const override;

  /// @brief Instantiate the building of the tracking geometry from an external builder
  ///        The tracking geometry builder should contain the surfaces built
  ///        from the same GeoModelTree as the GeoModelDetector itself.
  /// @param gctx: Geometry context to access the aligned transform of the volumes
  /// @param builder: Reference to the TrackingGeometry builder interpreting the GeoModelTree
  /// @return A constructed tracking geometry
  std::shared_ptr<const Acts::TrackingGeometry> buildTrackingGeometry(
      const Acts::GeometryContext& gctx,
      const Acts::ITrackingGeometryBuilder& builder) {
    if (!m_trackingGeometry) {
      m_trackingGeometry = builder.trackingGeometry(gctx);
    }
    return m_trackingGeometry;
  }

 private:
  Config m_cfg{};
};

}  // namespace ActsExamples
