// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ITrackingGeometryBuilder.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorObjectFactory.hpp"

namespace ActsExamples {

/// @brief Tracking Geometry Builder implementation to connect the GeoModel detector description with the
///        translation to a tracking geometry for a mockup muon spectrometer
///        detector
///@note This is a mockup implementation of muon detector that uses the gen3 Acts geometry with a blueprint creation

class GeoModelMuonMockupBuilder : public Acts::ITrackingGeometryBuilder {
 public:
  using SensitiveSurfaces = std::vector<Acts::GeoModelSensitiveSurface>;
  using GeoModelVolumeFPVsVec =
      std::vector<Acts::GeoModelDetectorObjectFactory::GeoModelVolumeFPVTuple>;

  struct Config {
    /// The sensitive surfaces as a tuple of Detector Element and Surface
    SensitiveSurfaces sensitiveSurfaces = {};

    /// The bounding boxes as a tuple of a Volume, DetectorVolume and
    /// PVConstLInk
    GeoModelVolumeFPVsVec boundingBoxes = {};

    /// The station names to be built (e.g for barrel: BIL, BML etc)
    std::vector<std::string> stationNames = {};

    /// The number of sectors in the barrel
    std::size_t nSectors = 8;
  };

  GeoModelMuonMockupBuilder(const Config& cfg);
  ~GeoModelMuonMockupBuilder() = default;

  std::unique_ptr<const Acts::TrackingGeometry> trackingGeometry(
      const Acts::GeometryContext& gctx) const override;

 private:
  Config m_cfg;

  std::shared_ptr<Acts::Experimental::StaticBlueprintNode> buildBarrelNode(
      const GeoModelVolumeFPVsVec& boundingBoxes,
      const SensitiveSurfaces& sensitiveSurfaces,
      const std::string& name) const;
};
}  // namespace ActsExamples
