// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

#include <memory>
#include <regex>
#include <string>
#include <string_view>

namespace dd4hep {
class Detector;
}

class TGeoNode;

namespace Acts {
class GeometryContext;
class Logger;
class TrackingGeometry;
}  // namespace Acts

namespace ActsPlugins::DD4hep {

namespace detail {

inline const std::regex kPixelLayerFilter{
    "(?:PixelLayer|PixelEndcap[NP])(\\d)"};
inline const std::regex kShortStripLayerFilter{
    "(?:ShortStripLayer|ShortStripEndcap[NP])(\\d)"};
inline const std::regex kLongStripLayerFilter{
    "(?:LongStripLayer|LongStripEndcap[NP])(\\d)"};
inline const std::regex kPixelBarrelLayerFilter{"PixelLayer\\d"};
inline const std::regex kPixelNegativeEndcapLayerFilter{"PixelEndcapN\\d"};
inline const std::regex kPixelPositiveEndcapLayerFilter{"PixelEndcapP\\d"};
inline const std::regex kShortStripBarrelLayerFilter{"ShortStripLayer\\d"};
inline const std::regex kShortStripNegativeEndcapLayerFilter{
    "ShortStripEndcapN\\d"};
inline const std::regex kShortStripPositiveEndcapLayerFilter{
    "ShortStripEndcapP\\d"};
inline const std::regex kLongStripBarrelLayerFilter{"LongStripLayer\\d"};
inline const std::regex kLongStripNegativeEndcapLayerFilter{
    "LongStripEndcapN\\d"};
inline const std::regex kLongStripPositiveEndcapLayerFilter{
    "LongStripEndcapP\\d"};

inline const Acts::ExtentEnvelope kBlueprintEnvelope =
    Acts::ExtentEnvelope::Zero()
        .set(Acts::AxisDirection::AxisZ, {20., 20.})
        .set(Acts::AxisDirection::AxisR, {0., 20.});

inline const Acts::ExtentEnvelope kLayerEnvelope =
    Acts::ExtentEnvelope::Zero()
        .set(Acts::AxisDirection::AxisZ, {2., 2.})
        .set(Acts::AxisDirection::AxisR, {2., 2.});

inline int layerIndexFromName(std::string_view elemName,
                              const std::regex& layerFilter) {
  std::cmatch match;
  if (std::regex_search(elemName.begin(), elemName.end(), match, layerFilter) &&
      match.size() > 1) {
    return std::stoi(match[1].str());
  }

  if (std::regex groupedLayerNameFilter{"layer(\\d+)"};
      std::regex_search(elemName.begin(), elemName.end(), match,
                        groupedLayerNameFilter) &&
      match.size() > 1) {
    return std::stoi(match[1].str());
  }

  return 0;
}

}  // namespace detail

/// Build the Open Data Detector tracking geometry using the BarrelEndcap
/// construction path (BarrelEndcapAssembler wrapping ElementLayerAssembler).
std::unique_ptr<Acts::TrackingGeometry> buildOpenDataDetectorBarrelEndcap(
    const dd4hep::Detector& detector, const Acts::GeometryContext& gctx,
    const Acts::Logger& logger);

/// Build the Open Data Detector tracking geometry using the TGeo backend with
/// metadata extracted from DD4hep and explicit ODD layer-name patterns.
std::unique_ptr<Acts::TrackingGeometry>
buildOpenDataDetectorBarrelEndcapViaTGeo(const TGeoNode& rootNode,
                                         const Acts::GeometryContext& gctx,
                                         const Acts::Logger& logger);

/// Build the Open Data Detector tracking geometry using the DirectLayer
/// construction path (ElementLayerAssembler directly).
std::unique_ptr<Acts::TrackingGeometry> buildOpenDataDetectorDirectLayer(
    const dd4hep::Detector& detector, const Acts::GeometryContext& gctx,
    const Acts::Logger& logger);

/// Build the Open Data Detector tracking geometry using the DirectLayerGrouped
/// construction path (SensorLayerAssembler with groupBy).
std::unique_ptr<Acts::TrackingGeometry> buildOpenDataDetectorDirectLayerGrouped(
    const dd4hep::Detector& detector, const Acts::GeometryContext& gctx,
    const Acts::Logger& logger);

}  // namespace ActsPlugins::DD4hep
