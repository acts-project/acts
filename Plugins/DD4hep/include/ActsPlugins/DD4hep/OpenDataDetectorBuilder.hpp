// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

namespace dd4hep {
class Detector;
}

namespace Acts {
class GeometryContext;
class Logger;
class TrackingGeometry;
}  // namespace Acts

namespace ActsPlugins::DD4hep {

/// Build the Open Data Detector tracking geometry using the BarrelEndcap
/// construction path (BarrelEndcapAssembler wrapping ElementLayerAssembler).
std::unique_ptr<Acts::TrackingGeometry> buildOpenDataDetectorBarrelEndcap(
    const dd4hep::Detector& detector, const Acts::GeometryContext& gctx,
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
