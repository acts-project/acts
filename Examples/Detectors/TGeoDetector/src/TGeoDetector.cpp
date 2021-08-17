// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TGeoDetector/TGeoDetector.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/TGeoDetector/BuildTGeoDetector.hpp"

#include <boost/program_options.hpp>

namespace ActsExamples {

void TGeoDetector::addOptions(
    boost::program_options::options_description& opt) const {
  ActsExamples::Options::addTGeoGeometryOptions(opt);
}

auto TGeoDetector::finalize(
    const boost::program_options::variables_map& vm,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  Config config;
  config.surfaceLogLevel =
      Acts::Logging::Level(vm["geo-surface-loglevel"].template as<size_t>());
  config.layerLogLevel =
      Acts::Logging::Level(vm["geo-layer-loglevel"].template as<size_t>());
  config.volumeLogLevel =
      Acts::Logging::Level(vm["geo-volume-loglevel"].template as<size_t>());

  config.fileName = vm["geo-tgeo-filename"].template as<std::string>();

  config.buildBeamPipe = vm.count("geo-tgeo-beampipe-parameters") > 0;
  if (config.buildBeamPipe) {
    auto beamPipeParameters =
        vm["geo-tgeo-beampipe-parameters"].template as<Options::Reals<3>>();
    config.beamPipeRadius = beamPipeParameters[0];
    config.beamPipeHalflengthZ = beamPipeParameters[1];
    config.beamPipeLayerThickness = beamPipeParameters[2];
  }

  return finalize(config, std::move(mdecorator));
}

auto TGeoDetector::finalize(
    const Config& cfg,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  Acts::GeometryContext tGeoContext;
  TrackingGeometryPtr tgeoTrackingGeometry =
      ActsExamples::TGeo::buildTGeoDetector(cfg, tGeoContext, detectorStore,
                                            mdecorator);

  ContextDecorators tgeoContextDeocrators = {};
  // Return the pair of geometry and empty decorators
  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(
      std::move(tgeoTrackingGeometry), std::move(tgeoContextDeocrators));
}

}  // namespace ActsExamples