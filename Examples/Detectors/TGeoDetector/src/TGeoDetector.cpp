// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/TGeoDetector/TGeoDetector.hpp"

#include "ACTFW/Framework/IContextDecorator.hpp"
#include "ACTFW/TGeoDetector/BuildTGeoDetector.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"

void TGeoDetector::addOptions(
    boost::program_options::options_description& opt) const {
  FW::Options::addTGeoGeometryOptions(opt);
}

auto TGeoDetector::finalize(const boost::program_options::variables_map& vm,
                            std::shared_ptr<const Acts::IMaterialDecorator>)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  Acts::GeometryContext tGeoContext;
  TrackingGeometryPtr tgeoTrackingGeometry =
      FW::TGeo::buildTGeoDetector(vm, tGeoContext, detectorStore);

  ContextDecorators tgeoContextDeocrators = {};
  // return the pair of geometry and empty decorators
  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(
      std::move(tgeoTrackingGeometry), std::move(tgeoContextDeocrators));
}
