// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/DD4hepDetector/DD4hepDetector.hpp"

#include <boost/program_options.hpp>

#include "ACTFW/DD4hepDetector/DD4hepDetectorOptions.hpp"
#include "ACTFW/DD4hepDetector/DD4hepGeometryService.hpp"
#include "ACTFW/Framework/IContextDecorator.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

void DD4hepDetector::addOptions(
    boost::program_options::options_description& opt) const {
  FW::Options::addDD4hepOptions(opt);
}

auto DD4hepDetector::finalize(
    const boost::program_options::variables_map& vm,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  Acts::GeometryContext dd4HepContext;
  // read the detector config & dd4hep detector
  auto dd4HepDetectorConfig =
      FW::Options::readDD4hepConfig<po::variables_map>(vm);
  dd4HepDetectorConfig.matDecorator = mdecorator;
  auto geometrySvc =
      std::make_shared<FW::DD4hep::DD4hepGeometryService>(dd4HepDetectorConfig);
  TrackingGeometryPtr dd4tGeometry =
      geometrySvc->trackingGeometry(dd4HepContext);
  ContextDecorators dd4ContextDeocrators = {};
  // return the pair of geometry and empty decorators
  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(
      std::move(dd4tGeometry), std::move(dd4ContextDeocrators));
}
