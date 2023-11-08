// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"

#include <stdexcept>

#include <boost/program_options.hpp>

namespace ActsExamples {
namespace DD4hep {

DD4hepDetector::DD4hepDetector() = default;

DD4hepDetector::DD4hepDetector(
    std::shared_ptr<DD4hepGeometryService> _geometryService)
    : geometryService(std::move(_geometryService)) {}

DD4hepDetector::~DD4hepDetector() = default;

auto DD4hepDetector::finalize(
    ActsExamples::DD4hep::DD4hepGeometryService::Config config,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  Acts::GeometryContext dd4HepContext;
  config.matDecorator = std::move(mdecorator);
  geometryService =
      std::make_shared<ActsExamples::DD4hep::DD4hepGeometryService>(config);
  TrackingGeometryPtr dd4tGeometry =
      geometryService->trackingGeometry(dd4HepContext);
  if (!dd4tGeometry) {
    throw std::runtime_error{
        "Did not receive tracking geometry from DD4hep geometry service"};
  }
  ContextDecorators dd4ContextDecorators = {};
  // return the pair of geometry and empty decorators
  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(
      std::move(dd4tGeometry), std::move(dd4ContextDecorators));
}

}  // namespace DD4hep
}  // namespace ActsExamples
