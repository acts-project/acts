// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Plugins/DD4hep/DD4hepFieldAdapter.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"

#include <string>
#include <cstddef>
#include <memory>
#include <stdexcept>

#include <DD4hep/Detector.h>
#include <DD4hep/Fields.h>
#include <boost/program_options.hpp>

namespace ActsExamples {
namespace DD4hep {

DD4hepDetector::DD4hepDetector(const std::vector<std::string>& _compactFiles)
    : geometryService(nullptr), compactFiles(_compactFiles) {}

DD4hepDetector::DD4hepDetector(
    std::shared_ptr<DD4hepGeometryService> _geometryService,
    const std::vector<std::string>& _compactFiles)
    : geometryService(std::move(_geometryService)),
      compactFiles(_compactFiles) {}

auto DD4hepDetector::finalize(
    ActsExamples::DD4hep::DD4hepGeometryService::Config config,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  if (geometryService == nullptr) {
    throw std::runtime_error{
        "No DD4hep geometry service configured, can not build "
        "TrackingGeometry."};
  }

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

auto DD4hepDetector::finalize(const Acts::GeometryContext& gctx)
    -> std::tuple<DetectorPtr, ContextDecorators,
                  Acts::DD4hepDetectorElement::Store> {
  
  // Check if the xml files are present
  auto dd4hepDetector = &dd4hep::Detector::getInstance();
  for (const auto& file : compactFiles) {
    dd4hepDetector->fromCompact(file.c_str());
  }
  dd4hepDetector->volumeManager();
  dd4hepDetector->apply("DD4hepVolumeManager", 0, nullptr);
  auto world = dd4hepDetector->world();

  // Build the detector structure


  // Prepare the return objects
  DetectorPtr detector = nullptr;
  ContextDecorators contextDecorators = {};
  Acts::DD4hepDetectorElement::Store detectorElements = {};

  return {detector, contextDecorators, detectorElements};
}

std::shared_ptr<Acts::DD4hepFieldAdapter> DD4hepDetector::field() const {
  const auto& detector = geometryService->detector();

  return std::make_shared<Acts::DD4hepFieldAdapter>(detector.field());
}

}  // namespace DD4hep
}  // namespace ActsExamples
