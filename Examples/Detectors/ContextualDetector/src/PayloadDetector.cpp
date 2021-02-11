// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/ContextualDetector/PayloadDetector.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/ContextualDetector/PayloadDecorator.hpp"
#include "ActsExamples/ContextualDetector/PayloadDetectorElement.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/GenericDetector/BuildGenericDetector.hpp"
#include "ActsExamples/GenericDetector/GenericDetectorOptions.hpp"

#include <boost/program_options.hpp>

void PayloadDetector::addOptions(
    boost::program_options::options_description& opt) const {
  /// Add the generic geometry options
  ActsExamples::Options::addGenericGeometryOptions(opt);
  // specify the rotation setp
  opt.add_options()(
      "align-rotation-step",
      boost::program_options::value<double>()->default_value(0.25 * M_PI),
      "Rotation step of the RotationDecorator")(
      "align-loglevel",
      boost::program_options::value<size_t>()->default_value(3),
      "Output log level of the alignment decorator.");
}

auto PayloadDetector::finalize(
    const boost::program_options::variables_map& vm,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  // --------------------------------------------------------------------------------
  DetectorElement::ContextType nominalContext;
  // set geometry building logging level
  Acts::Logging::Level surfaceLogLevel =
      Acts::Logging::Level(vm["geo-surface-loglevel"].template as<size_t>());
  Acts::Logging::Level layerLogLevel =
      Acts::Logging::Level(vm["geo-layer-loglevel"].template as<size_t>());
  Acts::Logging::Level volumeLogLevel =
      Acts::Logging::Level(vm["geo-volume-loglevel"].template as<size_t>());

  bool buildProto =
      (vm["mat-input-type"].template as<std::string>() == "proto");

  /// return the generic detector - with payload context decorator
  TrackingGeometryPtr pTrackingGeometry =
      ActsExamples::Generic::buildDetector<DetectorElement>(
          nominalContext, detectorStore, 0, std::move(mdecorator), buildProto,
          surfaceLogLevel, layerLogLevel, volumeLogLevel);

  ContextDecorators pContextDecorators = {};

  // Alignment service
  Decorator::Config agcsConfig;
  agcsConfig.trackingGeometry = pTrackingGeometry;
  agcsConfig.rotationStep = vm["align-rotation-step"].template as<double>();

  Acts::Logging::Level decoratorLogLevel =
      Acts::Logging::Level(vm["align-loglevel"].template as<size_t>());

  // Create the service
  auto agcDecorator = std::make_shared<Decorator>(
      agcsConfig,
      Acts::getDefaultLogger("PayloadDecorator", decoratorLogLevel));
  pContextDecorators.push_back(agcDecorator);

  // return the pair of geometry and the alignment decorator(s)
  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(
      std::move(pTrackingGeometry), std::move(pContextDecorators));
}
