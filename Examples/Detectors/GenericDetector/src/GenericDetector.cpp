// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/GenericDetector/GenericDetector.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/GenericDetector/BuildGenericDetector.hpp"
#include "ActsExamples/GenericDetector/GenericDetectorElement.hpp"
#include "ActsExamples/GenericDetector/GenericDetectorOptions.hpp"

#include <boost/program_options.hpp>

void GenericDetector::addOptions(
    boost::program_options::options_description& opt) const {
  ActsExamples::Options::addGenericGeometryOptions(opt);
}

auto GenericDetector::finalize(
    const boost::program_options::variables_map& vm,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  Config cfg;

  cfg.buildLevel = vm["geo-generic-buildlevel"].template as<size_t>();
  // set geometry building logging level
  cfg.surfaceLogLevel =
      Acts::Logging::Level(vm["geo-surface-loglevel"].template as<size_t>());
  cfg.layerLogLevel =
      Acts::Logging::Level(vm["geo-layer-loglevel"].template as<size_t>());
  cfg.volumeLogLevel =
      Acts::Logging::Level(vm["geo-volume-loglevel"].template as<size_t>());

  cfg.buildProto = (vm["mat-input-type"].template as<std::string>() == "proto");

  return finalize(cfg, std::move(mdecorator));
}

auto GenericDetector::finalize(
    const Config& cfg,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  DetectorElement::ContextType nominalContext;
  /// Return the generic detector
  TrackingGeometryPtr gGeometry =
      ActsExamples::Generic::buildDetector<DetectorElement>(
          nominalContext, detectorStore, cfg.buildLevel, std::move(mdecorator),
          cfg.buildProto, cfg.surfaceLogLevel, cfg.layerLogLevel,
          cfg.volumeLogLevel);
  ContextDecorators gContextDeocrators = {};
  // return the pair of geometry and empty decorators
  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(
      std::move(gGeometry), std::move(gContextDeocrators));
}