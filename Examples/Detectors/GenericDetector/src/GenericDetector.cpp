// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/GenericDetector/GenericDetector.hpp"

#include <boost/program_options.hpp>

#include "ACTFW/Framework/IContextDecorator.hpp"
#include "ACTFW/GenericDetector/BuildGenericDetector.hpp"
#include "ACTFW/GenericDetector/GenericDetectorElement.hpp"
#include "ACTFW/GenericDetector/GenericDetectorOptions.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"

void GenericDetector::addOptions(
    boost::program_options::options_description& opt) const {
  FW::Options::addGenericGeometryOptions(opt);
}

auto GenericDetector::finalize(
    const boost::program_options::variables_map& vm,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  // --------------------------------------------------------------------------------
  DetectorElement::ContextType nominalContext;

  auto buildLevel = vm["geo-generic-buildlevel"].template as<size_t>();
  // set geometry building logging level
  Acts::Logging::Level surfaceLogLevel =
      Acts::Logging::Level(vm["geo-surface-loglevel"].template as<size_t>());
  Acts::Logging::Level layerLogLevel =
      Acts::Logging::Level(vm["geo-layer-loglevel"].template as<size_t>());
  Acts::Logging::Level volumeLogLevel =
      Acts::Logging::Level(vm["geo-volume-loglevel"].template as<size_t>());

  bool buildProto =
      (vm["mat-input-type"].template as<std::string>() == "proto");

  /// Return the generic detector
  TrackingGeometryPtr gGeometry = FW::Generic::buildDetector<DetectorElement>(
      nominalContext, detectorStore, buildLevel, std::move(mdecorator),
      buildProto, surfaceLogLevel, layerLogLevel, volumeLogLevel);
  ContextDecorators gContextDeocrators = {};
  // return the pair of geometry and empty decorators
  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(
      std::move(gGeometry), std::move(gContextDeocrators));
}
