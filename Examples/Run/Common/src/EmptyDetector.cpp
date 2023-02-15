// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detector/EmptyDetector.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/Options/EmptyDetectorOptions.hpp"

#include <boost/program_options.hpp>

void EmptyDetector::addOptions(
    boost::program_options::options_description& opt) const {
  ActsExamples::Options::addEmptyGeometryOptions(opt);
}

auto EmptyDetector::finalize(
    const boost::program_options::variables_map& vm,
    std::shared_ptr<const Acts::IMaterialDecorator> /*unused*/)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  using namespace Acts::UnitLiterals;

  // Create an empty cylinder with the chosen radius / halflength
  double r = vm["geo-empty-radius"].as<double>() * 1_m;
  double hz = vm["geo-empty-halfLength"].as<double>() * 1_m;

  // The cylinder volume bounds
  auto cvBounds = std::make_shared<Acts::CylinderVolumeBounds>(0.0, r, hz);

  // Create the world volume
  auto worldVolume = Acts::TrackingVolume::create(
      Acts::Transform3::Identity(), cvBounds, nullptr, "EmptyCylinder");

  // Create the tracking geometry
  auto tgGeometry =
      std::make_shared<Acts::TrackingGeometry>(std::move(worldVolume), nullptr);

  /// Empty decorators
  ContextDecorators eContextDeocrators = {};

  // And return the pair with empty decorators
  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(
      std::move(tgGeometry), std::move(eContextDeocrators));
}
