// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4Detector/Geant4Detector.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"

#include <boost/program_options.hpp>

void Geant4Detector::addOptions(
    boost::program_options::options_description& /*unused*/) const {
}

auto Geant4Detector::finalize(
    const boost::program_options::variables_map& /*unused*/,
    std::shared_ptr<const Acts::IMaterialDecorator> /*unused*/)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(nullptr, {});
}

auto Geant4Detector::construct()
    -> std::pair<DetectorPtr, ContextDecorators> {
  return std::make_pair<DetectorPtr, ContextDecorators>(nullptr, {});
}