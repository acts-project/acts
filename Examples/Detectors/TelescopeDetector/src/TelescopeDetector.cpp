// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/TelescopeDetector/BuildTelescopeDetector.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetectorElement.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetectorOptions.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <boost/program_options.hpp>

void TelescopeDetector::addOptions(
    ActsExamples::Options::Description& desc) const {
  ActsExamples::Options::addTelescopeGeometryOptions(desc);
}

auto TelescopeDetector::finalize(
    const boost::program_options::variables_map& vm,
    std::shared_ptr<const Acts::IMaterialDecorator> /*mdecorator*/)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  // --------------------------------------------------------------------------------
  DetectorElement::ContextType nominalContext;
  auto positions = vm["geo-tele-positions"]
                       .template as<ActsExamples::Options::VariableReals>()
                       .values;
  auto offsets =
      vm["geo-tele-offsets"].template as<ActsExamples::Options::Reals<2>>();
  auto pSize =
      vm["geo-tele-size"].template as<ActsExamples::Options::Reals<2>>();
  // Translate the thickness in unit of mm
  auto thickness = vm["geo-tele-thickness"].template as<double>() * 0.001;
  auto binValue = vm["geo-tele-alignaxis"].template as<size_t>();
  if (binValue > 2) {
    throw std::invalid_argument("The axis value could only be 0, 1, or 2.");
  }
  // Sort the provided distances
  std::sort(positions.begin(), positions.end());

  /// Return the telescope detector
  TrackingGeometryPtr gGeometry = ActsExamples::Telescope::buildDetector(
      nominalContext, detectorStore, positions, offsets, pSize, thickness,
      static_cast<Acts::BinningValue>(binValue));
  ContextDecorators gContextDecorators = {};
  // return the pair of geometry and empty decorators
  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(
      std::move(gGeometry), std::move(gContextDecorators));
}
