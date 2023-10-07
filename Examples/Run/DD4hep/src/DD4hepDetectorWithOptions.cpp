// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detector/DD4hepDetectorWithOptions.hpp"

#include "ActsExamples/Detector/DD4hepDetectorOptions.hpp"

#include <boost/program_options.hpp>

namespace ActsExamples {

void DD4hepDetectorWithOptions::addOptions(
    boost::program_options::options_description& opt) const {
  ActsExamples::Options::addDD4hepOptions(opt);
}

auto DD4hepDetectorWithOptions::finalize(
    const boost::program_options::variables_map& vm,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  // read the detector config & dd4hep detector
  auto dd4HepDetectorConfig =
      ActsExamples::Options::readDD4hepConfig<po::variables_map>(vm);
  return m_detector.finalize(dd4HepDetectorConfig, mdecorator);
}

}  // namespace ActsExamples
