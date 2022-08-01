// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Io/Svg/SvgTrackingGeometryWriter.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"

#include <iostream>

namespace ActsExamples {

namespace Options {

/// Common Svg writing options
//////
/// @param opt The options object, where string based options are attached
void addSvgWriterOptions(boost::program_options::options_description& opt) {
  namespace po = boost::program_options;
  opt.add_options()("svg-base", po::value<std::string>()->default_value(""),
                    "Base name of the output files.");
}

/// read the evgen options and return a Config file
ActsExamples::SvgTrackingGeometryWriter::Config
readSvgTrackingGeometryWriterConfig(
    const boost::program_options::variables_map& vm) {
  namespace po = boost::program_options;
  ActsExamples::SvgTrackingGeometryWriter::Config svgTgConfig;
  svgTgConfig.baseName = vm["svg-base"].as<std::string>();
  return svgTgConfig;
}

}  // namespace Options
}  // namespace ActsExamples
