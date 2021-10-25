// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/DigitizationOptions.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/Digitization/SmearingConfig.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <numeric>
#include <string>

#include <boost/program_options.hpp>

void ActsExamples::Options::addDigitizationOptions(Description& desc) {
  using boost::program_options::bool_switch;
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("digi-config-file", value<std::string>()->default_value(""),
      "Configuration (.json) file for digitization description, overwrites "
      "smearing options input on command line.");
  opt("dump-digi-config", value<std::string>()->default_value(""),
      "Path to .json file in which to dump digitization configuration.");
  opt("digi-merge", bool_switch(), "Turn on hit merging");
  opt("digi-merge-nsigma", value<double>()->default_value(1.0),
      "Defines how close smeared parameters have to be when merging");
  opt("digi-merge-common-corner", bool_switch(),
      "Merge clusters which share a corner (8-cell connectivity)");
}
