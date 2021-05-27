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

  // each volume configuration is one logical block
  //
  //   --digi-smear-volume-id=8
  //   --digi-smear-indices=0:1:5 # loc0, loc1, and time
  //   --digi-smear-types=0:0:3   # loc{0,1} uses gaussian, time uses uniform
  //   # parameter 0: loc0 gaussian width
  //   # parameter 1: loc1 gaussian width
  //   # parameter 2-4: time pitch,min,max
  //   --digi-smear-parameters=10:20:2.5:-25:25
  //
  // which can be repeated as often as needed
  //
  //   --digi-smear-volume-id=11
  //   --digi-smear-indices=1       # loc1
  //   --digi-smear-types=0         # loc1 uses gaussian
  //   --digi-smear-parameters=12.5 # loc1 gaussian width
  //
  auto opt = desc.add_options();
  opt("digi-config-file", value<std::string>()->default_value(""),
      "Configuration (.json) file for digitization description, overwrites "
      "smearing options input on command line.");
  opt("dump-digi-config", value<std::string>()->default_value(""),
      "Path to .json file in which to dump digitization configuration.");
  opt("digi-smear", bool_switch(), "Smearing: Switching hit smearing on");
  opt("digi-smear-volume", value<std::vector<int>>(),
      "Smearing Input: sensitive volume identifiers.");
  opt("digi-smear-indices", value<std::vector<VariableIntegers>>(),
      "Smearing Input: smear parameter indices for this volume.");
  opt("digi-smear-types", value<std::vector<VariableIntegers>>(),
      "Smearing Input: smear function types as 0 (gauss), 1 (truncated gauss), "
      "2 (clipped gauss), 3 (uniform), 4 (digital).");
  opt("digi-smear-parameters", value<std::vector<VariableReals>>(),
      "Smearing Input: smear parameters depending on the smearing type, 1 "
      "parameter for simple gauss, 3 for all others (1 parameter, 2 range "
      "values.");
  opt("digi-merge", bool_switch(), "Turn on hit merging");
  opt("digi-merge-nsigma", value<double>()->default_value(1.0),
      "Defines how close smeared parameters have to be when merging");
  opt("digi-merge-common-corner", bool_switch(),
      "Merge clusters which share a corner (8-cell connectivity)");
}
