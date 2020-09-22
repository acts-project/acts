// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/DigitizationOptions.hpp"
#include "ActsExamples/Digitization/HitSmearers.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include <boost/program_options.hpp>

#include <string>

void ActsExamples::Options::addDigitizationOptions(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::value;
  using namespace std;

  auto opt = desc.add_options();
  opt("digi-input-hits", value<string>()->default_value(""),
      "Name of the input hit collection.");
  opt("digi-output-measurements", value<string>()->default_value(""),
      "Name of the output measurement collection.");
  opt("digi-config-file", value<string>()->default_value(""),
      "Configuration (.json) file for digitization description, overwrites "
      "options input on command line.");
  opt("digi-smear-volume-name",
      value<read_strings>()->multitoken()->default_value({}),
      "Input: sensitive volume names.");
  opt("digi-smear-dimensions",
      value<read_series>()->multitoken()->default_value({}),
      "Input: smear parameters for this volume.");
  opt("digi-smear-binvalues",
      value<read_series>()->multitoken()->default_value({}),
      "Input: smear binning values for this volume.");
  opt("digi-smear-sigmas", value<read_range>()->multitoken()->default_value({}),
      "Input: smear sigma values for this volume.");
}

ActsExamples::SmearingAlgorithm::Config
ActsExamples::Options::readDigitizationConfig(
    const boost::program_options::variables_map& variables) {
  using namespace Acts::UnitLiterals;
  SmearingAlgorithm::Config smearCfg;

  return smearCfg;
}
