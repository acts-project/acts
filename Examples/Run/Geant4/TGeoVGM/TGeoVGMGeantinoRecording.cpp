// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Options/CommonOptions.hpp"
#include "../GeantinoRecordingBase.hpp"
#include "TGeoVGMDetectorConstruction.hpp"

#include <boost/program_options.hpp>

using namespace ActsExamples;

int main(int argc, char* argv[]) {
  // Setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addOutputOptions(desc);
  Options::addGeant4Options(desc);
  desc.add_options()(
      "geo-tgeo-filename",
      boost::program_options::value<std::string>()->default_value(""),
      "Root file name.");
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  auto g4detector = std::make_unique<TGeoVGMDetectorConstruction>(
      vm["geo-tgeo-filename"].as<std::string>());

  return runGeantinoRecording(vm, std::move(g4detector));
}
