// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SpacePointMakerOptions.hpp"

#include "ActsExamples/Io/Json/JsonGeometryList.hpp"

#include <string>

#include <boost/program_options.hpp>

void ActsExamples::Options::addSpacePointMakerOptions(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("geo-selection-config-file", value<std::string>()->default_value(""),
      "Json file for space point geometry selection");
}

ActsExamples::SpacePointMaker::Config
ActsExamples::Options::readSpacePointMakerConfig(
    const ActsExamples::Options::Variables& variables) {
  SpacePointMaker::Config cfg;
  std::string path{variables["geo-selection-config-file"].as<std::string>()};
  if (not path.empty()) {
    cfg.geometrySelection = ActsExamples::readJsonGeometryList(path);
  }
  return cfg;
}
