// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Options/VertexingOptions.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <fstream>
#include <numeric>
#include <string>

#include <boost/program_options.hpp>

void ActsExamples::Options::addVertexingOptions(Description& desc) {
  using boost::program_options::bool_switch;
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("vertexing-eta-max", value<double>()->default_value(2.5),
      "Input track selection cut for primary vertexing - maximum absolute eta");
  opt("vertexing-rho-max", value<double>()->default_value(4.),
      "Input track selection cut for primary vertexing - maximum rho in mm");
  opt("vertexing-pt-min", value<double>()->default_value(500.),
      "Input track selection cut for primary vertexing - minimum pt in MeV");
}
