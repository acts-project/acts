// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFitting/TrackFittingOptions.hpp"

using namespace boost::program_options;

void ActsExamples::Options::addFittingOptions(
    boost::program_options::options_description& opt) {
  opt.add_options()("fit-directed-navigation", bool_switch(),
                    "Fit tracks with DirectNavigator");
  opt.add_options()("fit-pick-track", value<int>()->default_value(-1),
                    "Pick a single track by track number (-1 for all tracks)");
}
