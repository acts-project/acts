// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <utility>

#include <boost/program_options.hpp>

namespace ActsExamples {

namespace Options {

namespace po = boost::program_options;

/// The telescope geometry options
///
/// @tparam options_t Type of the options object (bound to boost API)
///
/// @param opt The provided object, where root specific options are attached
template <typename options_t>
void addTelescopeGeometryOptions(options_t& opt) {
  opt.add_options()(
      "geo-tele-longshifts",
      po::value<read_range>()->multitoken()->default_value(
          {0, 30, 60, 120, 150, 180}),
      "Telescope detector Input: the layers shifts in the longidutinal "
      "direction [mm], i.e. shifts along x, y or z axis")(
      "geo-tele-transhifts",
      po::value<read_range>()->multitoken()->default_value({0, 0}),
      "Telescope detector Input: the layers shifts in the transverse plane "
      "[mm], i.e. shifts along (y, z), (x, z) or (x, y) axies. Same values for "
      "all layers")(
      "geo-tele-layerbounds",
      po::value<read_range>()->multitoken()->default_value({30, 15}),
      "Telescope detector Input: the half-x and half-y of each layer size "
      "[mm]")("geo-tele-matthickness", po::value<double>()->default_value(80),
              "Telescope detector Input: the silicon material thickness of "
              "each layer [um]. Same value for all layers")(
      "geo-tele-alignaxis", po::value<size_t>()->default_value(2),
      "Telescope detector Input: the detector is aligned along which "
      "direction: 0 - x axis, 1 - y axis, 2 - z aixs");
}
}  // namespace Options
}  // namespace ActsExamples
