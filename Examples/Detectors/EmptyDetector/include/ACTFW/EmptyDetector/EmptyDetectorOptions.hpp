// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/program_options.hpp>
#include <utility>

#include "Acts/Utilities/Units.hpp"

namespace FW {

namespace Options {

using namespace Acts::UnitLiterals;
namespace po = boost::program_options;

/// The empty geometry options, the are prefixes with geo-generic
///
/// @tparam options_t Type of the options object (bound to boost API)
///
/// @param opt The provided object, where root specific options are attached
template <typename options_t>
void addEmptyGeometryOptions(options_t& opt) {
  opt.add_options()("geo-empty-radius", po::value<double>()->default_value(2_m),
                    "Radius of the empty cylinder [in m]. ")(
      "geo-empty-halfLength", po::value<double>()->default_value(10_m),
      "Half length of the empty cylinder [in m].");
}
}  // namespace Options
}  // namespace FW
