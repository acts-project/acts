// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
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

/// The generic geometry options, the are prefixes with geo-generic
///
/// @tparam options_t Type of the options object (bound to boost API)
///
/// @param opt The provided object, where root specific options are attached
template <typename options_t>
void addGenericGeometryOptions(options_t& opt) {
  opt.add_options()("geo-generic-buildlevel",
                    po::value<size_t>()->default_value(3),
                    "The building level: 0 - pixel barrel only, 1 - pixel "
                    "detector only, 2 - full barrel only, 3 - full detector "
                    "(without stereo modules).");
}
}  // namespace Options
}  // namespace ActsExamples
