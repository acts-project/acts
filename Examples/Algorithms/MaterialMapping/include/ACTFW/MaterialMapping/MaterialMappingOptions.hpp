// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Utilities/Units.hpp>
#include <iostream>

#include "ACTFW/Utilities/Options.hpp"

namespace po = boost::program_options;

namespace FW {

namespace Options {

/// @brief Material mapping options, specially added
///
/// @tparam aopt_t Type of the options object (API bound to boost)
///
/// @param [in] opt_t The options object where the specific digitization
/// options are attached to
template <typename aopt_t>
void addMaterialMappingOptions(aopt_t& opt) {
  opt.add_options()("mat-mapping-collection",
                    po::value<std::string>()->default_value("material-tracks"),
                    "Collection name of the material tracks for reading.")(
      "mat-mapping-emptybins", po::value<bool>()->default_value(true),
      "Empty bin correction (recommended). Corrects for vaccuum/emtpy "
      "assigments.");
}

}  // namespace Options
}  // namespace FW
