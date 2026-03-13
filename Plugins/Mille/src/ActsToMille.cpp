// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Mille/ActsToMille.hpp"

namespace ActsPlugins::ActsToMille {

/// this is a placeholder, which will be replaced by actual functionality
void dumpToMille(const ActsAlignment::detail::TrackAlignmentState&,
                 MilleRecord* record) {
  // call into Mille to test for linkage
  if (record) {
    // write a single dummy measurement to test the I/O functionality
    record->addData(
        -0.005f,                       // measurement
        0.008f,                        // uncertainty
        std::vector<unsigned int>{1},  // indices of local parameters
        std::vector<double>{1.},       // derivates by local parameters
        std::vector<int>{42},          // indices of alignment parameters
        std::vector<double>{-1.});     // derivatives by alignment parameters
  }
}
}  // namespace ActsPlugins::ActsToMille
