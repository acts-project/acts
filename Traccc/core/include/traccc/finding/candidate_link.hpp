// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/utils/pair.hpp"

namespace traccc {

// A link that contains the index of corresponding measurement and the index of
// a link from a previous step of track finding
struct candidate_link {
  // Step on which this link was found
  unsigned int step;

  // Index of the previous candidate
  unsigned int previous_candidate_idx;

  // Measurement index
  unsigned int meas_idx;

  // Index to the initial seed
  unsigned int seed_idx;

  // How many times it skipped a surface
  unsigned int n_skipped;

  // Number of consecutive holes; reset on measurement
  unsigned int n_consecutive_skipped;

  // chi2
  traccc::scalar chi2;

  // chi2 sum
  traccc::scalar chi2_sum;

  // degrees of freedom
  unsigned int ndf_sum;
};

}  // namespace traccc
