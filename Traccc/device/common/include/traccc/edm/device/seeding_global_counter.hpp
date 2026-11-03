// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace traccc::device {

/// Total number of doublets and triplets found in seeding
struct seeding_global_counter {
  /// The total number of middle-bottom doublets
  unsigned int m_nMidBot;

  /// The total number of middle-top doublets
  unsigned int m_nMidTop;

  /// The total number of triplets
  unsigned int m_nTriplets;

};  // struct seeding_global_counter

}  // namespace traccc::device
