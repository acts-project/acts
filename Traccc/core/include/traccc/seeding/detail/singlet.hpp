// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/container.hpp"

namespace traccc {

/// location of spacepoint in internal spacepoint container
struct sp_location {
  /// index of the bin of the spacepoint grid
  unsigned int bin_idx;
  /// index of the spacepoint in the bin
  unsigned int sp_idx;
};

inline TRACCC_HOST_DEVICE bool operator==(const sp_location& lhs,
                                          const sp_location& rhs) {
  return (lhs.bin_idx == rhs.bin_idx && lhs.sp_idx == rhs.sp_idx);
}

inline TRACCC_HOST_DEVICE bool operator!=(const sp_location& lhs,
                                          const sp_location& rhs) {
  return (lhs.bin_idx != rhs.bin_idx || lhs.sp_idx != rhs.sp_idx);
}

}  // namespace traccc
