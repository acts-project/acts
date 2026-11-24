// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <string_view>

namespace Acts::Experimental::detail {

/// Why a check rejected a candidate, filled in by the check for the caller
/// to record.
struct GbtsCutFailure {
  /// The cut it failed, one of the names in `GbtsProfilerTypes.hpp`.
  std::string_view cut;
  /// The value cut on.
  float value{};
  /// The threshold it failed, zero for a cut without one.
  float threshold{};
};

}  // namespace Acts::Experimental::detail
