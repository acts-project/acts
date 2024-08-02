// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>

namespace Acts {

struct HashingTrainingConfig {
  /// Random seed for Annoy
  unsigned int AnnoySeed = 123456789;

  /// Number of features to use
  std::int32_t f = 1;
};

}  // namespace Acts
