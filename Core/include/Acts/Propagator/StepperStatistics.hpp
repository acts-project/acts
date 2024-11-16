// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>

namespace Acts {

struct StepperStatistics {
  std::size_t nAttemptedSteps = 0;
  std::size_t nFailedSteps = 0;
  std::size_t nSuccessfulSteps = 0;
  std::size_t nReverseSteps = 0;

  double pathLength = 0;
  double absolutePathLength = 0;
};

}  // namespace Acts
