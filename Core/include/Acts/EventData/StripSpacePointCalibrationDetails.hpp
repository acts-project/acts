// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <array>

namespace Acts {

struct StripSpacePointCalibrationDetails final {
  std::array<float, 3> outerCenter;
  std::array<float, 3> outerToInnerGapVector;
  std::array<float, 3> outerHalfVector;
  std::array<float, 3> innerHalfVector;
};

struct StripSpacePointCalibrationDetailsDerived final {
  std::array<float, 3> outerToInnerGapCrossInnerHalfVector;
  std::array<float, 3> outerToInnerGapCrossOuterHalfVector;
  std::array<float, 3> innerCrossOuterHalfVector;
  std::array<float, 3> outerCenter;
  std::array<float, 3> outerHalfVector;
};

}  // namespace Acts
