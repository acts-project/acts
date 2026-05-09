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
  // TODO use inner?
  std::array<float, 3> outerStripCenter;
  std::array<float, 3> stripSeparation;
  std::array<float, 3> outerStripHalfVector;
  std::array<float, 3> innerStripHalfVector;
};

struct StripSpacePointCalibrationDetailsDerived final {
  std::array<float, 3> stripSeparationCrossInnerHalfVector;
  std::array<float, 3> stripSeparationCrossOuterHalfVector;
  std::array<float, 3> innerCrossOuterStripHalfVector;
  // TODO use inner?
  std::array<float, 3> outerStripCenter;
  std::array<float, 3> outerStripHalfVector;
};

}  // namespace Acts
