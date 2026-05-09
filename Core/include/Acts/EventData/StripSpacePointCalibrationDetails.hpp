// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <Eigen/Dense>

namespace Acts {

struct StripSpacePointCalibrationDetails final {
  Eigen::Vector3f innerCrossOuterStripHalfVector;
  Eigen::Vector3f stripSeparationCrossOuterHalfVector;
  Eigen::Vector3f stripSeparationCrossInnerHalfVector;
  // TODO use inner?
  Eigen::Vector3f outerStripCenter;
  Eigen::Vector3f outerStripHalfVector;
};

}  // namespace Acts
