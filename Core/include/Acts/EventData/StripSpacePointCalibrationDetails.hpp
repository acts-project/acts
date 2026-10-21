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

/// Collection of outer strip space point details, used for the calibration.
struct OuterStripSpacePointCalibrationDetails final {
  /// Center of the outer strip.
  std::array<float, 3> outerCenter;
  /// Separation vector from the inner strip center to the outer strip center.
  std::array<float, 3> innerToOuterSeparation;
  /// Half vector of the outer strip, pointing from the center to one end of the
  /// strip.
  std::array<float, 3> outerHalfVector;
  /// Half vector of the inner strip, pointing from the center to one end of the
  /// strip.
  std::array<float, 3> innerHalfVector;
};

/// Derived collection of outer strip space point details, used for the
/// calibration.
struct OuterStripSpacePointCalibrationDetailsDerived final {
  /// Cross product of the separation vector from the inner strip center to the
  /// outer strip center with the inner half vector.
  std::array<float, 3> innerToOuterSeparationCrossInnerHalfVector;
  /// Cross product of the separation vector from the inner strip center to the
  /// outer strip center with the outer half vector.
  std::array<float, 3> innerToOuterSeparationCrossOuterHalfVector;
  /// Cross product of the inner half vector with the outer half vector.
  std::array<float, 3> innerCrossOuterHalfVector;
  /// Center of the outer strip.
  std::array<float, 3> outerCenter;
  /// Half vector of the outer strip, pointing from the center to one end of
  /// the strip.
  std::array<float, 3> outerHalfVector;
};

}  // namespace Acts
