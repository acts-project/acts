// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/ScatteringTail.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace Acts {

ScatteringTailComponent siliconScatteringTail(
    const MaterialSlab& slab, const ParticleHypothesis& /*particleHypothesis*/,
    double /*qOverP*/) {
  // two-Gaussian EM fit to 2M samples per thickness, linear in ln(x/X0)
  static constexpr std::array<double, 6> lnThickness = {
      -5.8091, -5.2983, -4.6052, -3.9120, -3.2189, -2.5257};
  static constexpr std::array<double, 6> weight = {0.0152, 0.0151, 0.0156,
                                                   0.0171, 0.0199, 0.0276};
  static constexpr std::array<double, 6> widthScale = {7.263, 6.945, 6.207,
                                                       5.373, 4.554, 3.628};
  const double x =
      std::clamp(std::log(static_cast<double>(slab.thicknessInX0())),
                 lnThickness.front(), lnThickness.back());
  std::size_t i = 0;
  while (i + 2 < lnThickness.size() && x > lnThickness[i + 1]) {
    ++i;
  }
  const double f = (x - lnThickness[i]) / (lnThickness[i + 1] - lnThickness[i]);
  return {
      .weight = weight[i] + f * (weight[i + 1] - weight[i]),
      .widthScale = widthScale[i] + f * (widthScale[i + 1] - widthScale[i])};
}

}  // namespace Acts
