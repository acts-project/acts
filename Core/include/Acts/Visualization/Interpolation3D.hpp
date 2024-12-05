// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <unsupported/Eigen/Splines>

namespace Acts::Interpolation3D {

/// @brief Helper function to interpolate points using a spline
/// from Eigen
///
/// @tparam input_vector_type
/// @param inputs input vector points
/// @param nPoints number of interpolation points
///
/// @return std::vector<Acts::Vector3> interpolated points
template <typename input_vector_type>
std::vector<Acts::Vector3> spline(const std::vector<input_vector_type>& inputs,
                                  std::size_t nPoints) {
  std::vector<Acts::Vector3> output;

  if (nPoints < 2) {
    // No interpolation done return simply the output vector
    for (const auto& input : inputs) {
      output.push_back(input.template head<3>());
    }

  } else {
    Eigen::MatrixXd points(3, inputs.size());
    for (std::size_t i = 0; i < inputs.size(); ++i) {
      points.col(i) = inputs[i].template head<3>().transpose();
    }
    Eigen::Spline<double, 3> spline3D =
        Eigen::SplineFitting<Eigen::Spline<double, 3>>::Interpolate(points, 2);

    double step = 1. / (nPoints - 1);
    for (std::size_t i = 0; i < nPoints; ++i) {
      double t = i * step;
      output.push_back(spline3D(t));
    }
  }
  return output;
}

}  // namespace Acts::Interpolation3D
