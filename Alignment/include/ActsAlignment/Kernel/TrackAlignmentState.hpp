// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <unordered_map>
#include <utility>

namespace Acts {
class Surface;
}

namespace ActsAlignment {

class AlignableStructure;

/// A single track's contribution to alignment. Measurement rows follow the
/// same order in the residual, covariance, projection, and derivative matrices.
struct TrackAlignmentState {
  /// Number of scalar measurements.
  std::size_t measurementDim = 0;

  /// Number of track parameters represented in the projection matrix.
  std::size_t trackParametersDim = 0;

  /// Number of alignment parameters represented on this track.
  std::size_t alignmentDof = 0;

  /// Covariance of the measurements (measurementDim by measurementDim).
  Acts::DynamicMatrix measurementCovariance;

  /// Fitted track covariance (trackParametersDim by trackParametersDim).
  Acts::DynamicMatrix trackParametersCovariance;

  /// Projects track parameters into measurements (measurementDim by
  /// trackParametersDim).
  Acts::DynamicMatrix projectionMatrix;

  /// Measurement minus prediction, with one entry per scalar measurement.
  Acts::DynamicVector residual;

  /// Covariance of the residuals (measurementDim by measurementDim).
  Acts::DynamicMatrix residualCovariance;

  /// Track chi-squared.
  double chi2 = 0;

  /// Residual derivatives (measurementDim by alignmentDof).
  Acts::DynamicMatrix alignmentToResidualDerivative;

  /// Gradient of track chi-squared with respect to the alignment parameters.
  Acts::DynamicVector alignmentToChi2Derivative;

  /// Hessian of track chi-squared (alignmentDof by alignmentDof).
  Acts::DynamicMatrix alignmentToChi2SecondDerivative;

  /// Global and track-local surface indices, respectively. Each index
  /// identifies a block of six alignment parameters. Surface pointers are
  /// non-owning.
  std::unordered_map<const Acts::Surface*, std::pair<std::size_t, std::size_t>>
      alignedSurfaces;

  /// Immediate owning structure for each grouped surface, when a hierarchy is
  /// used. Both pointers are non-owning.
  std::unordered_map<const Acts::Surface*, const AlignableStructure*>
      surfaceStructures;
};

}  // namespace ActsAlignment
