// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsAlignment/Kernel/TrackAlignmentState.hpp"

#include <Eigen/LU>

namespace ActsAlignment {

void finaliseTrackAlignState(TrackAlignmentState& alignState) {
  alignState.chi2 = alignState.residual.transpose() *
                    alignState.measurementCovariance.inverse() *
                    alignState.residual;
  alignState.alignmentToChi2Derivative =
      Acts::DynamicVector::Zero(alignState.alignmentDof);
  alignState.alignmentToChi2SecondDerivative = Acts::DynamicMatrix::Zero(
      alignState.alignmentDof, alignState.alignmentDof);
  // The covariance of residual
  alignState.residualCovariance = Acts::DynamicMatrix::Zero(
      alignState.measurementDim, alignState.measurementDim);
  alignState.residualCovariance = alignState.measurementCovariance -
                                  alignState.projectionMatrix *
                                      alignState.trackParametersCovariance *
                                      alignState.projectionMatrix.transpose();

  alignState.alignmentToChi2Derivative =
      2 * alignState.alignmentToResidualDerivative.transpose() *
      alignState.measurementCovariance.inverse() *
      alignState.residualCovariance *
      alignState.measurementCovariance.inverse() * alignState.residual;
  alignState.alignmentToChi2SecondDerivative =
      2 * alignState.alignmentToResidualDerivative.transpose() *
      alignState.measurementCovariance.inverse() *
      alignState.residualCovariance *
      alignState.measurementCovariance.inverse() *
      alignState.alignmentToResidualDerivative;
}

}  // namespace ActsAlignment
