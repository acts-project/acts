// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsAlignment/Kernel/detail/AlignmentEngine.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"

#include <cstddef>

namespace ActsAlignment::detail {

void resetAlignmentDerivative(Acts::AlignmentToBoundMatrix& alignToBound,
                              AlignmentMask mask) {
  if (!ACTS_CHECK_BIT(mask, AlignmentMask::Center0)) {
    alignToBound.col(Acts::eAlignmentCenter0) = Acts::AlignmentVector::Zero();
  }
  if (!ACTS_CHECK_BIT(mask, AlignmentMask::Center1)) {
    alignToBound.col(Acts::eAlignmentCenter1) = Acts::AlignmentVector::Zero();
  }
  if (!ACTS_CHECK_BIT(mask, AlignmentMask::Center2)) {
    alignToBound.col(Acts::eAlignmentCenter2) = Acts::AlignmentVector::Zero();
  }
  if (!ACTS_CHECK_BIT(mask, AlignmentMask::Rotation0)) {
    alignToBound.col(Acts::eAlignmentRotation0) = Acts::AlignmentVector::Zero();
  }
  if (!ACTS_CHECK_BIT(mask, AlignmentMask::Rotation1)) {
    alignToBound.col(Acts::eAlignmentRotation1) = Acts::AlignmentVector::Zero();
  }
  if (!ACTS_CHECK_BIT(mask, AlignmentMask::Rotation2)) {
    alignToBound.col(Acts::eAlignmentRotation2) = Acts::AlignmentVector::Zero();
  }
}

void finaliseTrackAlignState(TrackAlignmentState& alignState) {
  // Calculate the chi2 and chi2 derivatives based on the alignment matrixs
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

void solveAlignmentParameters(
    const std::vector<TrackAlignmentState>& trackAlignmentStates,
    AlignmentResult& alignResult, const Acts::Logger& logger) {
  // The total alignment degree of freedom
  alignResult.alignmentDof =
      alignResult.idxedAlignSurfaces.size() * Acts::eAlignmentSize;
  // Initialize derivative of chi2 w.r.t. alignment parameters for all tracks
  alignResult.sumChi2Derivative =
      Acts::DynamicVector::Zero(alignResult.alignmentDof);
  alignResult.sumChi2SecondDerivative = Acts::DynamicMatrix::Zero(
      alignResult.alignmentDof, alignResult.alignmentDof);
  alignResult.chi2 = 0;
  alignResult.measurementDim = 0;
  alignResult.numTracks = trackAlignmentStates.size();
  double sumChi2ONdf = 0;
  for (const auto& alignState : trackAlignmentStates) {
    for (const auto& [rowSurface, rows] : alignState.alignedSurfaces) {
      const auto& [dstRow, srcRow] = rows;
      // Fill the results into full chi2 derivative matrix
      alignResult.sumChi2Derivative.segment<Acts::eAlignmentSize>(
          dstRow * Acts::eAlignmentSize) +=
          alignState.alignmentToChi2Derivative.segment(
              srcRow * Acts::eAlignmentSize, Acts::eAlignmentSize);

      for (const auto& [colSurface, cols] : alignState.alignedSurfaces) {
        const auto& [dstCol, srcCol] = cols;
        alignResult.sumChi2SecondDerivative
            .block<Acts::eAlignmentSize, Acts::eAlignmentSize>(
                dstRow * Acts::eAlignmentSize, dstCol * Acts::eAlignmentSize) +=
            alignState.alignmentToChi2SecondDerivative.block(
                srcRow * Acts::eAlignmentSize, srcCol * Acts::eAlignmentSize,
                Acts::eAlignmentSize, Acts::eAlignmentSize);
      }
    }
    alignResult.chi2 += alignState.chi2;
    alignResult.measurementDim += alignState.measurementDim;
    sumChi2ONdf += alignState.chi2 / alignState.measurementDim;
  }
  alignResult.averageChi2ONdf = sumChi2ONdf / alignResult.numTracks;

  // Get the inverse of chi2 second derivative matrix (we need this to
  // calculate the covariance of the alignment parameters)
  // @TODO: use more stable method for solving the inverse
  std::size_t alignDof = alignResult.alignmentDof;
  Acts::DynamicMatrix sumChi2SecondDerivativeInverse =
      Acts::DynamicMatrix::Zero(alignDof, alignDof);
  sumChi2SecondDerivativeInverse =
      alignResult.sumChi2SecondDerivative.inverse();
  if (sumChi2SecondDerivativeInverse.hasNaN()) {
    ACTS_LOG_WITH_LOGGER(logger, Acts::Logging::DEBUG,
                         "Chi2 second derivative inverse has NaN");
  }

  // Initialize the alignment results
  alignResult.deltaAlignmentParameters = Acts::DynamicVector::Zero(alignDof);
  alignResult.alignmentCovariance =
      Acts::DynamicMatrix::Zero(alignDof, alignDof);
  // Solve the linear equation to get alignment parameters change
  alignResult.deltaAlignmentParameters =
      -alignResult.sumChi2SecondDerivative.fullPivLu().solve(
          alignResult.sumChi2Derivative);
  ACTS_LOG_WITH_LOGGER(logger, Acts::Logging::VERBOSE,
                       "sumChi2SecondDerivative = \n"
                           << alignResult.sumChi2SecondDerivative);
  ACTS_LOG_WITH_LOGGER(logger, Acts::Logging::VERBOSE,
                       "sumChi2Derivative = \n"
                           << alignResult.sumChi2Derivative);
  ACTS_LOG_WITH_LOGGER(logger, Acts::Logging::VERBOSE,
                       "alignResult.deltaAlignmentParameters \n");

  // Alignment parameters covariance
  alignResult.alignmentCovariance = 2 * sumChi2SecondDerivativeInverse;
  // chi2 change
  alignResult.deltaChi2 = 0.5 * alignResult.sumChi2Derivative.transpose() *
                          alignResult.deltaAlignmentParameters;
}

}  // namespace ActsAlignment::detail
