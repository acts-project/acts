// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Mille/ActsToMille.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Alignment.hpp"
#include "ActsPlugins/Mille/Helpers.hpp"

#include <cstddef>
#include <iostream>
#include <numeric>

#include <Eigen/src/Core/Matrix.h>

namespace ActsPlugins::ActsToMille {

void dumpToMille(const ActsAlignment::detail::TrackAlignmentState& state,
                 MilleRecord* record) {
  // check for a valid record
  if (record == nullptr) {
    std::cerr << " Missing Mille record " << std::endl;
    return;
  }

  // prepare the vectors to interface to Mille
  std::vector<unsigned int> localIndices(state.trackParametersDim, 0);
  std::vector<int> globalIndices(state.alignmentDof, 0.);
  std::vector<double> localDeriv(state.trackParametersDim, 0.);
  std::vector<double> globalDeriv(state.alignmentDof, 0.);

  // prepare the track parameter index array (always the same)
  std::iota(localIndices.begin(), localIndices.end(), 1);

  // map the alignment parameter labels.
  // Important: Millepede expects indices to start with 1
  std::vector<std::pair<int, int>> aliParLocalToGlobal;
  for (auto& [surf, indices] : state.alignedSurfaces) {
    auto& [globalSurfIndex, localStartIndex] = indices;
    for (std::size_t iPar = 0; iPar < Acts::eAlignmentSize; ++iPar) {
      int dest = Acts::eAlignmentSize * globalSurfIndex + iPar + 1;
      int src = Acts::eAlignmentSize * localStartIndex + iPar;
      aliParLocalToGlobal.emplace_back(src, dest);
    }
  }

  /// 1) write out the local measurements on the surfaces and their direct
  /// derivatives. This will populate the upper / left three quadrants of the
  /// alignment matrix, including direct correlations between alignment and
  /// track parameters.
  /// TODO: Add explicit diagonalisation for correlated (stereo) measurements.
  for (std::size_t iMeas = 0; iMeas < state.measurementDim; ++iMeas) {
    // arrange the global parameters correctly
    for (auto& [srcGlobal, destGlobal] : aliParLocalToGlobal) {
      // index for each global derivative
      globalIndices[srcGlobal] = destGlobal;
      // value for each global derivative
      globalDeriv[srcGlobal] =
          state.alignmentToResidualDerivative(iMeas, srcGlobal);
    }
    // local derivatives due to measurement uncertainties
    for (std::size_t iTrkPar = 0; iTrkPar < state.trackParametersDim;
         ++iTrkPar) {
      localDeriv[iTrkPar] = state.projectionMatrix(iMeas, iTrkPar);
    }
    // write a measurement to the ongoing Mille record.
    record->addData(
        // residual
        state.residual(iMeas),
        // sigma
        std::sqrt(state.measurementCovariance(iMeas, iMeas)),
        // local parameter indices
        localIndices,
        // local derivatives
        localDeriv, globalIndices, globalDeriv);
  }

  /// 2) Write out additional pseudo-measurements representing the (local) track
  /// parameter correlations from the Kalman fit (linearisation point).
  /// These enter the bottom right quadrant of the alignment matrix and
  /// represent the additional contributions to the track fit chi2
  /// arising from the correlations between the parameters on different
  /// surfaces encoded in the Kalman track model.
  /// Step 1) already added the terms arising from measurement uncertainties,
  /// now we need to extend this to the full Kalman covariance.
  /// This requires "unfitting" the tracks, making this a rather slow
  /// and numerically tricky step.

  /// compute the part of the weight matrix arising from Step 1).
  /// This is already present in the Mille record and should not
  /// be duplicated
  auto weightMatMeasurements = state.projectionMatrix.transpose() *
                               state.measurementCovariance.inverse() *
                               state.projectionMatrix;

  // regularise the (full) Kalman covariance. This is needed to stabilise
  // poorly constrained directions (usually: time)
  Acts::DynamicMatrix regularisedCov =
      regulariseCovariance(state.trackParametersCovariance);

  // now we can get the piece of the weight matrix not already covered by
  // the measurement uncertainties
  Acts::DynamicMatrix correlationTerm =
      getInverseComplement(regularisedCov, weightMatMeasurements);

  // Decompose the matrix we need to add into a sum of rank-1 matrices,
  // C_add = sum (lambda_i v_i v_i^T), which can be interpreted
  // as pseudo-measurements with sigma_i 1/sqrt(lambda_i) and local derivatives
  // v_i. This relies on C_add being symmetric positive (semi)definite.
  auto eigensolver =
      Eigen::SelfAdjointEigenSolver<Acts::DynamicMatrix>(correlationTerm);
  if (eigensolver.info() != Eigen::Success) {
    std::cout << " FAILED to find decompose correlation term" << std::endl;
    return;
  }
  auto eigenVals = eigensolver.eigenvalues();
  auto eigenVecs = eigensolver.eigenvectors();

  // no dependence on global parameters - these terms only enter the
  // track covariance sub-matrix of the alignment problem (bottom right
  // quadrant)
  globalDeriv.clear();
  globalIndices.clear();

  /// convert each EV to a pseudo-measurement
  for (std::size_t iMeas = 0; iMeas < state.trackParametersDim; ++iMeas) {
    // fill the local derivatives from the current eigenvector
    for (std::size_t iTrkPar = 0; iTrkPar < localDeriv.size(); ++iTrkPar) {
      localDeriv[iTrkPar] = eigenVecs(iTrkPar, iMeas);
    }
    // and write a pseudo-measurement to Mille.
    record->addData(
        // residual == 0 for pseudo-measurements
        0,
        // sigma = sqrt(EV)
        1. / std::sqrt(eigenVals(iMeas)),
        // local parameter indices
        localIndices,
        // local derivatives
        localDeriv, globalIndices, globalDeriv);
  }
  // track is fully written - end the record in Mille
  record->writeRecord();
}
}  // namespace ActsPlugins::ActsToMille
