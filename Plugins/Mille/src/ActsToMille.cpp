// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Mille/ActsToMille.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "ActsPlugins/Mille/Helpers.hpp"

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <map>
#include <numeric>
#include <set>

#include <Eigen/src/Core/Matrix.h>
#include <Mille/MilleDataStructures.h>

#include "Mille/MilleDecoder.h"

namespace ActsPlugins::ActsToMille {

namespace {

/// for the ACTS-internal indices, start counting at 0
unsigned long internalIndexSurfToParam(unsigned long surfaceIndex,
                                       unsigned long dofIndex) {
  return surfaceIndex * Acts::eAlignmentSize + dofIndex;
}
/// for global alignment parameters, start counting at 1 (fotran convention) for
/// Mille
unsigned long globalIndexSurfToParam(unsigned long surfaceIndex,
                                     unsigned long dofIndex) {
  return surfaceIndex * Acts::eAlignmentSize + dofIndex + 1;
}

}  // namespace

void dumpToMille(const ActsAlignment::detail::TrackAlignmentState& state,
                 MilleRecord& record, bool removeUnconstrainedTrackPar) {
  // prepare the vectors to interface to Mille
  std::vector<int> globalIndices(state.alignmentDof, 0.);
  std::vector<double> globalDeriv(state.alignmentDof, 0.);

  // map the alignment parameter labels.
  // Important: Millepede expects indices to start with 1
  std::vector<std::pair<int, int>> aliParLocalToGlobal;
  for (auto& [surf, indices] : state.alignedSurfaces) {
    auto& [globalSurfIndex, localStartIndex] = indices;
    for (std::size_t iPar = 0; iPar < Acts::eAlignmentSize; ++iPar) {
      aliParLocalToGlobal.emplace_back(
          internalIndexSurfToParam(localStartIndex, iPar),
          globalIndexSurfToParam(globalSurfIndex, iPar));
    }
  }

  // Analyse the track fit, and potentially discard unconstrained
  // parameters to stabilise the system
  std::set<std::size_t> skippedTrackParams = {};
  if (removeUnconstrainedTrackPar) {
    // collect (sorted by covariance) the indices of all track parameters
    std::map<double, int> trkParByCov;
    for (std::size_t k = 0; k < state.trackParametersDim; ++k) {
      trkParByCov.emplace(state.trackParametersCovariance(k, k), k);
    }

    // now, loop through the parameter list and look for huge jumps.
    double prev = 0;
    for (auto& [sigmaSquared, index] : trkParByCov) {
      // a jump of 1e6 is indicative that we are not in Kansas anymore
      if (prev != 0 && sigmaSquared > 1e6 * prev) {
        skippedTrackParams.insert(index);
      } else {
        prev = sigmaSquared;
      }
    }
  }
  // calculate number of remaining active track parameters
  std::size_t effectiveTrackParDim =
      state.trackParametersDim - skippedTrackParams.size();

  std::vector<unsigned int> localIndices(effectiveTrackParDim, 0);
  std::vector<double> localDeriv(effectiveTrackParDim, 0.);
  // prepare the track parameter index array (always the same)
  std::iota(localIndices.begin(), localIndices.end(), 1);

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
    // index that we show to Mille (differs from internal index if we
    // skip any parameters)
    std::size_t iPar = 0;
    // local derivatives due to measurement uncertainties
    for (std::size_t iTrkPar = 0; iTrkPar < state.trackParametersDim;
         ++iTrkPar) {
      if (skippedTrackParams.contains(iTrkPar)) {
        continue;
      }
      localDeriv[iPar] = state.projectionMatrix(iMeas, iTrkPar);
      ++iPar;
    }
    // write a measurement to the ongoing Mille record.
    record.addData(
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

  // for this, we need to calculate the "reduced" covariance
  // and projection matrices, accounting for skipped parameters
  Acts::DynamicMatrix updatedCovariance{effectiveTrackParDim,
                                        effectiveTrackParDim};
  Acts::DynamicMatrix updatedProjection{state.measurementDim,
                                        effectiveTrackParDim};
  std::size_t milleTPindex = 0;
  for (std::size_t internalTPindex = 0;
       internalTPindex < state.trackParametersDim; ++internalTPindex) {
    if (skippedTrackParams.contains(internalTPindex)) {
      continue;
    }

    // update projection matrix
    updatedProjection.col(milleTPindex) =
        state.projectionMatrix.col(internalTPindex);

    std::size_t secondMilleIndex = 0;
    for (std::size_t secondInternalIndex = 0;
         secondInternalIndex < state.trackParametersDim;
         ++secondInternalIndex) {
      if (skippedTrackParams.contains(secondInternalIndex)) {
        continue;
      }
      updatedCovariance(milleTPindex, secondMilleIndex) =
          state.trackParametersCovariance(internalTPindex, secondInternalIndex);
      ++secondMilleIndex;
    }
    ++milleTPindex;
  }

  const Acts::DynamicMatrix weightMatMeasurements =
      updatedProjection.transpose() * state.measurementCovariance.inverse() *
      updatedProjection;

  // regularise the (full) Kalman covariance. This is needed to stabilise
  // poorly constrained directions
  double regCondCutOff = 1e-10;
  double regHugeLeading = 100.;
  double regOnDiag = 1.e-10;
  if (removeUnconstrainedTrackPar) {
    // if we trim the poorly constrained directions ahead of time,
    // we can be a bit less aggressive in the regularisation
    regCondCutOff = -1;
    regHugeLeading = -1.;
    regOnDiag = 1.e-9;
  }
  const Acts::DynamicMatrix regularisedCov = regulariseCovariance(
      updatedCovariance, regCondCutOff, regHugeLeading, regOnDiag);

  // now we can get the piece of the weight matrix not already covered by
  // the measurement uncertainties
  const Acts::DynamicMatrix correlationTerm =
      getInverseComplement(regularisedCov, weightMatMeasurements);

  // Decompose the matrix we need to add into a sum of rank-1 matrices,
  // C_add = sum (lambda_i v_i v_i^T), which can be interpreted
  // as pseudo-measurements with sigma_i 1/sqrt(lambda_i) and local derivatives
  // v_i. This relies on C_add being symmetric positive (semi)definite.
  Eigen::SelfAdjointEigenSolver<Acts::DynamicMatrix> eigenSolver(
      correlationTerm);
  if (eigenSolver.info() != Eigen::Success) {
    std::cout << " FAILED to find decompose correlation term" << std::endl;
    return;
  }
  const Acts::DynamicVector eigenVals = eigenSolver.eigenvalues();
  const Acts::DynamicMatrix eigenVecs = eigenSolver.eigenvectors();

  // no dependence on global parameters - these terms only enter the
  // track covariance sub-matrix of the alignment problem (bottom right
  // quadrant)
  globalDeriv.clear();
  globalIndices.clear();

  /// convert each EV to a pseudo-measurement
  for (long iMeas = 0; iMeas < eigenVecs.rows();
       ++iMeas) {  // fill the local derivatives from the current eigenvector
    // skip negative EV
    if (eigenVals(iMeas) <= 0) {
      continue;
    }
    for (std::size_t iPar = 0; iPar < effectiveTrackParDim; ++iPar) {
      localDeriv[iPar] = eigenVecs(iPar, iMeas);
    }
    // and write a pseudo-measurement to Mille.
    record.addData(
        // residual == 0 for pseudo-measurements
        0,
        // EV == weight = 1/sigma^2
        1. / std::sqrt(eigenVals(iMeas)),
        // local parameter indices
        localIndices,
        // local derivatives
        localDeriv, globalIndices, globalDeriv);
  }
  // track is fully written - end the record in Mille
  record.writeRecord();
}

Mille::MilleDecoder::ReadResult unpackMilleRecord(
    Mille::IMilleReader& reader,
    ActsAlignment::detail::TrackAlignmentState& targetState,
    const std::unordered_map<const Acts::Surface*, std::size_t>&
        idxedAlignSurfaces) {
  /// book a decoder
  Mille::MilleDecoder decoder;
  // vector to hold the extracted measurements
  std::vector<Mille::MilleMeasurement> measurements;
  // attempt to decode the next record from the binary
  auto res = decoder.decode(reader, measurements);

  // if we are EoF or encountered an error, return the result.
  if (res != Mille::MilleDecoder::ReadResult::OK || measurements.empty()) {
    return res;
  }

  // initialise the components of the target state which we will touch
  targetState.measurementDim = 0;

  // Still here - we got a valid record! Let's write it into our target state.
  // In the following, we emulate what MillePede-II is doing internally.
  // This is somewhat approximate, as we do not run any of the cleaning /
  // conditioning performed by MP-II.

  // Step 1: Parameter discovery
  // Goal: Identify all existing parameters and assign internal indices.
  int firstLocal = 99999;
  int lastLocal = 0;
  std::set<int> seenGlobalLabels;
  std::set<int> seenSurfaceLabels;
  // a measurement with a residual of identical zero is typically a
  // correlation constraint encoded as pseudo-measurement
  targetState.measurementDim =
      std::count_if(measurements.begin(), measurements.end(),
                    [](const Mille::MilleMeasurement& measurement) {
                      return measurement.measurement != 0;
                    });

  // discover labels in use
  for (const Mille::MilleMeasurement& measurement : measurements) {
    if (measurement.localLabels.empty()) {
      continue;
    }
    auto [minLabel, maxLabel] = std::minmax_element(
        measurement.localLabels.begin(), measurement.localLabels.end());
    firstLocal = std::min(firstLocal, *minLabel);
    lastLocal = std::max(lastLocal, *maxLabel);
    seenGlobalLabels.insert(measurement.globalLabels.begin(),
                            measurement.globalLabels.end());
  }
  if (lastLocal < firstLocal) {
    return Mille::MilleDecoder::ReadResult::error;
  }
  targetState.trackParametersDim = lastLocal - firstLocal + 1;
  targetState.alignmentDof = seenGlobalLabels.size();

  for (int label : seenGlobalLabels) {
    if ((label - 1) % Acts::eAlignmentSize == 0) {
      seenSurfaceLabels.insert((label - 1) / Acts::eAlignmentSize);
    }
  };

  /// the trackAlignmentState uses an internal indexing for alignment
  /// parameters - so remap the indices to replicate this internal logic.
  std::map<int, int> globalToInternal;

  /// try to map indices from the global indexed surface list, if we have it.
  unsigned long iExtra = 0;
  for (auto [surface, index] : idxedAlignSurfaces) {
    if (seenSurfaceLabels.contains(index)) {
      targetState.alignedSurfaces.emplace(surface,
                                          std::make_pair(index, iExtra++));
    }
  }

  // now use this to fill our internal "global to local" mapping function
  for (auto [surface, indices] : targetState.alignedSurfaces) {
    auto [globIx, intIx] = indices;
    for (std::size_t iAli = 0; iAli < Acts::eAlignmentSize; ++iAli) {
      globalToInternal.emplace(globalIndexSurfToParam(globIx, iAli),
                               internalIndexSurfToParam(intIx, iAli));
    }
  }

  // Now we have the needed information to initialise our matrices
  targetState.measurementCovariance = Acts::DynamicMatrix::Zero(
      targetState.measurementDim, targetState.measurementDim);

  targetState.projectionMatrix = Acts::DynamicMatrix::Zero(
      targetState.measurementDim, targetState.trackParametersDim);

  targetState.alignmentToResidualDerivative = Acts::DynamicMatrix::Zero(
      targetState.measurementDim, targetState.alignmentDof);

  targetState.trackParametersCovariance = Acts::DynamicMatrix::Zero(
      targetState.trackParametersDim, targetState.trackParametersDim);

  targetState.residual = Acts::DynamicVector::Zero(targetState.measurementDim);

  /// Second loop - fill the matrices

  std::size_t iMeas = 0;
  for (const auto& measurement : measurements) {
    // need distinction: Measurement on surface vs. correlation term.
    // The reason is that ACTS only counts surface measurements, and stores
    // the correlation information directly in the track parameter covariance.
    // MillePede considers constraints to be additional measurements.
    // A MillePede pseudomeasurement has residual 0 and no global derivatives.
    bool isMeasurementOnSurface = (measurement.measurement != 0 ||
                                   !measurement.globalDerivatives.empty());
    // surface measurements populate the residual vector and measurement
    // covariance matrix
    if (isMeasurementOnSurface) {
      targetState.residual(iMeas) = measurement.measurement;
      targetState.measurementCovariance(iMeas, iMeas) =
          measurement.uncertainty * measurement.uncertainty;
    }
    // loop over all track parameters affecting this measurement
    for (std::size_t iLoc = 0; iLoc < measurement.localLabels.size(); ++iLoc) {
      // find out where to book it in the ACTS matrix
      unsigned int localIndex = measurement.localLabels[iLoc] - firstLocal;
      // if we are a surface measurement, fill the projection matrix
      if (isMeasurementOnSurface) {
        targetState.projectionMatrix(iMeas, localIndex) =
            measurement.localDerivatives[iLoc];
      }
      // now fill the covariance matrix by looping over all products of (local)
      // derivatives
      for (std::size_t jLoc = 0; jLoc < measurement.localLabels.size();
           ++jLoc) {
        // again determine where to book the column index
        unsigned int localIndex2 = measurement.localLabels[jLoc] - firstLocal;
        // and update the covariance.
        targetState.trackParametersCovariance(localIndex, localIndex2) +=
            measurement.localDerivatives[iLoc] *
            measurement.localDerivatives[jLoc] / measurement.uncertainty /
            measurement.uncertainty;
      }
    }
    // loop over global (= alignment) parameters affecting the measurement
    for (std::size_t iGlob = 0; iGlob < measurement.globalLabels.size();
         ++iGlob) {
      // find out where to book - here we need to map to the ACTS track-level
      // indexing scheme
      int internalAliIndex = globalToInternal[measurement.globalLabels[iGlob]];
      // and update the alignment-to-residual derivative matrix.
      targetState.alignmentToResidualDerivative(iMeas, internalAliIndex) =
          measurement.globalDerivatives[iGlob];
    }
    // increment the measurement-on-surface index every time we finish
    // processing one.
    if (isMeasurementOnSurface)
      if (isMeasurementOnSurface) {
        ++iMeas;
      }
  }

  /// (carefully) invert the covariance - upstairs, we filled it as a weight
  /// matrix
  auto solver = targetState.trackParametersCovariance.ldlt();
  targetState.trackParametersCovariance =
      solver.solve(Acts::DynamicMatrix::Identity(
          targetState.trackParametersDim, targetState.trackParametersDim));

  /// and calculate the dependent members
  /// (first and second derivatives, chi2) of the state.
  ActsAlignment::detail::finaliseTrackAlignState(targetState);

  return res;
}

void dumpAsMillepedeRes(const ActsAlignment::AlignmentResult& result,
                        std::ostream& out) {
  for (auto [surface, index] : result.idxedAlignSurfaces) {
    for (std::size_t i = 0; i < Acts::eAlignmentSize; ++i) {
      std::size_t row = Acts::eAlignmentSize * index + i;
      out << std::setw(8)
          << row + 1
          // column 2: parameter delta from fit
          << "   " << std::setw(12)
          << result.deltaAlignmentParameters(row)
          // column 3: optional pre-sigma for parameter delta
          << "   " << std::setw(4)
          << 0.
          // column 4: change of parameter delta w.r.t start value
          << "   " << std::setw(12)
          << result.deltaAlignmentParameters(row)
          // column 5: uncertainty of parameter delta from the fit
          << "   " << std::setw(12)
          << std::sqrt(result.alignmentCovariance(row, row))
          // column 6: count of measurements seeing this parameter.
          // Not available in ACTS solver, write a placeholder of 99
          << "    99 " << std::endl;
    }
  }
}

}  // namespace ActsPlugins::ActsToMille
