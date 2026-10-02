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
#include <atomic>
#include <cstddef>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>

#include <Eigen/src/Core/Matrix.h>
#include <Mille/MilleDataStructures.h>
#include <Mille/MilleRecord.h>

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
                 MilleRecord& record, bool removeUnconstrainedTrackPar,
                 const Acts::Logger& logger) {
  // spawn a local buffer to be able to assemble the record without lock
  // contention.
  std::unique_ptr<Mille::MilleRecord> milleLocalBuf = record.spawnLocalBuffer();

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
    // collect (sorted by covariance) the indices of all track parameters.
    // A multimap: parameters with equal variances must all be considered.
    std::multimap<double, std::size_t> trkParByCov;
    for (std::size_t k = 0; k < state.trackParametersDim; ++k) {
      trkParByCov.emplace(state.trackParametersCovariance(k, k), k);
    }

    // now, loop through the parameter list and look for huge jumps.
    std::size_t nKeptMeasured = 0;
    double prev = 0;
    for (const auto& [sigmaSquared, index] : trkParByCov) {
      // a jump of 1e6 is indicative that we are not in Kansas anymore
      if (prev != 0 && sigmaSquared > 1e6 * prev) {
        // A parameter that a measurement projects onto is constrained by it
        // (its variance is bounded by the measurement variance) and must stay
        // a local parameter: without it, the measurement would be written
        // with its alignment derivatives but without its dependence on the
        // track. The precision is deliberately 0 (exact): H only selects the
        // measured parameters (entries 0 or 1), and Mille drops only
        // exactly-zero local derivatives.
        if (!state.projectionMatrix.col(index).isZero(0.)) {
          // keep it, and leave prev unchanged: the decisions on the other
          // parameters are the same as without this check
          ++nKeptMeasured;
          continue;
        }
        skippedTrackParams.insert(index);
      } else {
        prev = sigmaSquared;
      }
    }
    if (nKeptMeasured > 0) {
      // warn once, the criterion can misfire on every track of a detector
      static std::atomic<bool> warned{false};
      if (logger.doPrint(Acts::Logging::WARNING) && !warned.exchange(true)) {
        ACTS_WARNING("removeUnconstrainedTrackPar: kept "
                     << nKeptMeasured
                     << " measured track parameter(s) of a track that the "
                        "variance-jump criterion flagged as unconstrained. The "
                        "criterion compares variances of different units; "
                        "consider removeUnconstrainedTrackPar = false. "
                        "Further occurrences are reported at DEBUG level.");
      } else {
        ACTS_DEBUG("removeUnconstrainedTrackPar: kept "
                   << nKeptMeasured << " measured track parameter(s)");
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
    milleLocalBuf->addData(
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

  // Gradient of the measurement chi2 at the smoothed track, -2 * g. The
  // smoothed track is the minimum of the full chi2, so the correlation term
  // has to cancel it: its pseudo-measurements carry the residuals rho with
  // sum_i lambda_i rho_i v_i = -g. With zero residuals, pede's local fit would
  // move off the smoothed track by C * g, and the alignment gradient would
  // pick up a spurious projector -2 A^T V^-1 H C g. Whenever the smoothed
  // track has kinks (multiple scattering absorbing a misalignment), this biases
  // the result towards zero.
  const Acts::DynamicVector measGradient =
      updatedProjection.transpose() * state.measurementCovariance.inverse() *
      state.residual;

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
    const double pseudoResidual =
        -eigenVecs.col(iMeas).dot(measGradient) / eigenVals(iMeas);
    // and write a pseudo-measurement to Mille.
    milleLocalBuf->addData(
        // residual keeping the smoothed track at the local-fit minimum
        pseudoResidual,
        // EV == weight = 1/sigma^2
        1. / std::sqrt(eigenVals(iMeas)),
        // local parameter indices
        localIndices,
        // local derivatives
        localDeriv, globalIndices, globalDeriv);
  }
  // track is fully written - end the record in Mille
  // NB: This will automatically propagate the local buffer content to
  // the parent instance passed by the caller.
  milleLocalBuf->writeRecord();
}

Mille::MilleDecoder::ReadResult unpackMilleRecord(
    Mille::IMilleReader& reader,
    ActsAlignment::detail::TrackAlignmentState& targetState,
    const std::unordered_map<const Acts::Surface*, std::size_t>&
        idxedAlignSurfaces,
    const Acts::Logger& logger) {
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

  // Still here - we got a valid record! Build it in a fresh state and only
  // hand it over at the end: nothing may carry over from a previous record
  // (callers reuse the target state across records), and on a read error the
  // target state stays untouched.
  // In the following, we emulate what MillePede-II is doing internally.
  // This is somewhat approximate, as we do not run any of the cleaning /
  // conditioning performed by MP-II.
  ActsAlignment::detail::TrackAlignmentState state;

  // Step 1: Parameter discovery
  // Goal: Identify all existing parameters and assign internal indices.
  int firstLocal = 99999;
  int lastLocal = 0;
  std::set<int> seenGlobalLabels;
  std::set<int> seenSurfaceLabels;
  // Every entry of the record is a measurement of the local fit, as in
  // MillePede: the surface measurements, and the pseudo-measurements that
  // encode the track model (for dumpToMille, the Kalman correlations). There
  // is no need to tell them apart, which the record could not do reliably
  // anyway. A pseudo-measurement has no alignment derivatives, and its
  // residual enters the right-hand side: the chi2 derivatives below are those
  // of the local fit chi2, minimised over the track parameters, for any track
  // model written to the record.
  state.measurementDim = measurements.size();

  // discover labels in use
  for (const Mille::MilleMeasurement& measurement : measurements) {
    if (measurement.localLabels.empty()) {
      // A measurement always depends on the track. An entry with alignment
      // derivatives but no local derivatives has lost its track parameters
      // when it was written, and cannot be placed in the track model.
      if (!measurement.globalLabels.empty()) {
        ACTS_ERROR(
            "Mille record with an entry that has alignment derivatives "
            "(first label "
            << measurement.globalLabels.front()
            << ") but no local derivatives: its track parameters were "
               "dropped when it was written. Skipping the record.");
        return Mille::MilleDecoder::ReadResult::error;
      }
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
    ACTS_ERROR("Mille record without local derivatives. Skipping the record.");
    return Mille::MilleDecoder::ReadResult::error;
  }
  state.trackParametersDim = lastLocal - firstLocal + 1;

  // A surface is on the track if any of its labels appears: Mille does not
  // write zero derivatives, so any single label (e.g. the first one) may be
  // missing from a record.
  for (int label : seenGlobalLabels) {
    seenSurfaceLabels.insert((label - 1) / Acts::eAlignmentSize);
  }

  /// the trackAlignmentState uses an internal indexing for alignment
  /// parameters - so remap the indices to replicate this internal logic.
  /// Surfaces are numbered in the order of their global labels, and the
  /// derivative matrix gets all parameters of every surface on the track,
  /// as in the TrackAlignmentState that the record was written from.
  std::map<int, int> globalToInternal;
  std::map<std::size_t, std::size_t> surfaceToInternal;
  for (int surfaceLabel : seenSurfaceLabels) {
    const std::size_t intIx = surfaceToInternal.size();
    surfaceToInternal.emplace(surfaceLabel, intIx);
    for (std::size_t iAli = 0; iAli < Acts::eAlignmentSize; ++iAli) {
      globalToInternal.emplace(globalIndexSurfToParam(surfaceLabel, iAli),
                               internalIndexSurfToParam(intIx, iAli));
    }
  }
  state.alignmentDof = Acts::eAlignmentSize * surfaceToInternal.size();

  /// link the surfaces back to the geometry, if we have the indexed list
  for (auto [surface, index] : idxedAlignSurfaces) {
    if (auto it = surfaceToInternal.find(index);
        it != surfaceToInternal.end()) {
      state.alignedSurfaces.emplace(surface, std::make_pair(index, it->second));
    }
  }
  // a label of a surface that the indexed list does not know
  if (!idxedAlignSurfaces.empty() &&
      state.alignedSurfaces.size() != surfaceToInternal.size()) {
    ACTS_ERROR("Mille record with alignment labels of "
               << surfaceToInternal.size() - state.alignedSurfaces.size()
               << " surface(s) missing from the indexed alignable surfaces. "
                  "Skipping the record.");
    return Mille::MilleDecoder::ReadResult::error;
  }

  // Now we have the needed information to initialise our matrices
  state.measurementCovariance =
      Acts::DynamicMatrix::Zero(state.measurementDim, state.measurementDim);

  state.projectionMatrix =
      Acts::DynamicMatrix::Zero(state.measurementDim, state.trackParametersDim);

  state.alignmentToResidualDerivative =
      Acts::DynamicMatrix::Zero(state.measurementDim, state.alignmentDof);

  state.trackParametersCovariance = Acts::DynamicMatrix::Zero(
      state.trackParametersDim, state.trackParametersDim);

  state.residual = Acts::DynamicVector::Zero(state.measurementDim);

  /// Second loop - fill the matrices

  for (std::size_t iMeas = 0; iMeas < measurements.size(); ++iMeas) {
    const Mille::MilleMeasurement& measurement = measurements[iMeas];
    // every entry populates the residual vector and measurement covariance
    // matrix
    state.residual(iMeas) = measurement.measurement;
    state.measurementCovariance(iMeas, iMeas) =
        measurement.uncertainty * measurement.uncertainty;
    // loop over all track parameters affecting this measurement
    for (std::size_t iLoc = 0; iLoc < measurement.localLabels.size(); ++iLoc) {
      // find out where to book it in the ACTS matrix
      unsigned int localIndex = measurement.localLabels[iLoc] - firstLocal;
      // fill the projection matrix
      state.projectionMatrix(iMeas, localIndex) =
          measurement.localDerivatives[iLoc];
      // now fill the covariance matrix by looping over all products of (local)
      // derivatives
      for (std::size_t jLoc = 0; jLoc < measurement.localLabels.size();
           ++jLoc) {
        // again determine where to book the column index
        unsigned int localIndex2 = measurement.localLabels[jLoc] - firstLocal;
        // and update the covariance.
        state.trackParametersCovariance(localIndex, localIndex2) +=
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
      // (every label is mapped: its surface was registered above)
      const int internalAliIndex =
          globalToInternal.at(measurement.globalLabels[iGlob]);
      // and update the alignment-to-residual derivative matrix.
      state.alignmentToResidualDerivative(iMeas, internalAliIndex) =
          measurement.globalDerivatives[iGlob];
    }
  }

  /// (carefully) invert the covariance - upstairs, we filled it as a weight
  /// matrix
  auto solver = state.trackParametersCovariance.ldlt();
  state.trackParametersCovariance = solver.solve(Acts::DynamicMatrix::Identity(
      state.trackParametersDim, state.trackParametersDim));

  /// and calculate the dependent members
  /// (first and second derivatives, chi2) of the state.
  ActsAlignment::detail::finaliseTrackAlignState(state);

  targetState = std::move(state);
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
