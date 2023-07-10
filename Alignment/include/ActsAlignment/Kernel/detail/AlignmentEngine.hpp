// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//Alignment/include/ActsAlignment/Kernel/detail
#pragma once

#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"
#include "TH1F.h"
#include "TCanvas.h"

#include <unordered_map>

namespace ActsAlignment {
namespace detail {

using namespace Acts;

struct TrackAlignmentState {
  size_t measurementDim = 0;
  size_t trackParametersDim = 0;
  size_t alignmentDof = 0;
  ActsDynamicMatrix measurementCovariance;
  ActsDynamicMatrix trackParametersCovariance;
  ActsDynamicMatrix projectionMatrix;
  ActsDynamicVector residual;
  ActsDynamicMatrix residualCovariance;
  double chi2 = 0;
  ActsDynamicMatrix alignmentToResidualDerivative;
  ActsDynamicVector alignmentToChi2Derivative;
  ActsDynamicMatrix alignmentToChi2SecondDerivative;
  std::unordered_map<const Surface*, std::pair<size_t, size_t>> alignedSurfaces;
  ActsDynamicVector misalignmentParameters;
  ActsDynamicMatrix misalignmentCovariance;

  TrackAlignmentState()
      : misalignmentParameters(ActsDynamicVector::Zero(alignmentDof)) {}
};

void resetAlignmentDerivative(Acts::AlignmentToBoundMatrix& alignToBound,
                              AlignmentMask mask,
                              const ActsDynamicVector& misalignmentParameters);

template <typename traj_t, typename parameters_t = BoundTrackParameters>
TrackAlignmentState trackAlignmentState(
    const GeometryContext& gctx, const traj_t& multiTraj,
    const size_t& entryIndex,
    const std::pair<ActsDynamicMatrix, std::unordered_map<size_t, size_t>>&
        globalTrackParamsCov,
    const std::unordered_map<const Surface*, size_t>& idxedAlignSurfaces,
    const AlignmentMask& alignMask,
    const ActsDynamicMatrix& misalignmentCovariance) {
  using CovMatrix = typename parameters_t::CovarianceMatrix;

  TrackAlignmentState alignState;

  std::vector<std::pair<size_t, bool>> measurementStates;
  size_t nAlignSurfaces = 0;

  multiTraj.visitBackwards(entryIndex, [&](const auto& ts) {
    if (!ts.hasSmoothed()) {
      return true;
    }

    if (!ts.typeFlags().test(TrackStateFlag::MeasurementFlag)) {
      return true;
    }

    bool isAlignable = false;
    const auto surface = &ts.referenceSurface();
    auto it = idxedAlignSurfaces.find(surface);
    if (it != idxedAlignSurfaces.end()) {
      isAlignable = true;
      alignState.alignedSurfaces[surface].first = it->second;
      nAlignSurfaces++;
    }

    measurementStates.emplace_back(ts.index(), isAlignable);
    alignState.measurementDim += ts.calibratedSize();
    return true;
  });

  if (nAlignSurfaces == 0) {
    return alignState;
  }

  alignState.alignmentDof = eAlignmentSize * nAlignSurfaces;
  alignState.trackParametersDim = eBoundSize * measurementStates.size();

  alignState.measurementCovariance = ActsDynamicMatrix::Zero(
      alignState.measurementDim, alignState.measurementDim);
  alignState.projectionMatrix = ActsDynamicMatrix::Zero(
      alignState.measurementDim, alignState.trackParametersDim);
  alignState.alignmentToResidualDerivative = ActsDynamicMatrix::Zero(
      alignState.measurementDim, alignState.alignmentDof);
  alignState.trackParametersCovariance = ActsDynamicMatrix::Zero(
      alignState.trackParametersDim, alignState.trackParametersDim);
  alignState.residual = ActsDynamicVector::Zero(alignState.measurementDim);

  const auto& [sourceTrackParamsCov, stateRowIndices] = globalTrackParamsCov;

  size_t iMeasurement = alignState.measurementDim;
  size_t iParams = alignState.trackParametersDim;
  size_t iSurface = nAlignSurfaces;
  for (const auto& [rowStateIndex, isAlignable] : measurementStates) {
    const auto& state = multiTraj.getTrackState(rowStateIndex);
    const size_t measdim = state.calibratedSize();
    iMeasurement -= measdim;
    iParams -= eBoundSize;

    const ActsDynamicMatrix measCovariance =
        state.effectiveCalibratedCovariance();
    alignState.measurementCovariance.block(iMeasurement, iMeasurement, measdim,
                                           measdim) = measCovariance;

    const ActsDynamicMatrix H = state.effectiveProjector();
    alignState.projectionMatrix.block(iMeasurement, iParams, measdim,
                                      eBoundSize) = H;

    alignState.residual.segment(iMeasurement, measdim) =
        state.effectiveCalibrated() - H * state.smoothed();

    if (isAlignable) {
      iSurface -= 1;
      const auto surface = &state.referenceSurface();
      alignState.alignedSurfaces.at(surface).second = iSurface;

      const FreeVector freeParams =
          Acts::MultiTrajectoryHelpers::freeSmoothed(gctx, state);
      const Vector3 direction = freeParams.segment<3>(eFreeDir0);
      FreeVector pathDerivative = FreeVector::Zero();
      pathDerivative.head<3>() = direction;

      AlignmentToBoundMatrix alignToBound =
          surface->alignmentToBoundDerivative(gctx, freeParams, pathDerivative);

      resetAlignmentDerivative(alignToBound, alignMask,
                               alignState.misalignmentParameters);

      alignState.alignmentToResidualDerivative.block(
          iMeasurement, iSurface * eAlignmentSize, measdim,
          eAlignmentSize) = -H * alignToBound;
    }

    for (unsigned int iColState = 0; iColState < measurementStates.size();
         iColState++) {
      size_t colStateIndex = measurementStates.at(iColState).first;
      CovMatrix correlation =
          sourceTrackParamsCov.block<eBoundSize, eBoundSize>(
              stateRowIndices.at(rowStateIndex),
              stateRowIndices.at(colStateIndex));

      size_t iCol =
          alignState.trackParametersDim - (iColState + 1) * eBoundSize;
      alignState.trackParametersCovariance.block<eBoundSize, eBoundSize>(
          iParams, iCol) = correlation;
    }
  }

  alignState.chi2 = alignState.residual.transpose() *
                    alignState.measurementCovariance.inverse() *
                    alignState.residual;
  alignState.alignmentToChi2Derivative =
      2 * alignState.alignmentToResidualDerivative.transpose() *
      alignState.measurementCovariance.inverse() * alignState.residual;
  alignState.alignmentToChi2SecondDerivative =
      2 * alignState.alignmentToResidualDerivative.transpose() *
      alignState.measurementCovariance.inverse() *
      alignState.alignmentToResidualDerivative;

  alignState.residualCovariance = alignState.measurementCovariance -
                                  alignState.projectionMatrix *
                                      alignState.trackParametersCovariance *
                                      alignState.projectionMatrix.transpose();

  alignState.alignmentToChi2SecondDerivative.block(
      0, 0, alignState.alignmentDof, alignState.alignmentDof) +=
      alignState.misalignmentCovariance;

  return alignState;
}

void resetAlignmentDerivative(Acts::AlignmentToBoundMatrix& alignToBound,
                              AlignmentMask mask,
                              const ActsDynamicVector& misalignmentParameters) {
  const size_t alignmentDof = misalignmentParameters.size();

  for (size_t i = 0; i < alignmentDof; ++i) {
    if (mask.isFixed(i)) {
      alignToBound.col(i).setZero();
    }
  }
}

void plotUnbiasedResiduals(const TrackAlignmentState& alignState, TH1F* histogram) {
  ActsDynamicVector unbiasedResiduals = alignState.residual -
      alignState.projectionMatrix * alignState.misalignmentParameters;
  for (size_t i = 0; i < alignState.measurementDim; ++i) {
    histogram->Fill(unbiasedResiduals(i));
  }
}

