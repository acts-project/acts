// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/CovarianceEngine.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/GenericCurvilinearTrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Propagator/detail/JacobianEngine.hpp"
#include "Acts/Utilities/Result.hpp"

#include <optional>
#include <system_error>
#include <utility>

namespace Acts {

/// Some type defs
using Jacobian = BoundMatrix;
using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;
using CurvilinearState =
    std::tuple<CurvilinearTrackParameters, Jacobian, double>;

Result<BoundState> detail::boundState(
    const GeometryContext& geoContext, BoundSquareMatrix& covarianceMatrix,
    BoundMatrix& jacobian, FreeMatrix& transportJacobian,
    FreeVector& derivatives, BoundToFreeMatrix& jacToGlobal,
    FreeVector& parameters, const ParticleHypothesis& particleHypothesis,
    bool covTransport, double accumulatedPath, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) {
  // Covariance transport
  std::optional<BoundSquareMatrix> cov = std::nullopt;
  if (covTransport) {
    // Initialize the jacobian from start local to final local
    jacobian = BoundMatrix::Identity();
    // Calculate the jacobian and transport the covarianceMatrix to final local.
    // Then reinitialize the transportJacobian, derivatives and the
    // jacToGlobal
    transportCovarianceToBound(geoContext, covarianceMatrix, jacobian,
                               transportJacobian, derivatives, jacToGlobal,
                               parameters, surface, freeToBoundCorrection);
  }
  if (covarianceMatrix != BoundSquareMatrix::Zero()) {
    cov = covarianceMatrix;
  }

  // Create the bound parameters
  Result<BoundVector> bv =
      detail::transformFreeToBoundParameters(parameters, surface, geoContext);
  if (!bv.ok()) {
    return bv.error();
  }
  // Create the bound state
  return std::make_tuple(
      BoundTrackParameters(surface.getSharedPtr(), *bv, std::move(cov),
                           particleHypothesis),
      jacobian, accumulatedPath);
}

CurvilinearState detail::curvilinearState(
    BoundSquareMatrix& covarianceMatrix, BoundMatrix& jacobian,
    FreeMatrix& transportJacobian, FreeVector& derivatives,
    BoundToFreeMatrix& jacToGlobal, const FreeVector& parameters,
    const ParticleHypothesis& particleHypothesis, bool covTransport,
    double accumulatedPath) {
  const Vector3& direction = parameters.segment<3>(eFreeDir0);

  // Covariance transport
  std::optional<BoundSquareMatrix> cov = std::nullopt;
  if (covTransport) {
    // Initialize the jacobian from start local to final local
    jacobian = BoundMatrix::Identity();
    // Calculate the jacobian and transport the covarianceMatrix to final local.
    // Then reinitialize the transportJacobian, derivatives and the
    // jacToGlobal
    transportCovarianceToCurvilinear(covarianceMatrix, jacobian,
                                     transportJacobian, derivatives,
                                     jacToGlobal, direction);
  }
  if (covarianceMatrix != BoundSquareMatrix::Zero()) {
    cov = covarianceMatrix;
  }

  // Create the curvilinear parameters
  Vector4 pos4 = Vector4::Zero();
  pos4[ePos0] = parameters[eFreePos0];
  pos4[ePos1] = parameters[eFreePos1];
  pos4[ePos2] = parameters[eFreePos2];
  pos4[eTime] = parameters[eFreeTime];
  CurvilinearTrackParameters curvilinearParams(
      pos4, direction, parameters[eFreeQOverP], std::move(cov),
      particleHypothesis);
  // Create the curvilinear state
  return std::make_tuple(std::move(curvilinearParams), jacobian,
                         accumulatedPath);
}

void detail::transportCovarianceToBound(
    const GeometryContext& geoContext, BoundSquareMatrix& boundCovariance,
    BoundMatrix& fullTransportJacobian, FreeMatrix& freeTransportJacobian,
    FreeVector& freeToPathDerivatives, BoundToFreeMatrix& boundToFreeJacobian,
    FreeVector& freeParameters, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) {
  // Calculate the full jacobian from local parameters at the start surface to
  // current bound parameters
  fullTransportJacobian = boundToBoundTransportJacobian(
      geoContext, freeParameters, boundToFreeJacobian, freeTransportJacobian,
      freeToPathDerivatives, surface);

  bool correction = false;
  if (freeToBoundCorrection) {
    BoundToFreeMatrix startBoundToFinalFreeJacobian =
        freeTransportJacobian * boundToFreeJacobian;
    FreeSquareMatrix freeCovariance = startBoundToFinalFreeJacobian *
                                      boundCovariance *
                                      startBoundToFinalFreeJacobian.transpose();

    auto transformer =
        detail::CorrectedFreeToBoundTransformer(freeToBoundCorrection);
    auto correctedRes =
        transformer(freeParameters, freeCovariance, surface, geoContext);

    if (correctedRes.has_value()) {
      auto correctedValue = correctedRes.value();
      BoundVector boundParams = std::get<BoundVector>(correctedValue);
      // 1. Update the free parameters with the corrected bound parameters
      freeParameters = detail::transformBoundToFreeParameters(
          surface, geoContext, boundParams);

      // 2. Update the bound covariance
      boundCovariance = std::get<BoundSquareMatrix>(correctedValue);

      correction = true;
    }
  }

  if (!correction) {
    // Apply the actual covariance transport to get covariance of the current
    // bound parameters
    boundCovariance = fullTransportJacobian * boundCovariance *
                      fullTransportJacobian.transpose();
  }

  // Reinitialize jacobian components:
  // ->The transportJacobian is reinitialized to Identity
  // ->The derivatives is reinitialized to Zero
  // ->The boundToFreeJacobian is initialized to that at the current surface
  reinitializeJacobians(geoContext, freeTransportJacobian,
                        freeToPathDerivatives, boundToFreeJacobian,
                        freeParameters, surface);
}

void detail::transportCovarianceToCurvilinear(
    BoundSquareMatrix& boundCovariance, BoundMatrix& fullTransportJacobian,
    FreeMatrix& freeTransportJacobian, FreeVector& freeToPathDerivatives,
    BoundToFreeMatrix& boundToFreeJacobian, const Vector3& direction) {
  // Calculate the full jacobian from local parameters at the start surface to
  // current curvilinear parameters
  fullTransportJacobian = boundToCurvilinearTransportJacobian(
      direction, boundToFreeJacobian, freeTransportJacobian,
      freeToPathDerivatives);

  // Apply the actual covariance transport to get covariance of the current
  // curvilinear parameters
  boundCovariance = fullTransportJacobian * boundCovariance *
                    fullTransportJacobian.transpose();

  // Reinitialize jacobian components:
  // ->The free transportJacobian is reinitialized to Identity
  // ->The path derivatives is reinitialized to Zero
  // ->The boundToFreeJacobian is reinitialized to that at the current
  // curvilinear surface
  reinitializeJacobians(freeTransportJacobian, freeToPathDerivatives,
                        boundToFreeJacobian, direction);
}

}  // namespace Acts
