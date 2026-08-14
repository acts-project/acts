// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/SympyCovarianceEngine.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Propagator/detail/JacobianEngine.hpp"
#include "Acts/Propagator/detail/SympyJacobianEngine.hpp"

#include "codegen/sympy_cov_math.hpp"

namespace Acts::detail {

namespace {

/// Transport a bound covariance, picking the cheaper kernel where it applies
///
/// d(q/p)/d(q/p) is exactly one unless the track crossed material: the vacuum
/// step kernel never writes that entry of the bound-to-free jacobian, so it
/// keeps the value initialisation gave it, bit for bit. The comparison is
/// therefore structural rather than numerical, and the vacuum path skips the
/// multiplications the dense one needs.
///
/// @param jacobian the full bound-to-bound transport jacobian
/// @param in the covariance to transport
/// @param [out] out the transported covariance
void applyBoundCovarianceTransport(const BoundMatrix& jacobian,
                                   const BoundMatrix& in, BoundMatrix& out) {
  if (jacobian(eBoundQOverP, eBoundQOverP) == 1) {
    transportCovarianceToBoundImpl(in.data(), jacobian.data(), out.data());
  } else {
    transportCovarianceToBoundDenseImpl(in.data(), jacobian.data(), out.data());
  }
}

}  // namespace

/// Some type defs
using Jacobian = BoundMatrix;
using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;

Result<BoundState> sympy::boundState(
    const GeometryContext& geoContext, const Surface& surface,
    BoundMatrix& boundCovariance, BoundMatrix& fullTransportJacobian,
    FreeMatrix& freeTransportJacobian, FreeVector& freeToPathDerivatives,
    BoundToFreeMatrix& boundToFreeJacobian,
    const std::optional<FreeMatrix>& additionalFreeCovariance,
    FreeVector& freeParameters, const ParticleHypothesis& particleHypothesis,
    bool covTransport, double accumulatedPath,
    const FreeToBoundCorrection& freeToBoundCorrection) {
  // Create the bound parameters
  Result<BoundVector> bv =
      transformFreeToBoundParameters(freeParameters, surface, geoContext);
  if (!bv.ok()) {
    return bv.error();
  }

  // Covariance transport
  std::optional<BoundMatrix> cov = std::nullopt;
  if (covTransport) {
    // Calculate the jacobian and transport the covarianceMatrix to final local.
    // Then reinitialize the transportJacobian, derivatives and the
    // boundToFreeJacobian
    transportCovarianceToBound(
        geoContext, surface, boundCovariance, fullTransportJacobian,
        freeTransportJacobian, freeToPathDerivatives, boundToFreeJacobian,
        additionalFreeCovariance, freeParameters, freeToBoundCorrection);
    cov = boundCovariance;
  }

  // Create the bound state
  return std::make_tuple(
      BoundTrackParameters(surface.getSharedPtr(), *bv, std::move(cov),
                           particleHypothesis),
      fullTransportJacobian, accumulatedPath);
}

BoundState sympy::curvilinearState(
    BoundMatrix& boundCovariance, BoundMatrix& fullTransportJacobian,
    FreeMatrix& freeTransportJacobian, FreeVector& freeToPathDerivatives,
    BoundToFreeMatrix& boundToFreeJacobian,
    const std::optional<FreeMatrix>& additionalFreeCovariance,
    const FreeVector& freeParameters,
    const ParticleHypothesis& particleHypothesis, bool covTransport,
    double accumulatedPath) {
  const Vector3& direction = freeParameters.segment<3>(eFreeDir0);

  // Covariance transport
  std::optional<BoundMatrix> cov = std::nullopt;
  if (covTransport) {
    // Calculate the jacobian and transport the covarianceMatrix to final local.
    // Then reinitialize the transportJacobian, derivatives and the
    // boundToFreeJacobian
    transportCovarianceToCurvilinear(boundCovariance, fullTransportJacobian,
                                     freeTransportJacobian,
                                     freeToPathDerivatives, boundToFreeJacobian,
                                     additionalFreeCovariance, direction);
    cov = boundCovariance;
  }

  // Create the curvilinear parameters
  Vector4 pos4 = Vector4::Zero();
  pos4[ePos0] = freeParameters[eFreePos0];
  pos4[ePos1] = freeParameters[eFreePos1];
  pos4[ePos2] = freeParameters[eFreePos2];
  pos4[eTime] = freeParameters[eFreeTime];
  BoundTrackParameters curvilinearParams =
      BoundTrackParameters::createCurvilinear(
          pos4, direction, freeParameters[eFreeQOverP], std::move(cov),
          particleHypothesis);
  // Create the curvilinear state
  return {std::move(curvilinearParams), fullTransportJacobian, accumulatedPath};
}

void sympy::transportCovarianceToBound(
    const GeometryContext& geoContext, const Surface& surface,
    BoundMatrix& boundCovariance, BoundMatrix& fullTransportJacobian,
    FreeMatrix& freeTransportJacobian, FreeVector& freeToPathDerivatives,
    BoundToFreeMatrix& boundToFreeJacobian,
    const std::optional<FreeMatrix>& additionalFreeCovariance,
    FreeVector& freeParameters,
    const FreeToBoundCorrection& freeToBoundCorrection) {
  FreeToBoundMatrix freeToBoundJacobian;

  // Calculate the full jacobian from local parameters at the start surface to
  // current bound parameters
  sympy::boundToBoundTransportJacobian(
      geoContext, surface, freeParameters, boundToFreeJacobian,
      freeTransportJacobian, freeToBoundJacobian, freeToPathDerivatives,
      fullTransportJacobian);

  bool correction = false;
  if (freeToBoundCorrection) {
    BoundToFreeMatrix startBoundToFinalFreeJacobian =
        freeTransportJacobian * boundToFreeJacobian;
    FreeMatrix freeCovariance = startBoundToFinalFreeJacobian *
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
      freeParameters =
          transformBoundToFreeParameters(surface, geoContext, boundParams);

      // 2. Update the bound covariance
      boundCovariance = std::get<BoundMatrix>(correctedValue);

      correction = true;
    }
  }

  if (!correction) {
    // Apply the actual covariance transport to get covariance of the current
    // bound parameters
    BoundMatrix newBoundCovariance;
    applyBoundCovarianceTransport(fullTransportJacobian, boundCovariance,
                                  newBoundCovariance);
    boundCovariance = newBoundCovariance;
  }

  if (additionalFreeCovariance) {
    boundCovariance += freeToBoundJacobian * (*additionalFreeCovariance) *
                       freeToBoundJacobian.transpose();
  }

  // Reinitialize jacobian components:
  // ->The transportJacobian is reinitialized to Identity
  // ->The derivatives is reinitialized to Zero
  // ->The boundToFreeJacobian is initialized to that at the current surface
  reinitializeJacobians(geoContext, surface, freeTransportJacobian,
                        freeToPathDerivatives, boundToFreeJacobian,
                        freeParameters);
}

void sympy::transportCovarianceToCurvilinear(
    BoundMatrix& boundCovariance, BoundMatrix& fullTransportJacobian,
    FreeMatrix& freeTransportJacobian, FreeVector& freeToPathDerivatives,
    BoundToFreeMatrix& boundToFreeJacobian,
    const std::optional<FreeMatrix>& additionalFreeCovariance,
    const Vector3& direction) {
  FreeToBoundMatrix freeToBoundJacobian;

  // Calculate the full jacobian from local parameters at the start surface to
  // current curvilinear parameters
  sympy::boundToCurvilinearTransportJacobian(
      direction, boundToFreeJacobian, freeTransportJacobian,
      freeToBoundJacobian, freeToPathDerivatives, fullTransportJacobian);

  // Apply the actual covariance transport to get covariance of the current
  // curvilinear parameters
  BoundMatrix newBoundCovariance;
  applyBoundCovarianceTransport(fullTransportJacobian, boundCovariance,
                                newBoundCovariance);
  boundCovariance = newBoundCovariance;

  if (additionalFreeCovariance) {
    boundCovariance += freeToBoundJacobian * (*additionalFreeCovariance) *
                       freeToBoundJacobian.transpose();
  }

  // Reinitialize jacobian components:
  // ->The free transportJacobian is reinitialized to Identity
  // ->The path derivatives is reinitialized to Zero
  // ->The boundToFreeJacobian is reinitialized to that at the current
  // curvilinear surface
  reinitializeJacobians(freeTransportJacobian, freeToPathDerivatives,
                        boundToFreeJacobian, direction);
}

}  // namespace Acts::detail
