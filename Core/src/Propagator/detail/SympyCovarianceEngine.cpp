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

/// Transport a bound covariance, taking the cheaper vacuum kernel where it
/// applies. d(q/p)/d(q/p) is left untouched by a vacuum step, so it is exactly
/// one there and the test is structural rather than numerical.
///
/// @param jacobian the full bound-to-bound transport jacobian
/// @param in the covariance to transport
/// @param [out] out the transported covariance
void applyBoundCovarianceTransport(const BoundMatrix& jacobian,
                                   const BoundMatrix& in, BoundMatrix& out) {
  const auto j = std::span<const double, 36>(jacobian.data(), 36);
  const auto c = std::span<const double, 36>(in.data(), 36);
  const auto o = std::span<double, 36>(out.data(), 36);
  if (jacobian(eBoundQOverP, eBoundQOverP) == 1) {
    transportCovarianceToBoundVacuumImpl(c, j, o);
  } else {
    transportCovarianceToBoundDenseImpl(c, j, o);
  }
}

}  // namespace

Result<void> sympy::transportCovarianceToBound(
    const GeometryContext& geoContext, const Surface& surface,
    BoundMatrix& boundCovariance, BoundMatrix& fullTransportJacobian,
    FreeVector& freeToPathDerivatives, BoundToFreeMatrix& boundToFreeJacobian,
    const std::optional<FreeMatrix>& additionalFreeCovariance,
    FreeVector& freeParameters,
    const FreeToBoundCorrection& freeToBoundCorrection) {
  FreeToBoundMatrix freeToBoundJacobian;

  // Calculate the full jacobian from local parameters at the start surface to
  // current bound parameters
  sympy::boundToBoundTransportJacobian(
      geoContext, surface, freeParameters, boundToFreeJacobian,
      freeToBoundJacobian, freeToPathDerivatives, fullTransportJacobian);

  bool correction = false;
  if (freeToBoundCorrection) {
    FreeMatrix freeCovariance =
        boundToFreeJacobian * boundCovariance * boundToFreeJacobian.transpose();

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
  // ->The derivatives are reinitialized to Zero
  // ->The boundToFreeJacobian is initialized to that at the current surface
  return reinitializeJacobians(geoContext, surface, freeToPathDerivatives,
                               boundToFreeJacobian, freeParameters);
}

void sympy::transportCovarianceToCurvilinear(
    BoundMatrix& boundCovariance, BoundMatrix& fullTransportJacobian,
    FreeVector& freeToPathDerivatives, BoundToFreeMatrix& boundToFreeJacobian,
    const std::optional<FreeMatrix>& additionalFreeCovariance,
    const Vector3& direction) {
  FreeToBoundMatrix freeToBoundJacobian;

  // Calculate the full jacobian from local parameters at the start surface to
  // current curvilinear parameters
  sympy::boundToCurvilinearTransportJacobian(
      direction, boundToFreeJacobian, freeToBoundJacobian,
      freeToPathDerivatives, fullTransportJacobian);

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
  // ->The path derivatives are reinitialized to Zero
  // ->The boundToFreeJacobian is reinitialized to that at the current
  // curvilinear surface
  reinitializeJacobians(freeToPathDerivatives, boundToFreeJacobian, direction);
}

}  // namespace Acts::detail
