// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/CovarianceEngine.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/GenericCurvilinearTrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Propagator/detail/JacobianEngine.hpp"
#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Utilities/JacobianHelpers.hpp"
#include "Acts/Utilities/Result.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <system_error>
#include <type_traits>
#include <utility>

namespace Acts {
namespace {
/// Some type defs
using Jacobian = BoundMatrix;

using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;
using CurvilinearState =
    std::tuple<CurvilinearTrackParameters, Jacobian, double>;

/// @brief Evaluate the projection Jacobian from free to curvilinear parameters
/// without transport jacobian.
///
/// @param [in] direction Normalised direction vector
///
/// @return Projection Jacobian
FreeToBoundMatrix freeToCurvilinearJacobian(const Vector3& direction) {
  // Optimized trigonometry on the propagation direction
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)
  // prepare the jacobian to curvilinear
  FreeToBoundMatrix jacToCurv = FreeToBoundMatrix::Zero();
  if (std::abs(z) < s_curvilinearProjTolerance) {
    auto [cosPhi, sinPhi, cosTheta, sinTheta] =
        VectorHelpers::evaluateTrigonomics(direction);
    // We normally operate in curvilinear coordinates defined as follows
    jacToCurv(eBoundLoc0, eFreePos0) = -sinPhi;
    jacToCurv(eBoundLoc0, eFreePos1) = cosPhi;
    // jacToCurv(eBoundLoc0, eFreePos2) = 0;
    jacToCurv(eBoundLoc1, eFreePos0) = -cosPhi * cosTheta;
    jacToCurv(eBoundLoc1, eFreePos1) = -sinPhi * cosTheta;
    jacToCurv(eBoundLoc1, eFreePos2) = sinTheta;
  } else {
    // Under grazing incidence to z, the above coordinate system definition
    // becomes numerically unstable, and we need to switch to another one
    const double c = std::hypot(y, z);
    const double invC = 1. / c;
    // jacToCurv(eBoundLoc0, eFreePos0) = 0;
    jacToCurv(eBoundLoc0, eFreePos1) = -z * invC;
    jacToCurv(eBoundLoc0, eFreePos2) = y * invC;
    jacToCurv(eBoundLoc1, eFreePos0) = c;
    jacToCurv(eBoundLoc1, eFreePos1) = -x * y * invC;
    jacToCurv(eBoundLoc1, eFreePos2) = -x * z * invC;
  }
  // Time parameter
  jacToCurv(eBoundTime, eFreeTime) = 1.;
  // Directional and momentum parameters for curvilinear
  jacToCurv.block<2, 3>(eBoundPhi, eFreeDir0) =
      freeToSphericalDirectionJacobian(direction);
  jacToCurv(eBoundQOverP, eFreeQOverP) = 1.;

  return jacToCurv;
}

/// @brief This function calculates the full jacobian from local parameters at
///        the start surface to bound parameters at the final surface
///
/// @note Modifications of the jacobian related to the
/// projection onto a surface is considered. Since a variation of the start
/// parameters within a given uncertainty would lead to a variation of the end
/// parameters, these need to be propagated onto the target surface. This an
/// approximated approach to treat the (assumed) small change.
///
/// @param [in] geoContext The geometry Context
/// @param [in] freeParameters Free, nominal parametrisation
/// @param [in] boundToFreeJacobian The projection jacobian from start local
///        to start free parameters
/// @param [in] freeTransportJacobian The transport jacobian from start free to
///        final free parameters
/// @param [in] freeToPathDerivatives Path length derivatives of the final free
///        parameters
/// @param [out] fullTransportJacobian The full jacobian from start local to
///        bound parameters at the final surface
/// @param [in] surface The final surface onto which the projection should be
///        performed
void boundToBoundJacobian(const GeometryContext& geoContext,
                          const FreeVector& freeParameters,
                          const BoundToFreeMatrix& boundToFreeJacobian,
                          const FreeMatrix& freeTransportJacobian,
                          FreeToBoundMatrix& freeToBoundJacobian,
                          const FreeVector& freeToPathDerivatives,
                          BoundMatrix& fullTransportJacobian,
                          const Surface& surface) {
  // Calculate the derivative of path length at the final surface or the
  // point-of-closest approach w.r.t. free parameters
  const FreeToPathMatrix freeToPath =
      surface.freeToPathDerivative(geoContext, freeParameters);
  // Calculate the jacobian from free to bound at the final surface
  freeToBoundJacobian = surface.freeToBoundJacobian(geoContext, freeParameters);
  // Calculate the full jacobian from the local/bound parameters at the start
  // surface to local/bound parameters at the final surface
  // @note jac(locA->locB) = jac(gloB->locB)*(1+
  // pathCorrectionFactor(gloB))*jacTransport(gloA->gloB) *jac(locA->gloA)
  fullTransportJacobian =
      freeToBoundJacobian *
      (FreeMatrix::Identity() + freeToPathDerivatives * freeToPath) *
      freeTransportJacobian * boundToFreeJacobian;
}

/// @brief This function calculates the full jacobian from local parameters at
/// the start surface to final curvilinear parameters
///
/// @note Modifications of the jacobian related to the
/// projection onto a curvilinear surface is considered. Since a variation of
/// the start parameters within a given uncertainty would lead to a variation of
/// the end parameters, these need to be propagated onto the target surface.
/// This is an approximated approach to treat the (assumed) small change.
///
/// @param [in] direction Normalised direction vector
/// @param [in] boundToFreeJacobian The projection jacobian from local start
///        to global final parameters
/// @param [in] freeTransportJacobian The transport jacobian from start free to
///        final free parameters
/// @param [in] freeToPathDerivatives Path length derivatives of the final free
///        parameters
/// @param [in, out] jacFull The full jacobian from start local to curvilinear
///        parameters
///
/// @note The parameter @p surface is only required if projected to bound
/// parameters. In the case of curvilinear parameters the geometry and the
/// position is known and the calculation can be simplified
void boundToCurvilinearJacobian(const Vector3& direction,
                                const BoundToFreeMatrix& boundToFreeJacobian,
                                const FreeMatrix& freeTransportJacobian,
                                FreeToBoundMatrix& freeToBoundJacobian,
                                const FreeVector& freeToPathDerivatives,
                                BoundMatrix& fullTransportJacobian) {
  // Calculate the jacobian from global to local at the curvilinear surface
  freeToBoundJacobian = freeToCurvilinearJacobian(direction);

  // Update the jacobian to include the derivative of the path length at the
  // curvilinear surface w.r.t. the free parameters
  freeToBoundJacobian.topLeftCorner<6, 3>() +=
      (freeToBoundJacobian * freeToPathDerivatives) *
      (-1.0 * direction).transpose();

  // Calculate the full jocobian from the local parameters at the start surface
  // to curvilinear parameters
  // @note jac(locA->locB) = jac(gloB->locB)*(1+
  // pathCorrectionFactor(gloB))*jacTransport(gloA->gloB) *jac(locA->gloA)
  fullTransportJacobian =
      blockedMult(freeToBoundJacobian,
                  blockedMult(freeTransportJacobian, boundToFreeJacobian));
}

/// @brief This function reinitialises the state members required for the
/// covariance transport
///
/// @param [in] geoContext The geometry context
/// @param [in, out] freeTransportJacobian The transport jacobian from start
///        free to final free parameters
/// @param [in, out] freeToPathDerivatives Path length derivatives of the free,
///        nominal parameters
/// @param [in, out] boundToFreeJacobian Projection jacobian of the last bound
///        parametrisation to free parameters
/// @param [in] freeParameters Free, nominal parametrisation
/// @param [in] surface The reference surface of the local parametrisation
Result<void> reinitializeJacobians(const GeometryContext& geoContext,
                                   FreeMatrix& freeTransportJacobian,
                                   FreeVector& freeToPathDerivatives,
                                   BoundToFreeMatrix& boundToFreeJacobian,
                                   const FreeVector& freeParameters,
                                   const Surface& surface) {
  using VectorHelpers::phi;
  using VectorHelpers::theta;

  // Reset the jacobians
  freeTransportJacobian = FreeMatrix::Identity();
  freeToPathDerivatives = FreeVector::Zero();

  // Get the local position
  const Vector3 position = freeParameters.segment<3>(eFreePos0);
  const Vector3 direction = freeParameters.segment<3>(eFreeDir0);
  auto lpResult = surface.globalToLocal(geoContext, position, direction);
  if (!lpResult.ok()) {
    return lpResult.error();
  }
  // Transform from free to bound parameters
  Result<BoundVector> boundParameters = detail::transformFreeToBoundParameters(
      freeParameters, surface, geoContext);
  if (!boundParameters.ok()) {
    return boundParameters.error();
  }
  // Reset the jacobian from local to global
  boundToFreeJacobian =
      surface.boundToFreeJacobian(geoContext, *boundParameters);
  return Result<void>::success();
}

/// @brief This function reinitialises the state members required for the
/// covariance transport
///
/// @param [in, out] freeTransportJacobian The transport jacobian from start
///        free to final free parameters
/// @param [in, out] freeToPathDerivatives Path length derivatives of the free,
///        nominal parameters
/// @param [in, out] boundToFreeJacobian Projection jacobian of the last bound
///        parametrisation to free parameters
/// @param [in] direction Normalised direction vector
void reinitializeJacobians(FreeMatrix& freeTransportJacobian,
                           FreeVector& freeToPathDerivatives,
                           BoundToFreeMatrix& boundToFreeJacobian,
                           const Vector3& direction) {
  // Reset the jacobians
  freeTransportJacobian = FreeMatrix::Identity();
  freeToPathDerivatives = FreeVector::Zero();
  boundToFreeJacobian = BoundToFreeMatrix::Zero();

  auto [cosPhi, sinPhi, cosTheta, sinTheta] =
      VectorHelpers::evaluateTrigonomics(direction);

  boundToFreeJacobian(eFreePos0, eBoundLoc0) = -sinPhi;
  boundToFreeJacobian(eFreePos0, eBoundLoc1) = -cosPhi * cosTheta;
  boundToFreeJacobian(eFreePos1, eBoundLoc0) = cosPhi;
  boundToFreeJacobian(eFreePos1, eBoundLoc1) = -sinPhi * cosTheta;
  boundToFreeJacobian(eFreePos2, eBoundLoc1) = sinTheta;
  boundToFreeJacobian(eFreeTime, eBoundTime) = 1;
  boundToFreeJacobian.block<3, 2>(eFreeDir0, eBoundPhi) =
      sphericalToFreeDirectionJacobian(direction);
  boundToFreeJacobian(eFreeQOverP, eBoundQOverP) = 1;
}
}  // namespace

namespace detail {

Result<BoundState> boundState(
    const GeometryContext& geoContext, BoundSquareMatrix& covarianceMatrix,
    BoundMatrix& fullTransportJacobian, FreeMatrix& freeTransportJacobian,
    FreeVector& freeToPathDerivatives, BoundToFreeMatrix& boundToFreeJacobian,
    std::optional<FreeMatrix>& additionalFreeCovariance, FreeVector& parameters,
    const ParticleHypothesis& particleHypothesis, bool covTransport,
    double accumulatedPath, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) {
  // Covariance transport
  std::optional<BoundSquareMatrix> cov = std::nullopt;
  if (covTransport) {
    // Initialize the jacobian from start local to final local
    fullTransportJacobian = BoundMatrix::Identity();
    // Calculate the jacobian and transport the covarianceMatrix to final local.
    // Then reinitialize the freeTransportJacobian, freeToPathDerivatives and
    // the boundToFreeJacobian
    transportCovarianceToBound(
        geoContext, covarianceMatrix, fullTransportJacobian,
        freeTransportJacobian, freeToPathDerivatives, boundToFreeJacobian,
        additionalFreeCovariance, parameters, surface, freeToBoundCorrection);
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
      fullTransportJacobian, accumulatedPath);
}

CurvilinearState curvilinearState(
    BoundSquareMatrix& covarianceMatrix, BoundMatrix& fullTransportJacobian,
    FreeMatrix& freeTransportJacobian, FreeVector& freeToPathDerivatives,
    BoundToFreeMatrix& boundToFreeJacobian,
    std::optional<FreeMatrix>& additionalFreeCovariance,
    const FreeVector& parameters, const ParticleHypothesis& particleHypothesis,
    bool covTransport, double accumulatedPath) {
  const Vector3& direction = parameters.segment<3>(eFreeDir0);

  // Covariance transport
  std::optional<BoundSquareMatrix> cov = std::nullopt;
  if (covTransport) {
    // Initialize the jacobian from start local to final local
    fullTransportJacobian = BoundMatrix::Identity();
    // Calculate the jacobian and transport the covarianceMatrix to final local.
    // Then reinitialize the freeTransportJacobian, freeToPathDerivatives and
    // the boundToFreeJacobian
    transportCovarianceToCurvilinear(covarianceMatrix, fullTransportJacobian,
                                     freeTransportJacobian,
                                     freeToPathDerivatives, boundToFreeJacobian,
                                     additionalFreeCovariance, direction);
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
  return std::make_tuple(std::move(curvilinearParams), fullTransportJacobian,
                         accumulatedPath);
}

void transportCovarianceToBound(
    const GeometryContext& geoContext, BoundSquareMatrix& boundCovariance,
    BoundMatrix& fullTransportJacobian, FreeMatrix& freeTransportJacobian,
    FreeVector& freeToPathDerivatives, BoundToFreeMatrix& boundToFreeJacobian,
    std::optional<FreeMatrix>& additionalFreeCovariance,
    FreeVector& freeParameters, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) {
  FreeToBoundMatrix freeToBoundJacobian;

  // Calculate the full jacobian from local parameters at the start surface to
  // current bound parameters
  boundToBoundJacobian(geoContext, freeParameters, boundToFreeJacobian,
                       freeTransportJacobian, freeToBoundJacobian,
                       freeToPathDerivatives, fullTransportJacobian, surface);

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

  if (additionalFreeCovariance) {
    boundCovariance += freeToBoundJacobian * (*additionalFreeCovariance) *
                       freeToBoundJacobian.transpose();
  }

  // Reinitialize jacobian components:
  // ->The transportJacobian is reinitialized to Identity
  // ->The derivatives is reinitialized to Zero
  // ->The boundToFreeJacobian is initialized to that at the current surface
  reinitializeJacobians(geoContext, freeTransportJacobian,
                        freeToPathDerivatives, boundToFreeJacobian,
                        freeParameters, surface);
}

void transportCovarianceToCurvilinear(
    BoundSquareMatrix& boundCovariance, BoundMatrix& fullTransportJacobian,
    FreeMatrix& freeTransportJacobian, FreeVector& freeToPathDerivatives,
    BoundToFreeMatrix& boundToFreeJacobian,
    std::optional<FreeMatrix>& additionalFreeCovariance,
    const Vector3& direction) {
  FreeToBoundMatrix freeToBoundJacobian;

  // Calculate the full jacobian from local parameters at the start surface to
  // current curvilinear parameters
  boundToCurvilinearJacobian(direction, boundToFreeJacobian,
                             freeTransportJacobian, freeToBoundJacobian,
                             freeToPathDerivatives, fullTransportJacobian);

  // Apply the actual covariance transport to get covariance of the current
  // curvilinear parameters
  boundCovariance = fullTransportJacobian * boundCovariance *
                    fullTransportJacobian.transpose();

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

Acts::Result<Acts::BoundTrackParameters> boundToBoundConversion(
    const GeometryContext& gctx, const BoundTrackParameters& boundParameters,
    const Surface& targetSurface, const Vector3& bField) {
  const auto& sourceSurface = boundParameters.referenceSurface();

  Acts::FreeVector freePars = Acts::detail::transformBoundToFreeParameters(
      sourceSurface, gctx, boundParameters.parameters());

  auto res = Acts::detail::transformFreeToBoundParameters(freePars,
                                                          targetSurface, gctx);

  if (!res.ok()) {
    return res.error();
  }
  Acts::BoundVector parOut = *res;

  std::optional<Acts::BoundMatrix> covOut = std::nullopt;

  if (boundParameters.covariance().has_value()) {
    Acts::BoundToFreeMatrix boundToFreeJacobian =
        sourceSurface.boundToFreeJacobian(gctx, boundParameters.parameters());

    Acts::FreeMatrix freeTransportJacobian = FreeMatrix::Identity();

    FreeVector freeToPathDerivatives = FreeVector::Zero();
    freeToPathDerivatives.head<3>() = freePars.segment<3>(eFreeDir0);

    freeToPathDerivatives.segment<3>(eFreeDir0) =
        bField.cross(freePars.segment<3>(eFreeDir0));

    BoundMatrix boundToBoundJac = detail::boundToBoundTransportJacobian(
        gctx, freePars, boundToFreeJacobian, freeTransportJacobian,
        freeToPathDerivatives, targetSurface);

    covOut = boundToBoundJac * (*boundParameters.covariance()) *
             boundToBoundJac.transpose();
  }

  return Acts::BoundTrackParameters{targetSurface.getSharedPtr(), parOut,
                                    covOut,
                                    boundParameters.particleHypothesis()};
}

}  // namespace detail
}  // namespace Acts
