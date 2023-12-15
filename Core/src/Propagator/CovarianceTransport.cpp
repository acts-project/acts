// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/CovarianceTransport.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Propagator/detail/JacobianEngine.hpp"

Acts::CovarianceCache::CovarianceCache(const GeometryContext& gctx,
                                       const Surface& surface,
                                       const Vector3& position,
                                       const BoundVector& boundParameters,
                                       const BoundSquareMatrix& boundCovariance)
    : applyTransport(true), atSurface(&surface), atPosition(position) {
  covariance.template emplace<BoundSquareMatrix>(boundCovariance);
  boundToFreeJacobian = surface.boundToFreeJacobian(gctx, boundParameters);
}

Acts::CovarianceCache::CovarianceCache(const Vector3& position,
                                       const Vector3& direction,
                                       const BoundSquareMatrix& boundCovariance)
    : applyTransport(true), atPosition(position) {
  covariance.emplace<BoundSquareMatrix>(boundCovariance);
  boundToFreeJacobian = detail::curvilinearToFreeJacobian(direction);
}

Acts::CovarianceCache::CovarianceCache(const FreeVector& freeParameters,
                                       const FreeSquareMatrix& freeCovariance)
    : applyTransport(true), atPosition(freeParameters.segment<3>(eFreePos0)) {
  covariance.emplace<FreeSquareMatrix>(freeCovariance);
  anglesToDirectionJacobian =
      detail::anglesToDirectionJacobian(freeParameters.segment<3>(eFreeDir0));
  directionToAnglesJacobian =
      detail::directionToAnglesJacobian(freeParameters.segment<3>(eFreeDir0));
}

std::tuple<Acts::VariantCovariance, Acts::VariantTransportJacobian>
Acts::transportCovarianceToBound(const GeometryContext& gctx,
                                 const Surface& surface,
                                 const FreeVector& freeParameters,
                                 CovarianceCache& cCache) {
  if (cCache.boundToFreeJacobian.has_value()) {
    // Create the full transport jacobian: bound/curvilinear to bound
    const auto ftJacobian = detail::boundToBoundTransportJacobian(
        gctx, freeParameters, cCache.boundToFreeJacobian.value(),
        cCache.freeTransportJacobian, cCache.freeToPathDerivatives, surface);
    // Perform the transport
    const auto& covariance = std::get<BoundSquareMatrix>(cCache.covariance);
    BoundSquareMatrix newCovariance =
        ftJacobian * covariance * ftJacobian.transpose();
    return {newCovariance, ftJacobian};
  } else {
    // Create the full transport jacobian: free to bound
    const auto ftJacobian = detail::freeToBoundTransportJacobian(
        gctx, freeParameters, cCache.directionToAnglesJacobian.value(),
        cCache.anglesToDirectionJacobian.value(), cCache.freeTransportJacobian,
        cCache.freeToPathDerivatives, surface);
    // Perform the transport
    const auto& covariance = std::get<FreeSquareMatrix>(cCache.covariance);
    BoundSquareMatrix newCovariance =
        ftJacobian * covariance * ftJacobian.transpose();
    return {newCovariance, ftJacobian};
  }
}

std::tuple<Acts::VariantCovariance, Acts::VariantTransportJacobian>
Acts::transportCovarianceToCurvilinear(const Vector3& direction,
                                       CovarianceCache& cCache) {
  if (cCache.boundToFreeJacobian.has_value()) {
    // Create the full transport jacobian: bound/curvilinear to
    // curvilinear
    auto ftJacobian = detail::boundToCurvilinearTransportJacobian(
        direction, cCache.boundToFreeJacobian.value(),
        cCache.freeTransportJacobian, cCache.freeToPathDerivatives);
    // Perform the transport
    const auto& covariance = std::get<BoundSquareMatrix>(cCache.covariance);
    BoundSquareMatrix newCovariance =
        ftJacobian * covariance * ftJacobian.transpose();
    return {newCovariance, ftJacobian};
  } else {
    // Create the full transport jacobian: free to curvilinear
    const auto ftJacobian = detail::freeToCurvilinearTransportJacobian(
        direction, cCache.directionToAnglesJacobian.value(),
        cCache.anglesToDirectionJacobian.value(), cCache.freeTransportJacobian,
        cCache.freeToPathDerivatives);
    // Perform the transport
    const auto& covariance = std::get<FreeSquareMatrix>(cCache.covariance);
    BoundSquareMatrix newCovariance =
        ftJacobian * covariance * ftJacobian.transpose();
    return {newCovariance, ftJacobian};
  }
}

std::tuple<Acts::VariantCovariance, Acts::VariantTransportJacobian>
Acts::transportCovarianceToFree(CovarianceCache& cCache) {
  if (cCache.boundToFreeJacobian.has_value()) {
    // Create the full transport jacobian: bound/curvilinear to free
    const auto ftJacobian = detail::boundToFreeTransportJacobian(
        cCache.boundToFreeJacobian.value(), cCache.freeTransportJacobian);
    // Perform the transport
    const auto& covariance = std::get<BoundSquareMatrix>(cCache.covariance);
    FreeSquareMatrix newCovariance =
        ftJacobian * covariance * ftJacobian.transpose();
    return {newCovariance, ftJacobian};
  } else {
    // Create the full transport jacobian: free to free
    const auto ftJacobian = detail::freeToFreeTransportJacobian(
        cCache.directionToAnglesJacobian.value(),
        cCache.anglesToDirectionJacobian.value(), cCache.freeTransportJacobian);
    // Perform the transport
    const auto& covariance = std::get<FreeSquareMatrix>(cCache.covariance);
    FreeSquareMatrix newCovariance =
        ftJacobian * covariance * ftJacobian.transpose();
    return {newCovariance, ftJacobian};
  }
}

Acts::Result<Acts::BoundTrackParameters> Acts::boundToBoundConversion(
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
