// This file is part of the Acts project.
//
// Copyright (C) 2021-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"

#include "Acts/Surfaces/Surface.hpp"

/// Get the non-linearity corrected bound parameters and its covariance
std::optional<std::tuple<Acts::BoundVector, Acts::BoundSymMatrix>>
Acts::detail::CorrectedFreeToBoundTransformer::operator()(
    const Acts::FreeVector& freeParams,
    const Acts::FreeSymMatrix& freeCovariance, const Acts::Surface& surface,
    const Acts::GeometryContext& geoContext, NavigationDirection navDir) {
  // Get the incidence angle
  Vector3 dir = freeParams.segment<3>(eFreeDir0);
  Vector3 normal = surface.normal(geoContext);
  ActsScalar absCosIncidenceAng = std::abs(dir.dot(normal));
  // No correction if the incidentAngle is small enough (not necessary ) or too
  // large (correction could be invalid).
  if (absCosIncidenceAng < cosIncidentAngleMinCutoff or
      absCosIncidenceAng > cosIncidentAngleMaxCutoff) {
    return std::nullopt;
  }

  size_t sampleSize = 2 * eFreeSize + 1;
  std::vector<std::tuple<FreeVector, ActsScalar, ActsScalar>> sampledFreeParams;
  sampledFreeParams.reserve(sampleSize);

  // Initialize the covariance sqrt root matrix
  FreeSymMatrix covSqrt = FreeSymMatrix::Zero();
  // Get the covariance sqrt root matrix
  Eigen::JacobiSVD<FreeSymMatrix> svd(
      freeCovariance, Eigen::ComputeThinU | Eigen::ComputeThinV);
  // U*S*V^-1
  auto S = svd.singularValues();
  auto U = svd.matrixU();
  FreeMatrix D = FreeMatrix::Zero();
  for (unsigned i = 0; i < eFreeSize; ++i) {
    if (S(i) == 0) {
      continue;
    }
    D(i, i) = std::sqrt(S(i));
  }

  FreeMatrix UP = FreeMatrix::Zero();
  for (unsigned i = 0; i < eFreeSize; ++i) {
    for (unsigned j = 0; j < eFreeSize; ++j) {
      UP(i, j) = U(i, j);
    }
  }
  covSqrt = UP * D;

  double kappa = alpha * alpha * eFreeSize;
  double gamma = std::sqrt(kappa);
  double lambda = kappa - eFreeSize;

  // Sample the free parameters
  // 1. the baseline parameter
  sampledFreeParams.push_back({freeParams, lambda / kappa,
                               lambda / kappa + (1.0 - alpha * alpha + beta)});
  // 2. the shifted parameters
  for (unsigned i = 0; i < eFreeSize; ++i) {
    sampledFreeParams.push_back(
        {freeParams + covSqrt.col(i) * gamma, 0.5 / kappa, 0.5 / kappa});
    sampledFreeParams.push_back(
        {freeParams - covSqrt.col(i) * gamma, 0.5 / kappa, 0.5 / kappa});
  }

  // Initialize the mean of the bound parameters
  BoundVector bpMean = BoundVector::Zero();
  // Initialize the bound covariance
  BoundSymMatrix bv = BoundSymMatrix::Zero();

  // The transformed bound parameters
  std::vector<std::pair<BoundVector, ActsScalar>> transformedBoundParams;

  // 1. The nominal one
  const auto& [params_, mweight_, cweight_] = sampledFreeParams[0];
  // Transform the free to bound
  auto nominalRes =
      detail::transformFreeToBoundParameters(params_, surface, geoContext);
  if (not nominalRes.ok()) {
    return std::nullopt;
  }
  auto nominalBound = nominalRes.value();
  transformedBoundParams.push_back({nominalBound, cweight_});
  bpMean = bpMean + mweight_ * nominalBound;

  // 2. Loop over the rest sample points of the free parameters to get the
  // weighted bound parameters
  for (unsigned i = 1; i < sampledFreeParams.size(); ++i) {
    const auto& [params, mweight, cweight] = sampledFreeParams[i];
    FreeVector correctedFreeParams = params;

    // Reintersect to get the corrected free params without boundary check
    SurfaceIntersection intersection =
        surface.intersect(geoContext, params.segment<3>(eFreePos0),
                          navDir * params.segment<3>(eFreeDir0), false);
    correctedFreeParams.segment<3>(eFreePos0) =
        intersection.intersection.position;

    // Transform the free to bound
    auto result = detail::transformFreeToBoundParameters(correctedFreeParams,
                                                         surface, geoContext);
    if (not result.ok()) {
      return std::nullopt;
    }

    auto bp = result.value();
    transformedBoundParams.push_back({bp, cweight});
    bpMean = bpMean + mweight * bp;
  }

  if (transformedBoundParams.empty()) {
    return std::nullopt;
  }

  // Get the weighted bound covariance
  for (unsigned isample = 0; isample < sampleSize; ++isample) {
    BoundVector bSigma = transformedBoundParams[isample].first - bpMean;

    bv = bv +
         transformedBoundParams[isample].second * bSigma * bSigma.transpose();
  }

  return std::make_tuple(bpMean, bv);
}
