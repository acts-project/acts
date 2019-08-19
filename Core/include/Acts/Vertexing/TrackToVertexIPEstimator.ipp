// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PerigeeSurface.hpp"

template <typename input_track_t, typename propagator_t, typename action_list_t,
          typename aborter_list_t>
Acts::Result<std::unique_ptr<Acts::ImpactParametersAndSigma>>
Acts::TrackToVertexIPEstimator<
    input_track_t, propagator_t, action_list_t,
    aborter_list_t>::estimate(const BoundParameters& track,
                              const Vertex<input_track_t>& vtx) const {
  // estimating the d0 and its significance by propagating the trajectory state
  // towards
  // the vertex position. By this time the vertex should NOT contain this
  // trajectory anymore

  const Vector3D& lp = vtx.position();

  const std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(lp);

  // Set the propagation direction to be backward as needed below
  auto options = m_cfg.pOptions;
  options.direction = backward;

  // Do the propagation to linPoint
  auto result = m_cfg.propagator->propagate(track, *perigeeSurface, options);

  if (!result.ok()) {
    return result.error();
  }

  const auto& propRes = *result;
  const auto& params = propRes.endParameters->parameters();
  const double d0 = params[ParID_t::eLOC_D0];
  const double z0 = params[ParID_t::eLOC_Z0];
  const double phi = params[ParID_t::ePHI];
  const double theta = params[ParID_t::eTHETA];

  ActsSymMatrixD<2> vrtXYCov = vtx.covariance().template block<2, 2>(0, 0);

  // Covariance of perigee parameters after propagation to perigee surface
  const auto& perigeeCov = *(propRes.endParameters->covariance());

  ActsVectorD<2> d0JacXY(-std::sin(phi), std::cos(phi));

  std::unique_ptr<ImpactParametersAndSigma> newIPandSigma =
      std::make_unique<ImpactParametersAndSigma>();
  newIPandSigma->IPd0 = d0;
  double d0_PVcontrib = d0JacXY.transpose() * (vrtXYCov * d0JacXY);
  if (d0_PVcontrib >= 0) {
    newIPandSigma->sigmad0 = std::sqrt(
        d0_PVcontrib + perigeeCov(ParID_t::eLOC_D0, ParID_t::eLOC_D0));
    newIPandSigma->PVsigmad0 = std::sqrt(d0_PVcontrib);
  } else {
    newIPandSigma->sigmad0 =
        std::sqrt(perigeeCov(ParID_t::eLOC_D0, ParID_t::eLOC_D0));
    newIPandSigma->PVsigmad0 = 0;
  }

  ActsSymMatrixD<2> covPerigeeZ0Theta;
  covPerigeeZ0Theta(0, 0) = perigeeCov(ParID_t::eLOC_Z0, ParID_t::eLOC_Z0);
  covPerigeeZ0Theta(0, 1) = perigeeCov(ParID_t::eLOC_Z0, ParID_t::eTHETA);
  covPerigeeZ0Theta(1, 0) = perigeeCov(ParID_t::eTHETA, ParID_t::eLOC_Z0);
  covPerigeeZ0Theta(1, 1) = perigeeCov(ParID_t::eTHETA, ParID_t::eTHETA);

  double vtxZZCov = vtx.covariance()(eZ, eZ);

  ActsVectorD<2> z0JacZ0Theta(std::sin(theta), z0 * std::cos(theta));

  if (vtxZZCov >= 0) {
    newIPandSigma->IPz0SinTheta = z0 * std::sin(theta);
    newIPandSigma->sigmaz0SinTheta = std::sqrt(
        z0JacZ0Theta.transpose() * (covPerigeeZ0Theta * z0JacZ0Theta) +
        std::sin(theta) * vtxZZCov * std::sin(theta));

    newIPandSigma->PVsigmaz0SinTheta =
        std::sqrt(std::sin(theta) * vtxZZCov * std::sin(theta));
    newIPandSigma->IPz0 = z0;
    newIPandSigma->sigmaz0 = std::sqrt(vtxZZCov + perigeeCov(eZ, eZ));
    newIPandSigma->PVsigmaz0 = std::sqrt(vtxZZCov);
  } else {
    ACTS_WARNING(
        "Contribution to z0_err from PV is negative: Error in PV "
        "error matrix! Removing contribution from PV");
    newIPandSigma->IPz0SinTheta = z0 * std::sin(theta);
    double sigma2z0sinTheta =
        (z0JacZ0Theta.transpose() * (covPerigeeZ0Theta * z0JacZ0Theta));
    newIPandSigma->sigmaz0SinTheta = std::sqrt(sigma2z0sinTheta);
    newIPandSigma->PVsigmaz0SinTheta = 0;

    newIPandSigma->IPz0 = z0;
    newIPandSigma->sigmaz0 = std::sqrt(perigeeCov(eZ, eZ));
    newIPandSigma->PVsigmaz0 = 0;
  }

  return std::move(newIPandSigma);
}
