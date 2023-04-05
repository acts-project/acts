// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Charge.hpp"
#include "Acts/EventData/SingleBoundTrackParameters.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/CovarianceTransport.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <Eigen/src/Core/util/Memory.h>
#include <boost/graph/graph_traits.hpp>
#include <edm4hep/Track.h>
#include <edm4hep/TrackState.h>

#include "edm4hep/MutableTrack.h"

namespace Acts {
namespace EDM4hepUtil {

static constexpr std::int32_t EDM4HEP_ACTS_POSITION_TYPE = 42;

namespace detail {
struct Parameters {
  Acts::ActsVector<5> values;
  double time;
  std::optional<Acts::ActsSymMatrix<5>> covariance;
  std::shared_ptr<const Acts::Surface> surface;
};

ActsSymMatrix<5> jacobianToEdm4hep(double theta, double qOverP, double Bz) {
  // Calculate jacobian from our internal parametrization (d0, z0, phi, theta,
  // q/p) to the LCIO / edm4hep (see:
  // https://bib-pubdb1.desy.de/record/81214/files/LC-DET-2006-004%5B1%5D.pdf)
  // one (d0, z0, phi, tan(lambda), omega). Top left 3x3 matrix in the
  // jacobian is 1. Bottom right 2x2 matrix is:
  //
  // [  d                                 ]
  // [------(cot(theta))         0        ]
  // [dtheta                              ]
  // [                                    ]
  // [  d   /  B*q/p   \   d  /  B*q/p   \]
  // [------|----------|  ----|----------|]
  // [dtheta\sin(theta)/  dq/p\sin(theta)/]
  //
  // =
  //
  // [     2                        ]
  // [- cot (theta) - 1       0     ]
  // [                              ]
  // [-B*q/p*cos(theta)       B     ]
  // [------------------  ----------]
  // [      2             sin(theta)]
  // [   sin (theta)                ]

  ActsSymMatrix<5> J;
  J.setZero();
  J(0, 0) = 1;
  J(1, 1) = 1;
  J(2, 2) = 1;
  double cotTheta = std::tan(M_PI_2 + theta);
  J(3, 3) = -cotTheta * cotTheta - 1;  // d(tanLambda) / dTheta
  J(4, 4) = Bz / std::sin(theta);      // dOmega / d(qop)
  double sinTheta = std::sin(theta);
  J(4, 3) = -Bz * qOverP *
            std::cos(theta / (sinTheta * sinTheta));  // dOmega / dTheta
  return J;
}

template <typename charge_t>
Parameters convertTrackParametersToEdm4hep(
    const Acts::GeometryContext& gctx, double Bz,
    const SingleBoundTrackParameters<charge_t>& params) {
  Acts::Vector3 global = params.referenceSurface().localToGlobal(
      gctx, params.parameters().template head<2>(), params.momentum());

  std::shared_ptr<const Acts::Surface> refSurface =
      params.referenceSurface().getSharedPtr();

  Acts::BoundVector targetPars = params.parameters();
  std::optional<Acts::FreeVector> freePars;

  auto makeFreePars = [&]() {
    return Acts::detail::transformBoundToFreeParameters(
        params.referenceSurface(), gctx, params.parameters());
  };

  // If the reference surface is a perigee surface, we use that. Otherwise
  // we create a new perigee surface at the global positon of the track
  // parameters.
  if (dynamic_cast<const Acts::PerigeeSurface*>(refSurface.get()) == nullptr) {
    refSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(global);

    // We need to convert to the target parameters
    // Keep the free parameters around we might need them for the covariance
    // conversion
    freePars = makeFreePars();
    targetPars = Acts::detail::transformFreeToBoundParameters(freePars.value(),
                                                              *refSurface, gctx)
                     .value();
  }

  Parameters result;
  result.surface = refSurface;

  // Only run covariance conversion if we have a covariance input
  if (params.covariance()) {
    auto boundToFree =
        refSurface->boundToFreeJacobian(gctx, params.parameters());
    Acts::FreeMatrix freeCov =
        boundToFree * params.covariance().value() * boundToFree.transpose();

    // ensure we have free pars
    if (!freePars.has_value()) {
      freePars = makeFreePars();
    }

    Acts::CovarianceCache covCache{freePars.value(), freeCov};
    auto [varNewCov, varNewJac] = Acts::transportCovarianceToBound(
        gctx, *refSurface, freePars.value(), covCache);
    auto targetCov = std::get<Acts::BoundSymMatrix>(varNewCov);

    Acts::ActsSymMatrix<5> J = jacobianToEdm4hep(targetPars[eBoundTheta],
                                                 targetPars[eBoundQOverP], Bz);
    Acts::ActsSymMatrix<5> cIn = targetCov.template topLeftCorner<5, 5>();
    result.covariance = J * cIn * J.transpose();
  }

  result.values[0] = targetPars[Acts::eBoundLoc0];
  result.values[1] = targetPars[Acts::eBoundLoc1];
  result.values[2] = targetPars[Acts::eBoundPhi];
  result.values[3] = std::tan(M_PI_2 - targetPars[Acts::eBoundTheta]);
  result.time = targetPars[Acts::eBoundTime];

  result.values[4] = params.charge() * Bz / params.transverseMomentum();

  return result;
}

template <typename charge_t>
SingleBoundTrackParameters<charge_t> convertTrackParametersFromEdm4hep(
    double Bz, const Parameters& params) {
  // [     d      /                     pi\                                  ]
  // [------------|-atan(\tan\lambda) + --|                 0                ]
  // [d\tan\lambda\                     2 /                                  ]
  // [                                                                       ]
  // [     d      /         \Omega        \     d   /         \Omega        \]
  // [------------|-----------------------|  -------|-----------------------|]
  // [d\tan\lambda|     __________________|  d\Omega|     __________________|]
  // [            |    /            2     |         |    /            2     |]
  // [            \B*\/  \tan\lambda  + 1 /         \B*\/  \tan\lambda  + 1 /]
  //
  // =
  //
  // [         -1                                     ]
  // [   ----------------                 0           ]
  // [              2                                 ]
  // [   \tan\lambda  + 1                             ]
  // [                                                ]
  // [  -\Omega*\tan\lambda               1           ]
  // [-----------------------  -----------------------]
  // [                    3/2       __________________]
  // [  /           2    \         /            2     ]
  // [B*\\tan\lambda  + 1/     B*\/  \tan\lambda  + 1 ]

  ActsSymMatrix<5> J;
  J.setIdentity();
  J(3, 3) = -1 / (params.values[3] * params.values[3] + 1);
  J(4, 3) = -1 * params.values[3] * params.values[4] /
            (Bz * std::pow(params.values[3] * params.values[3] + 1, 3. / 2.));
  J(4, 4) = 1 / (Bz * std::sqrt(params.values[3] * params.values[3] + 1));

  BoundVector targetPars;
  // @TODO: Double check this
  BoundMatrix cov;
  cov.setIdentity();
  cov.template topLeftCorner<5, 5>() =
      J * params.covariance.value() * J.transpose();

  targetPars[eBoundLoc0] = params.values[0];
  targetPars[eBoundLoc1] = params.values[1];
  targetPars[eBoundPhi] = params.values[2];
  targetPars[eBoundTheta] = M_PI_2 - std::atan(params.values[3]);
  targetPars[eBoundQOverP] =
      params.values[4] * std::sin(targetPars[eBoundTheta]) / Bz;
  targetPars[eBoundTime] = params.time;

  return {params.surface, targetPars, cov};
}
}  // namespace detail

template <typename track_container_t, typename track_state_container_t,
          template <typename> class holder_t>
void writeTrack(
    const Acts::GeometryContext& gctx,
    Acts::TrackProxy<track_container_t, track_state_container_t, holder_t, true>
        track,
    edm4hep::MutableTrack to, double Bz) {
  to.setChi2(track.chi2());
  to.setNdf(track.nDoF());

  std::vector<edm4hep::TrackState> outTrackStates;
  outTrackStates.reserve(track.nTrackStates());

  auto setParameters = [Bz](edm4hep::TrackState& trackState,
                            const detail::Parameters& params) {
    trackState.D0 = params.values[0];
    trackState.Z0 = params.values[1];
    trackState.phi = params.values[2];
    trackState.tanLambda = params.values[3];
    trackState.omega = params.values[4];
    trackState.time = params.time;

    if (params.covariance) {
      const auto& c = params.covariance.value();

      trackState.covMatrix = {
          static_cast<float>(c(0, 0)), static_cast<float>(c(1, 0)),
          static_cast<float>(c(1, 1)), static_cast<float>(c(2, 0)),
          static_cast<float>(c(2, 1)), static_cast<float>(c(2, 2)),
          static_cast<float>(c(3, 0)), static_cast<float>(c(3, 1)),
          static_cast<float>(c(3, 2)), static_cast<float>(c(3, 3)),
          static_cast<float>(c(4, 0)), static_cast<float>(c(4, 1)),
          static_cast<float>(c(4, 2)), static_cast<float>(c(4, 3)),
          static_cast<float>(c(4, 4))};
    }
  };

  for (const auto& state : track.trackStates()) {
    auto typeFlags = state.typeFlags();
    if (!typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
      continue;
    }

    edm4hep::TrackState& trackState = outTrackStates.emplace_back();
    trackState.location = edm4hep::TrackState::AtOther;

    // This makes the hard assumption that |q| = 1
    SingleBoundTrackParameters<SinglyCharged> params{
        state.referenceSurface().getSharedPtr(), state.parameters(),
        state.covariance()};

    // Convert to LCIO track parametrization expected by EDM4hep
    detail::Parameters converted =
        detail::convertTrackParametersToEdm4hep(gctx, Bz, params);

    // Write the converted parameters to the EDM4hep track state
    setParameters(trackState, converted);

    // Converted parameters are relative to an ad-hoc perigee surface created at
    // the hit location
    auto center = converted.surface->center(gctx);
    trackState.referencePoint.x = center.x();
    trackState.referencePoint.y = center.y();
    trackState.referencePoint.z = center.z();
  }
  outTrackStates.front().location = edm4hep::TrackState::AtLastHit;
  outTrackStates.back().location = edm4hep::TrackState::AtFirstHit;

  // Add a track state that represents the IP parameters
  auto& ipState = outTrackStates.emplace_back();

  // Convert the track parameters at the IP
  SingleBoundTrackParameters<SinglyCharged> trackParams{
      track.referenceSurface().getSharedPtr(), track.parameters(),
      track.covariance()};

  // Convert to LCIO track parametrization expected by EDM4hep
  auto converted =
      detail::convertTrackParametersToEdm4hep(gctx, Bz, trackParams);
  setParameters(ipState, converted);
  ipState.location = edm4hep::TrackState::AtIP;

  // Write the converted parameters to the EDM4hep track state
  // The reference point is at the location of the reference surface of the
  // track itself, but if that's not a perigee surface, another ad-hoc perigee
  // at the position will be created.
  auto center = converted.surface->center(gctx);
  ipState.referencePoint.x = center.x();
  ipState.referencePoint.y = center.y();
  ipState.referencePoint.z = center.z();

  for (auto& trackState : outTrackStates) {
    to.addToTrackStates(trackState);
  }
}
}  // namespace EDM4hepUtil
}  // namespace Acts
