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
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/CovarianceTransport.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <edm4hep/Track.h>
#include <edm4hep/TrackState.h>

#include "edm4hep/MutableTrack.h"

namespace Acts {
namespace EDM4hepUtil {

static constexpr std::int32_t EDM4HEP_ACTS_POSITION_TYPE = 42;

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

  auto setParameters = [Bz](edm4hep::TrackState& trackState, const auto& params,
                            const auto& cov) {
    // Conversion:
    // https://bib-pubdb1.desy.de/record/81214/files/LC-DET-2006-004%5B1%5D.pdf
    trackState.D0 = params[Acts::eBoundLoc0];
    trackState.Z0 = params[Acts::eBoundLoc1];
    trackState.phi = params[Acts::eBoundPhi];
    trackState.tanLambda = std::tan(M_PI_2 - params[Acts::eBoundTheta]);
    trackState.time = params[Acts::eBoundTime];

    double p = SinglyCharged{}.extractMomentum(params[Acts::eBoundQOverP]);
    double pt = std::sin(params[Acts::eBoundTheta]) * p;
    double q = SinglyCharged{}.extractCharge(params[Acts::eBoundQOverP]);

    trackState.omega = q * Bz / pt;

    // Calculate jacobian from our internal parametrization (d0, z0, phi, theta,
    // q/p) to the LCIO / edm4hep one (d0, z0, phi, tan(lambda), omega). Top
    // left 3x3 matrix in the jacobian is 1.
    // Bottom right 2x2 matrix is:
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
    double cotTheta = std::tan(M_PI_2 + params[Acts::eBoundTheta]);
    J(3, 3) = -cotTheta * cotTheta - 1;  // d(tanLambda) / dTheta
    J(4, 4) = Bz / std::sin(params[Acts::eBoundTheta]);  // dOmega / d(qop)
    double sinTheta = std::sin(params[eBoundTheta]);
    J(4, 3) = -Bz * params[Acts::eBoundQOverP] *
              std::cos(params[eBoundTheta] /
                       (sinTheta * sinTheta));  // dOmega / dTheta

    Acts::ActsSymMatrix<5> cIn = cov.template topLeftCorner<5, 5>();
    Acts::ActsSymMatrix<5> cOut = J * cIn * J.transpose();

    trackState.covMatrix = {
        static_cast<float>(cOut(0, 0)), static_cast<float>(cOut(1, 0)),
        static_cast<float>(cOut(1, 1)), static_cast<float>(cOut(2, 0)),
        static_cast<float>(cOut(2, 1)), static_cast<float>(cOut(2, 2)),
        static_cast<float>(cOut(3, 0)), static_cast<float>(cOut(3, 1)),
        static_cast<float>(cOut(3, 2)), static_cast<float>(cOut(3, 3)),
        static_cast<float>(cOut(4, 0)), static_cast<float>(cOut(4, 1)),
        static_cast<float>(cOut(4, 2)), static_cast<float>(cOut(4, 3)),
        static_cast<float>(cOut(4, 4))};
  };

  for (const auto& state : track.trackStates()) {
    auto typeFlags = state.typeFlags();
    if (!typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
      continue;
    }

    // EDM4hep wants IP parameters relative to the surface center
    Acts::BoundVector params = state.parameters();
    Acts::BoundMatrix cov = state.covariance();

    const Acts::Surface& refSurface = state.referenceSurface();

    Acts::Vector3 direction =
        makeDirectionUnitFromPhiTheta(params[eBoundPhi], params[eBoundTheta]);
    // This makes the hard assumption that |q| = 1
    ActsScalar absMomentum =
        SinglyCharged{}.extractMomentum(params[eBoundQOverP]);
    Acts::Vector3 momentum = direction * absMomentum;

    Acts::Vector3 global =
        refSurface.localToGlobal(gctx, params.head<2>(), momentum);

    auto pseudoPerigee =
        Acts::Surface::makeShared<Acts::PerigeeSurface>(global);

    auto boundToFree = pseudoPerigee->boundToFreeJacobian(gctx, params);

    Acts::FreeVector freePars =
        Acts::detail::transformBoundToFreeParameters(refSurface, gctx, params);
    Acts::FreeMatrix freeCov = boundToFree * cov * boundToFree.transpose();

    Acts::BoundVector targetPars = Acts::detail::transformFreeToBoundParameters(
                                       freePars, *pseudoPerigee, gctx)
                                       .value();

    Acts::CovarianceCache covCache{freePars, freeCov};
    auto [varNewCov, varNewJac] = Acts::transportCovarianceToBound(
        gctx, *pseudoPerigee, freePars, covCache);
    auto targetCov = std::get<Acts::BoundSymMatrix>(varNewCov);

    edm4hep::TrackState& trackState = outTrackStates.emplace_back();
    trackState.location = edm4hep::TrackState::AtOther;

    setParameters(trackState, targetPars, targetCov);
    auto center = pseudoPerigee->center(gctx);
    trackState.referencePoint.x = center.x();
    trackState.referencePoint.y = center.y();
    trackState.referencePoint.z = center.z();
  }
  outTrackStates.front().location = edm4hep::TrackState::AtLastHit;
  outTrackStates.back().location = edm4hep::TrackState::AtFirstHit;

  // add a track state that represents the IP parameters
  auto& ipState = outTrackStates.emplace_back();
  setParameters(ipState, track.parameters(), track.covariance());
  ipState.location = edm4hep::TrackState::AtIP;
  ipState.referencePoint = {0, 0, 0};

  for (auto& trackState : outTrackStates) {
    to.addToTrackStates(trackState);
  }
}
}  // namespace EDM4hepUtil
}  // namespace Acts
