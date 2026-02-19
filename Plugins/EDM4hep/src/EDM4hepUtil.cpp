// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/EDM4hep/EDM4hepUtil.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"
#include "Acts/Propagator/detail/JacobianEngine.hpp"

#include <numbers>

#include <edm4hep/EDM4hepVersion.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/MutableSimTrackerHit.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/TrackState.h>

using namespace Acts;
using namespace Acts::detail;

namespace ActsPlugins::EDM4hepUtil {
namespace detail {

SquareMatrix<6> jacobianToEdm4hep(double theta, double qOverP, double Bz) {
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
  // [- csc (theta)           0     ]
  // [                              ]
  // [-B*q/p*cos(theta)       B     ]
  // [------------------  ----------]
  // [      2             sin(theta)]
  // [   sin (theta)                ]

  SquareMatrix<6> J;
  J.setIdentity();
  double sinTheta = std::sin(theta);
  J(3, 3) = -1.0 / (sinTheta * sinTheta);
  J(4, 4) = Bz / sinTheta;  // dOmega / d(qop)
  J(4, 3) = -Bz * qOverP * std::cos(theta) /
            (sinTheta * sinTheta);  // dOmega / dTheta
  return J;
}

SquareMatrix<6> jacobianFromEdm4hep(double tanLambda, double omega, double Bz) {
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

  SquareMatrix<6> J;
  J.setIdentity();
  J(3, 3) = -1 / (tanLambda * tanLambda + 1);
  J(4, 3) = -1 * omega * tanLambda /
            (Bz * std::pow(tanLambda * tanLambda + 1, 3. / 2.));
  J(4, 4) = 1 / (Bz * std::hypot(tanLambda, 1));

  return J;
}

void packCovariance(const SquareMatrix<6>& from, float* to) {
  for (int i = 0; i < from.rows(); i++) {
    for (int j = 0; j <= i; j++) {
      std::size_t k = (i + 1) * i / 2 + j;
      to[k] = from(i, j);
    }
  }
}

void unpackCovariance(const float* from, SquareMatrix<6>& to) {
  auto k = [](std::size_t i, std::size_t j) { return (i + 1) * i / 2 + j; };
  for (int i = 0; i < to.rows(); i++) {
    for (int j = 0; j < to.cols(); j++) {
      to(i, j) = from[j <= i ? k(i, j) : k(j, i)];
    }
  }
}

Parameters convertTrackParametersToEdm4hep(const GeometryContext& gctx,
                                           double Bz,
                                           const BoundTrackParameters& params) {
  Vector3 global = params.referenceSurface().localToGlobal(
      gctx, params.parameters().template head<2>(), params.direction());

  std::shared_ptr<const Surface> refSurface =
      params.referenceSurface().getSharedPtr();
  BoundVector targetPars = params.parameters();
  std::optional<BoundMatrix> targetCov = params.covariance();

  // If the reference surface is a perigee surface, we use that. Otherwise
  // we create a new perigee surface at the global position of the track
  // parameters.
  if (dynamic_cast<const PerigeeSurface*>(refSurface.get()) == nullptr) {
    refSurface = Surface::makeShared<PerigeeSurface>(global);

    // We need to convert to the target parameters
    // Keep the free parameters around we might need them for the covariance
    // conversion

    auto perigeeParams =
        boundToBoundConversion(gctx, params, *refSurface, Vector3{0, 0, Bz})
            .value();
    targetPars = perigeeParams.parameters();
    targetCov = perigeeParams.covariance();
  }

  Parameters result;
  result.surface = refSurface;

  // Only run covariance conversion if we have a covariance input
  if (targetCov) {
    SquareMatrix<6> J = jacobianToEdm4hep(targetPars[eBoundTheta],
                                          targetPars[eBoundQOverP], Bz);
    result.covariance = J * targetCov.value() * J.transpose();
  }

  result.values[0] = targetPars[eBoundLoc0];
  result.values[1] = targetPars[eBoundLoc1];
  result.values[2] = targetPars[eBoundPhi];
  result.values[3] = std::tan(std::numbers::pi / 2. - targetPars[eBoundTheta]);
  result.values[4] =
      targetPars[eBoundQOverP] / std::sin(targetPars[eBoundTheta]) * Bz;
  result.values[5] = targetPars[eBoundTime];

  result.particleHypothesis = params.particleHypothesis();

  return result;
}

BoundTrackParameters convertTrackParametersFromEdm4hep(
    double Bz, const Parameters& params) {
  BoundVector targetPars;

  SquareMatrix<6> J =
      jacobianFromEdm4hep(params.values[3], params.values[4], Bz);

  std::optional<BoundMatrix> cov;
  if (params.covariance.has_value()) {
    cov = J * params.covariance.value() * J.transpose();
  }

  targetPars[eBoundLoc0] = params.values[0];
  targetPars[eBoundLoc1] = params.values[1];
  targetPars[eBoundPhi] = params.values[2];
  targetPars[eBoundTheta] = std::numbers::pi / 2. - std::atan(params.values[3]);
  targetPars[eBoundQOverP] =
      params.values[4] * std::sin(targetPars[eBoundTheta]) / Bz;
  targetPars[eBoundTime] = params.values[5];

  return {params.surface, targetPars, cov, params.particleHypothesis};
}

}  // namespace detail

#if EDM4HEP_VERSION_MAJOR >= 1 || \
    (EDM4HEP_VERSION_MAJOR == 0 && EDM4HEP_VERSION_MINOR == 99)
edm4hep::MCParticle getParticle(const edm4hep::SimTrackerHit& hit) {
  return hit.getParticle();
}

void setParticle(edm4hep::MutableSimTrackerHit& hit,
                 const edm4hep::MCParticle& particle) {
  hit.setParticle(particle);
}
#else
edm4hep::MCParticle getParticle(const edm4hep::SimTrackerHit& hit) {
  return hit.getMCParticle();
}

void setParticle(edm4hep::MutableSimTrackerHit& hit,
                 const edm4hep::MCParticle& particle) {
  hit.setMCParticle(particle);
}
#endif

std::size_t SimHitAssociation::size() const {
  return m_internalToEdm4hep.size();
}

void SimHitAssociation::reserve(std::size_t size) {
  m_internalToEdm4hep.reserve(size);
}

void SimHitAssociation::add(std::size_t internalIndex,
                            const edm4hep::SimTrackerHit& edm4hepHit) {
  m_internalToEdm4hep.push_back(edm4hepHit);
  // m_edm4hepToInternal.at(edm4hepHit.id()) = internalIndex;
  m_edm4hepToInternal.emplace(edm4hepHit.id(), internalIndex);
}

edm4hep::SimTrackerHit SimHitAssociation::lookup(
    std::size_t internalIndex) const {
  return m_internalToEdm4hep.at(internalIndex);
}

std::size_t SimHitAssociation::lookup(const edm4hep::SimTrackerHit& hit) const {
  return m_edm4hepToInternal.at(hit.id());
}

}  // namespace ActsPlugins::EDM4hepUtil
