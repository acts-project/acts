// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/DenseMaterialInteractor.hpp"

#include "Acts/Material/Interactions.hpp"

void Acts::DenseMaterialInteractor::reset(Acts::ActsScalar q,
                                          Acts::ActsScalar momentum,
                                          Acts::ActsScalar time,
                                          Vector3 direction,
                                          result_type& result) const {
  result.s = 0.;
  result.qop0 = (q == 0. ? 1. : q) / momentum;
  result.t0 = time;
  result.sInX0 = 0.;
  result.dir0 = direction;
}

Acts::BoundSymMatrix Acts::DenseMaterialInteractor::evaluateMultipleScattering(
    int pdg, float mass, float q, Acts::ActsScalar momentum,
    Acts::ActsScalar time, result_type& result) const {
  /// The content of this method follows the derivations in
  /// ATL-SOFT-PUB-2008-003 Eq. (18). While in this matrix the scattering
  /// contribution from a given layer expressed at a point at a path length d
  /// away from the actual interaction, the formalism is modified in this
  /// context. The actual contribution is evaluated within this actor at the end
  /// of the volume material or whenever a covariance transport occured. Hence,
  /// d becomes zero and some addends vanish. Since this paper did not consider
  /// time propagation, the corresponding term has been added using the
  /// derivation for q/p (named lambda in the paper).

  float qOverP = result.qop0;
  float x0 = result.sInX0;

  // Evaluate the standard deviation in theta
  const auto theta0 = computeMultipleScatteringTheta0(x0, pdg, mass, qOverP, q);
  const auto varTheta0 = theta0 * theta0;
  const ActsScalar s2 = result.s * result.s;
  const ActsScalar invSinTheta =
      result.dir0.norm() / VectorHelpers::perp(result.dir0);

  const ActsScalar qop = (q == 0. ? 1. : q) / momentum;
  const ActsScalar deltaQop = qop - result.qop0;
  const ActsScalar deltaT = time - result.t0;

  // Build the covariance matrix
  BoundSymMatrix msCovariance = BoundSymMatrix::Zero();
  msCovariance(eBoundLoc0, eBoundLoc0) = s2 / 3.;
  msCovariance(eBoundLoc0, eBoundPhi) = result.s * invSinTheta * 0.5;
  msCovariance(eBoundLoc1, eBoundLoc1) = s2 / 3.;
  msCovariance(eBoundLoc1, eBoundTheta) = -result.s * 0.5;
  msCovariance(eBoundTheta, eBoundLoc0) = result.s * invSinTheta * 0.5;
  msCovariance(eBoundPhi, eBoundPhi) = invSinTheta * invSinTheta;
  msCovariance(eBoundTheta, eBoundLoc1) = -result.s * 0.5;
  msCovariance(eBoundTheta, eBoundTheta) = 1.;
  msCovariance(eBoundQOverP, eBoundQOverP) =
      3. * deltaQop * deltaQop * varTheta0;
  msCovariance(eBoundTime, eBoundTime) = 3. * deltaT * deltaT * varTheta0;
  msCovariance *= varTheta0;

  return msCovariance;
}