// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/CompositeSpacePoint.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <concepts>
#include <iostream>
#include <ranges>
#include <vector>

namespace Acts {
namespace CombinatorialSeedSolver {

/// @brief A Combinatorial Seed Solver for seed estimation from combinatoric hits from four layers (e.g Muon NSW seeding)

///***The layers equations for the four layers (Si,Di) can be */
/// S1 + lamda*D1 = M
/// S2 + alpha*D2 = M + A*Dm
/// S3 + gamma*D3 = M + G*Dm
/// S4 + kappa*D4 = M + K*Dm
/// where {Si,Di}  are the position and direction of the strip
/// {M,Dm} the position and direction of the muon's trajectory on the 1st plane
/// solving the system we can reduce it to a 2x2 system for kappa and gamma --->
/// https://gitlab.cern.ch/atlas-nextgen/work-package-2.5/analyticalsegment
/// Bx1*kappa + Bx2*gamma = Y2x
/// By1*kappa + By2*gamma = Y2y

/// @brief defines the betaMatrix calculated from the combinatoric hits
/// @tparam spacePointContainer the space point container
/// @param spacePoints the space points of the combinatorics
/// @return the 2x2 beta matrix
template <typename spacePointContainer>
  requires Acts::Experimental::CompositeSpacePointContainer<spacePointContainer>
SquareMatrix2 betaMatrix(const spacePointContainer& spacePoints) {
  SquareMatrix2 bMatrix{SquareMatrix2::Identity()};

  Vector3 b1Aterm = (spacePoints[3]->sensorDirection().dot(
                        spacePoints[1]->sensorDirection())) *
                        spacePoints[1]->sensorDirection() -
                    spacePoints[3]->sensorDirection();
  Vector3 b1Gterm = (spacePoints[3]->sensorDirection().dot(
                        spacePoints[0]->sensorDirection())) *
                        spacePoints[0]->sensorDirection() -
                    (spacePoints[3]->sensorDirection().dot(
                        spacePoints[1]->sensorDirection())) *
                        spacePoints[1]->sensorDirection();

  Vector3 b2Aterm = spacePoints[2]->sensorDirection() -
                    (spacePoints[2]->sensorDirection().dot(
                        spacePoints[1]->sensorDirection())) *
                        spacePoints[1]->sensorDirection();
  Vector3 b2Kterm = (spacePoints[2]->sensorDirection().dot(
                        spacePoints[1]->sensorDirection())) *
                        spacePoints[1]->sensorDirection() -
                    (spacePoints[2]->sensorDirection().dot(
                        spacePoints[0]->sensorDirection())) *
                        spacePoints[0]->sensorDirection();

  // get the distances of the layers along z direction
  double A = (spacePoints[0]->localPosition().z() -
              spacePoints[1]->localPosition().z());
  double G = (spacePoints[0]->localPosition().z() -
              spacePoints[2]->localPosition().z());
  double K = (spacePoints[0]->localPosition().z() -
              spacePoints[3]->localPosition().z());

  // define B matrix
  Vector3 b1 = A * b1Aterm + G * b1Gterm;
  Vector3 b2 = K * b2Kterm + A * b2Aterm;

  bMatrix.col(0) = Vector2(b1.x(), b1.y());
  bMatrix.col(1) = Vector2(b2.x(), b2.y());

  return bMatrix;
}

/// @brief calculates the parameters lamda,alpha,gamma,kappa of the system
/// @tparam spacePointContainer the space point container
/// @param betaMatrix the betaMatrix for the system
/// @param spacePoints the space points of the combinatorics
/// @return an array of the calculated parameters
template <typename spacePointContainer>
  requires Acts::Experimental::CompositeSpacePointContainer<spacePointContainer>
std::array<double, 4> defineParameters(SquareMatrix2 betaMatrix,
                                       const spacePointContainer& spacePoints) {
  double A = (spacePoints[0]->localPosition().z() -
              spacePoints[1]->localPosition().z());
  double G = (spacePoints[0]->localPosition().z() -
              spacePoints[2]->localPosition().z());
  double K = (spacePoints[0]->localPosition().z() -
              spacePoints[3]->localPosition().z());

  // Define y2 for the linear system
  Vector3 y0 =
      K * (spacePoints[2]->localPosition() - spacePoints[0]->localPosition()) +
      G * (spacePoints[0]->localPosition() - spacePoints[3]->localPosition());
  Vector3 y1 =
      A * (spacePoints[3]->localPosition() - spacePoints[2]->localPosition()) +
      G * (spacePoints[1]->localPosition() - spacePoints[3]->localPosition()) +
      K * (spacePoints[2]->localPosition() - spacePoints[1]->localPosition());
  Vector3 y2 =
      (K - G) *
          (spacePoints[0]->localPosition() - spacePoints[1]->localPosition()) -
      (y1.dot(spacePoints[1]->sensorDirection())) *
          spacePoints[1]->sensorDirection() +
      (y0.dot(spacePoints[0]->sensorDirection())) *
          spacePoints[0]->sensorDirection() +
      A * (spacePoints[3]->localPosition() - spacePoints[2]->localPosition());

  Vector2 solution = betaMatrix.inverse() * y2.block<2, 1>(0, 0);
  double kappa = solution.x();
  double gamma = solution.y();

  double lamda = (y0.dot(spacePoints[0]->sensorDirection()) +
                  K * gamma *
                      (spacePoints[0]->sensorDirection().dot(
                          spacePoints[2]->sensorDirection())) -
                  G * kappa *
                      (spacePoints[0]->sensorDirection().dot(
                          spacePoints[3]->sensorDirection()))) /
                 (K - G);
  double alpha = (y1.dot(spacePoints[1]->sensorDirection()) +
                  (A - G) * kappa *
                      spacePoints[3]->sensorDirection().dot(
                          spacePoints[1]->sensorDirection()) +
                  (K - A) * gamma *
                      spacePoints[2]->sensorDirection().dot(
                          spacePoints[1]->sensorDirection())) /
                 (K - G);

  return std::array<double, 4>({lamda, alpha, gamma, kappa});
}

// Solve the equations for the seed position and direction (M,DM)
/// @brief solves the equation system to calculate the seed
/// @tparam spacePointContainr the space point container
/// @param spacePoints the space points of the combinatorics
/// @param parameters the lamda,alpha,gamma,kappa paramaters of the four layers
/// @return the pair of the seed position and direction
template <typename spacePointContainer>
  requires Acts::Experimental::CompositeSpacePointContainer<spacePointContainer>
std::pair<Vector3, Vector3> seedSolution(
    const spacePointContainer& spacePoints,
    const std::array<double, 4>& parameters) {
  // estimate the seed positionInChamber from the 1st equation of the system of
  // the layer equations
  Vector3 seedPosition = spacePoints[0]->localPosition() +
                         parameters[0] * spacePoints[0]->sensorDirection();

  // estimate the seed direction from the 2nd equation of the system of the
  // layer equations
  Vector3 seedDirection =
      ((spacePoints[1]->localPosition() +
        parameters[1] * spacePoints[1]->sensorDirection() - seedPosition))
          .normalized();

  // calculate the position at z=0
  Intersection3D intersectionZ0 = PlanarHelper::intersectPlane(
      seedPosition, seedDirection, Vector3::UnitZ(), 0.0);
  Vector3 seedPositionZ0 = intersectionZ0.position();

  return std::make_pair(seedPositionZ0,
                        seedDirection.z() > 0 ? seedDirection : -seedDirection);
};

}  // namespace CombinatorialSeedSolver

}  // namespace Acts