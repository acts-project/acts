// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/CombinatorialSeedSolver.hpp"

namespace Acts::Experimental::detail {

template <Experimental::CompositeSpacePointPtr Point_t>
std::array<std::size_t, 3> separateLayers(
    const std::array<Point_t, 3>& layerTriplet) {
  // make sure we have exacxtly one 2D and two 1D measurements available
  auto n2D = std::ranges::count_if(layerTriplet, [](const auto& sp) {
    unsigned int dimension = sp->measuresLoc0() + sp->measuresLoc1();
    return dimension == 2;
  });

  auto n1D = std::ranges::count_if(layerTriplet, [](const auto& sp) {
    unsigned int dimension = sp->measuresLoc0() + sp->measuresLoc1();
    return dimension == 1;
  });

  if (n2D != 1 || n1D != 2) {
    throw std::runtime_error(
        "Triplet must contain exactly one 2D and two 1D space points");
  }

  std::array<std::size_t, 3> retIndices{
      std::numeric_limits<std::size_t>::max()};

  for (std::size_t idx = 0; idx < layerTriplet.size(); ++idx) {
    const Point_t& spacePoint = layerTriplet[idx];
    unsigned int dimension =
        spacePoint->measuresLoc0() + spacePoint->measuresLoc1();
    if (dimension == 1) {
      retIndices[0] == std::numeric_limits<std::size_t>::max()
          ? retIndices[0] = idx
          : retIndices[1] = idx;
    } else {
      retIndices[2] = idx;
    }
  }

  return retIndices;
}

}  // namespace Acts::Experimental::detail

namespace Acts::Experimental::CombinatorialSeedSolver {

using Acts::Experimental::detail::separateLayers;

template <Experimental::CompositeSpacePointPtr Point_t>
SquareMatrix2 betaMatrix(const std::array<Point_t, 4>& layerQuartett) {
  SquareMatrix2 bMatrix{SquareMatrix2::Identity()};

  Vector3 b1Aterm = (layerQuartett[3]->sensorDirection().dot(
                        layerQuartett[1]->sensorDirection())) *
                        layerQuartett[1]->sensorDirection() -
                    layerQuartett[3]->sensorDirection();
  Vector3 b1Gterm = (layerQuartett[3]->sensorDirection().dot(
                        layerQuartett[0]->sensorDirection())) *
                        layerQuartett[0]->sensorDirection() -
                    (layerQuartett[3]->sensorDirection().dot(
                        layerQuartett[1]->sensorDirection())) *
                        layerQuartett[1]->sensorDirection();

  Vector3 b2Aterm = layerQuartett[2]->sensorDirection() -
                    (layerQuartett[2]->sensorDirection().dot(
                        layerQuartett[1]->sensorDirection())) *
                        layerQuartett[1]->sensorDirection();
  Vector3 b2Kterm = (layerQuartett[2]->sensorDirection().dot(
                        layerQuartett[1]->sensorDirection())) *
                        layerQuartett[1]->sensorDirection() -
                    (layerQuartett[2]->sensorDirection().dot(
                        layerQuartett[0]->sensorDirection())) *
                        layerQuartett[0]->sensorDirection();

  // get the distances of the layers along z direction
  double A = (layerQuartett[0]->localPosition().z() -
              layerQuartett[1]->localPosition().z());
  double G = (layerQuartett[0]->localPosition().z() -
              layerQuartett[2]->localPosition().z());
  double K = (layerQuartett[0]->localPosition().z() -
              layerQuartett[3]->localPosition().z());

  // define B matrix
  Vector3 b1 = A * b1Aterm + G * b1Gterm;
  Vector3 b2 = K * b2Kterm + A * b2Aterm;

  bMatrix.col(0) = Vector2(b1.x(), b1.y());
  bMatrix.col(1) = Vector2(b2.x(), b2.y());

  return bMatrix;
}

template <Experimental::CompositeSpacePointPtr Point_t>
std::array<double, 4> defineParameters(
    const SquareMatrix2& betaMatrix,
    const std::array<Point_t, 4>& layerQuartett) {
  double A = (layerQuartett[0]->localPosition().z() -
              layerQuartett[1]->localPosition().z());
  double G = (layerQuartett[0]->localPosition().z() -
              layerQuartett[2]->localPosition().z());
  double K = (layerQuartett[0]->localPosition().z() -
              layerQuartett[3]->localPosition().z());

  // Define y2 for the linear system
  Vector3 y0 = K * (layerQuartett[2]->localPosition() -
                    layerQuartett[0]->localPosition()) +
               G * (layerQuartett[0]->localPosition() -
                    layerQuartett[3]->localPosition());
  Vector3 y1 = A * (layerQuartett[3]->localPosition() -
                    layerQuartett[2]->localPosition()) +
               G * (layerQuartett[1]->localPosition() -
                    layerQuartett[3]->localPosition()) +
               K * (layerQuartett[2]->localPosition() -
                    layerQuartett[1]->localPosition());
  Vector3 y2 = (K - G) * (layerQuartett[0]->localPosition() -
                          layerQuartett[1]->localPosition()) -
               (y1.dot(layerQuartett[1]->sensorDirection())) *
                   layerQuartett[1]->sensorDirection() +
               (y0.dot(layerQuartett[0]->sensorDirection())) *
                   layerQuartett[0]->sensorDirection() +
               A * (layerQuartett[3]->localPosition() -
                    layerQuartett[2]->localPosition());

  Vector2 solution = betaMatrix.inverse() * y2.block<2, 1>(0, 0);
  double kappa = solution.x();
  double gamma = solution.y();

  double lambda = (y0.dot(layerQuartett[0]->sensorDirection()) +
                   K * gamma *
                       (layerQuartett[0]->sensorDirection().dot(
                           layerQuartett[2]->sensorDirection())) -
                   G * kappa *
                       (layerQuartett[0]->sensorDirection().dot(
                           layerQuartett[3]->sensorDirection()))) /
                  (K - G);
  double alpha = (y1.dot(layerQuartett[1]->sensorDirection()) +
                  (A - G) * kappa *
                      layerQuartett[3]->sensorDirection().dot(
                          layerQuartett[1]->sensorDirection()) +
                  (K - A) * gamma *
                      layerQuartett[2]->sensorDirection().dot(
                          layerQuartett[1]->sensorDirection())) /
                 (K - G);

  return std::array<double, 4>({lambda, alpha, gamma, kappa});
}

template <Experimental::CompositeSpacePointPtr Point_t>
std::pair<Vector3, Vector3> seedSolution(
    const std::array<Point_t, 4>& layerQuartett,
    const std::array<double, 4>& parameters) {
  // estimate the seed positionInChamber from the 1st equation of the system of
  // the layer equations
  Vector3 seedPosition = layerQuartett[0]->localPosition() +
                         parameters[0] * layerQuartett[0]->sensorDirection();

  // estimate the seed direction from the 2nd equation of the system of the
  // layer equations
  Vector3 seedDirection =
      ((layerQuartett[1]->localPosition() +
        parameters[1] * layerQuartett[1]->sensorDirection() - seedPosition))
          .normalized();

  // calculate the position at z=0
  Intersection3D intersectionZ0 = PlanarHelper::intersectPlane(
      seedPosition, seedDirection, Vector3::UnitZ(), 0.0);
  Vector3 seedPositionZ0 = intersectionZ0.position();

  return std::make_pair(seedPositionZ0,
                        copySign(seedDirection, seedDirection.z()));
}

template <Experimental::CompositeSpacePointPtr Point_t>
SquareMatrix2 betaMatrix(const std::array<Point_t, 3>& layerTriplet) {
  SquareMatrix2 bMatrix{SquareMatrix2::Identity()};

  // keep the 1D spacepoints and the 2D spacepoint seperately- for invariance
  // under the order of the layers
  std::array<std::size_t, 3> indices = separateLayers(layerTriplet);
  Point_t spacePoint2D{layerTriplet[indices[2]]};
  std::array<Point_t, 2> spacePoints1D{layerTriplet[indices[0]],
                                       layerTriplet[indices[1]]};

  double R =
      spacePoint2D->localPosition().z() - spacePoints1D[0]->localPosition().z();
  double L =
      spacePoint2D->localPosition().z() - spacePoints1D[1]->localPosition().z();
  double dotProd = spacePoints1D[0]->sensorDirection().dot(
      spacePoints1D[1]->sensorDirection());

  // construct the betaMatrix
  bMatrix.col(0) = Vector2(L, L * dotProd);
  bMatrix.col(1) = Vector2(-R * dotProd, -R);

  return bMatrix;
}

template <Experimental::CompositeSpacePointPtr Point_t>
std::array<double, 2> defineParameters(
    const SquareMatrix2& betaMatrix,
    const std::array<Point_t, 3>& layerTriplet) {
  std::array<std::size_t, 3> indices = separateLayers(layerTriplet);
  const Point_t& spacePoint2D{layerTriplet[indices[2]]};
  std::array<Point_t, 2> spacePoints1D{layerTriplet[indices[0]],
                                       layerTriplet[indices[1]]};

  const Vector3& M = spacePoint2D->localPosition();
  double R =
      spacePoint2D->localPosition().z() - spacePoints1D[0]->localPosition().z();
  double L =
      spacePoint2D->localPosition().z() - spacePoints1D[1]->localPosition().z();

  Vector3 y = R * spacePoints1D[1]->localPosition() -
              L * spacePoints1D[0]->localPosition() + (L - R) * M;
  Vector2 y2 = {y.dot(spacePoints1D[0]->sensorDirection()),
                y.dot(spacePoints1D[1]->sensorDirection())};

  Vector2 solution = betaMatrix.inverse() * y2.block<2, 1>(0, 0);
  double beta = solution.x();
  double delta = solution.y();

  return std::array<double, 2>({beta, delta});
}

template <Experimental::CompositeSpacePointPtr Point_t>
std::pair<Vector3, Vector3> seedSolution(
    const std::array<Point_t, 3>& layerTriplet,
    const std::array<double, 2>& parameters) {
  // separate 2D and 1D spacepoints
  std::array<std::size_t, 3> indices = separateLayers(layerTriplet);
  Point_t spacePoint2D{layerTriplet[indices[2]]};
  std::array<Point_t, 2> spacePoints1D{layerTriplet[indices[0]],
                                       layerTriplet[indices[1]]};

  // the position of the seed can be evaluated from the 2D measurement
  const Vector3& seedPosition = spacePoint2D->localPosition();
  // The direction of the seed can be estimated from the second layer equation
  Vector3 seedDirection =
      (spacePoints1D[0]->localPosition() +
       parameters[0] * spacePoints1D[0]->sensorDirection() - seedPosition)
          .normalized();
  // calculate the position at z=0
  Intersection3D intersectionZ0 = PlanarHelper::intersectPlane(
      seedPosition, seedDirection, Vector3::UnitZ(), 0.0);
  Vector3 seedPositionZ0 = intersectionZ0.position();

  return std::make_pair(seedPositionZ0,
                        copySign(seedDirection, seedDirection.z()));
}

}  // namespace Acts::Experimental::CombinatorialSeedSolver
