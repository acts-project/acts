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
#include "Acts/Utilities/MathHelpers.hpp"

#include <concepts>
#include <iostream>
#include <ranges>
#include <vector>

namespace Acts::Experimental::detail {

/// Helper function that splits the layers' spacepoints between the 2D and two
/// 1D ones when three layer combinatorics are available
///@tparam Point_t the space point type
///@param layerTriplet the space points of the combinatorics
///@return The tuple with the two 1D space points and one 2D spacepoint

template <typename Point_t>
std::tuple<std::array<Point_t, 2>, Point_t> separateLayers(
    const std::array<Point_t, 3>& layerTriplet);

}  // namespace Acts::Experimental::detail

namespace Acts::Experimental::CombinatorialSeedSolver {

/// A Combinatorial Seed Solver for seed estimation from combinatoric hits from
/// four or three layers
// (e.g Muon NSW seeding for ATLAS) with overloaded functions to implement the
// mathematics
/// The combinatoric layers are expected to be sorted in z in chamber's frame

/// ===============================
/// The 4-layer combinatorics functions
/// ===============================

///***The layers equations for the four layers (Si,Di) can be */
/// S1 + lambda*D1 = M
/// S2 + alpha*D2 = M + A*Dm
/// S3 + gamma*D3 = M + G*Dm
/// S4 + kappa*D4 = M + K*Dm
/// where {Si,Di}  are the position and direction of the strip
/// {M,Dm} the position and direction of the muon's trajectory on the 1st plane
/// solving the system we can reduce it to a 2x2 system for kappa and gamma --->
/// https://gitlab.cern.ch/atlas-nextgen/work-package-2.5/analyticalsegment
/// Bx1*kappa + Bx2*gamma = Y2x
/// By1*kappa + By2*gamma = Y2y

/// Defines the betaMatrix calculated from the combinatoric hits
/// @tparam Point_t the space point type
/// @param layerQuartett the space points of the combinatorics
/// @return the 2x2 beta matrix
template <Experimental::CompositeSpacePointPtr Point_t>
SquareMatrix2 betaMatrix(const std::array<Point_t, 4>& layerQuartett);

// Calculates the parameters lambda,alpha,gamma,kappa of the system
/// @tparam Point_t the space point type
/// @param betaMatrix the betaMatrix for the system
/// @param layerQuartett the space points of the combinatorics
/// @return an array of the calculated parameters
template <Experimental::CompositeSpacePointPtr Point_t>
std::array<double, 4> defineParameters(
    const SquareMatrix2& betaMatrix,
    const std::array<Point_t, 4>& layerQuartett);

/// Solves the equation system to calculate the seed position and direction
/// @tparam spacePointContainr the space point container
/// @param layerQuartett the space points of the combinatorics
/// @param parameters the lambda,alpha,gamma,kappa parameters of the four layers
/// @return the pair of the seed position and direction
template <Experimental::CompositeSpacePointPtr Point_t>
std::pair<Vector3, Vector3> seedSolution(
    const std::array<Point_t, 4>& layerQuartett,
    const std::array<double, 4>& parameters);

/// ===============================
/// The 3-layer combinatorics functions
/// ===============================

/// A combinatorial seed solver from three available layers where one 2D
/// measurement and two 1D measurement are available(e.g STgc from NSW)
///***The layers equations for the three layers (Si,Di) can be */
/// P = M //from the 2D measurement
/// S1 + beta*D1 = M + R*Dm //1st layer of the 1D measurement
/// S2 + delta*D2 = M + L*Dm //2nd layer of the 1D measurement
/// where {Si,Di}  are the position and direction of the strip
/// {M,Dm} the position and direction of the muon's trajectory on the 1st plane
/// solving the system we can reduce it to a 2x2 system for beta and delta --->
/// https://gitlab.cern.ch/atlas-nextgen/work-package-2.5/analyticalsegment
/// Bx1*beta + Bx2*delta = Yx
/// By1*beta + By2*delta = Yy

/// Defines the components of the beta matrix for a triplet of combinatorics
/// @tparam Point_t the space point type
/// @param layerTriplet the space points from the combinatorics
/// @return The 2x2 beta matrix
template <Experimental::CompositeSpacePointPtr Point_t>
SquareMatrix2 betaMatrix(const std::array<Point_t, 3>& layerTriplet);

/// Calculates the parameters beta, delta of the system of the three layers
/// combinatorics
/// @tparam Point_t the space point type
/// @param betaMatrix the betaMatrix for the system
/// @param layerTriplet the space points of the combinatorics
/// @return an array of the calculated parameters
template <Experimental::CompositeSpacePointPtr Point_t>
std::array<double, 2> defineParameters(
    const SquareMatrix2& betaMatrix,
    const std::array<Point_t, 3>& layerTriplet);

/// Solve the equations for the seed position and direction (M,DM) for the three
/// layers case
/// @tparam Point_t the spacepoint type
/// @param layerTriplet the spacepoints of the combinatorics
/// @param parameters the beta,delta parameters of the system
template <Experimental::CompositeSpacePointPtr Point_t>
std::pair<Vector3, Vector3> seedSolution(
    const std::array<Point_t, 3>& layerTriplet,
    const std::array<double, 2>& parameters);

}  // namespace Acts::Experimental::CombinatorialSeedSolver

#include "Acts/Seeding/CombinatorialSeedSolver.ipp"
