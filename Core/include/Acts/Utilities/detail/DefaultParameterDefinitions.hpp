// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
// STL include(s)
#include <cmath>

// Acts includes
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterTypes.hpp"

namespace Acts {
enum ParDef : unsigned int {
  eLOC_0 = 0,  ///< first coordinate in local surface frame
  eLOC_1 = 1,  ///< second coordinate in local surface frame
  eLOC_R = eLOC_0,
  eLOC_PHI = eLOC_1,
  eLOC_RPHI = eLOC_0,
  eLOC_Z = eLOC_1,
  eLOC_X = eLOC_0,
  eLOC_Y = eLOC_1,
  eLOC_D0 = eLOC_0,
  eLOC_Z0 = eLOC_1,
  ePHI = 2,    ///< phi direction of momentum in global frame
  eTHETA = 3,  ///< theta direction of momentum in global frame
  eQOP = 4,    ///< charge/momentum for charged tracks, for neutral tracks it is
               /// 1/momentum
  eT = 5,      /// < The time of the particle
  BoundParsDim  /// < The local dimensions
};

/// The dimensions of tracks in free coordinates
constexpr unsigned int FreeParsDim = 8;

/// The dimension of a space point
constexpr unsigned int SpacePointDim = 4;

using ParID_t = ParDef;
using ParValue_t = double;

///
/// Type namings with bound parameters
///

/// Vector of bound parameters
using BoundVector = ActsVector<ParValue_t, BoundParsDim>;
/// Row vector of bound parameters
using BoundRowVector = ActsRowVector<ParValue_t, BoundParsDim>;
/// Matrix of bound-to-bound parameters
using BoundMatrix = ActsMatrix<ParValue_t, BoundParsDim, BoundParsDim>;
/// Symmetrical matrix of bound-to-bound parameters
using BoundSymMatrix = ActsSymMatrix<ParValue_t, BoundParsDim>;

///
/// Type naming with free parameters
///

/// Vector of free track parameters
using FreeVector = ActsVector<ParValue_t, FreeParsDim>;
/// Matrix of free-to-free parameters
using FreeMatrix = ActsMatrix<ParValue_t, FreeParsDim, FreeParsDim>;
/// Symmetrical matrix of free-to-free parameters
using FreeSymMatrix = ActsSymMatrix<ParValue_t, FreeParsDim>;

///
/// Type namings with bound & free parameters
///

/// Matrix of bound-to-free parameters
using BoundToFreeMatrix = ActsMatrix<ParValue_t, FreeParsDim, BoundParsDim>;
/// Matrix of free-to-bound parameters
using FreeToBoundMatrix = ActsMatrix<ParValue_t, BoundParsDim, FreeParsDim>;

///
/// Type namings with space points
///

/// Vector with space point parameters
using SpacePointVector = ActsVector<ParValue_t, SpacePointDim>;
/// Symmetrical matrix of space point-to-space point
using SpacePointSymMatrix = ActsSymMatrix<ParValue_t, SpacePointDim>;

///
/// Type namings with space points & bound parameters
///

/// Matrix of space point-to-bound parameters
using SpacePointToBoundMatrix =
    ActsMatrix<ParValue_t, BoundParsDim, SpacePointDim>;
/// Matrix with bound parameters-to-space point
using BoundToSpacePointMatrix =
    ActsMatrix<ParValue_t, SpacePointDim, BoundParsDim>;

template <ParID_t>
struct par_type;

template <ParID_t par>
using par_type_t = typename par_type<par>::type;

template <>
struct par_type<ParDef::eLOC_0> {
  using type = local_parameter;
};

template <>
struct par_type<ParDef::eLOC_1> {
  using type = local_parameter;
};

template <>
struct par_type<ParDef::ePHI> {
  static constexpr double pMin() { return -M_PI; }
  static constexpr double pMax() { return M_PI; }
  using type = cyclic_parameter<double, pMin, pMax>;
};

template <>
struct par_type<ParDef::eTHETA> {
  static constexpr double pMin() { return 0; }
  static constexpr double pMax() { return M_PI; }
  using type = bound_parameter<double, pMin, pMax>;
};

template <>
struct par_type<ParDef::eQOP> {
  using type = unbound_parameter;
};

template <>
struct par_type<ParDef::eT> {
  using type = unbound_parameter;
};
}  // namespace Acts