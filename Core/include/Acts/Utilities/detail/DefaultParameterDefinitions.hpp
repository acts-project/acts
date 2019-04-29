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
#include "Acts/Utilities/ParameterTypes.hpp"
#include "Acts/Utilities/Definitions.hpp"

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
  TrackParsDim /// < The local dimensions
};

/// The dimensions of tracks in global coordinates
constexpr unsigned int GlobalParsDim = 7;
/// The dimension of a space point
constexpr unsigned int SpacePointDim = 3;

using ParID_t    = ParDef;
using ParValue_t = double;

///
/// Type namings with local parameters
///

/// Vector of local track parameters
using TrackVector         = ActsVector<ParValue_t, TrackParsDim>;
/// Row vector of local track parameters
using TrackRowVector      = ActsRowVector<ParValue_t, TrackParsDim>;
/// Matrix of local-to-local parameters
using TrackMatrix         = ActsMatrix<ParValue_t, TrackParsDim, TrackParsDim>;
/// Symmetical matrix of local-to-local parameters
using TrackSymMatrix      = ActsSymMatrix<ParValue_t, TrackParsDim>;

///
/// Type naming with global parameters
///

/// Vector of global track parameters
using GlobalVector        = ActsVector<ParValue_t, GlobalParsDim>;
/// Matrix of global-to-global parameters
using GlobalMatrix = ActsMatrix<ParValue_t, GlobalParsDim, GlobalParsDim>;

///
/// Type namings with local & global parameters
///

/// Matrix of local-to-global parameters
using TrackToGlobalMatrix = ActsMatrix<ParValue_t, GlobalParsDim, TrackParsDim>;
/// Matrix of global-to-local parameters
using GlobalToTrackMatrix = ActsMatrix<ParValue_t, TrackParsDim, GlobalParsDim>;

///
/// Type namings with space points
///

/// Vector with space point parameters
using SpacePointVector = ActsVector<ParValue_t, SpacePointDim>;
/// Symmetrical matrix of space point-to-space point
using SpacePointSymMatrix = ActsSymMatrix<ParValue_t, SpacePointDim>;

///
/// Type namings with space points & local parameters
///

/// Matrix of space point-to-local parameters
using SpacePointToTrackMatrix = ActsMatrix<ParValue_t, TrackParsDim, SpacePointDim>;
/// Matrix with local parameters-to-space point
using TrackToSpacePointMatrix = ActsMatrix<ParValue_t, SpacePointDim, TrackParsDim>;

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
}  // namespace Acts