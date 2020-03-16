// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <type_traits>

#include "Acts/Utilities/Definitions.hpp"

#ifdef ACTS_PARAMETER_DEFINITIONS_PLUGIN
#include ACTS_PARAMETER_DEFINITIONS_PLUGIN
#endif

namespace Acts {

/// typedefs for parameter identifier and parameter value must be present
static_assert(std::is_enum<Acts::ParID_t>::value,
              "'ParID_t' is not an enum type");
static_assert(std::is_floating_point<Acts::ParValue_t>::value,
              "'ParValue_t' is not floating point type");

/// parameter ID type must be convertible to size_t
static_assert(std::is_convertible<Acts::ParID_t, size_t>::value,
              "'ParID_t' is not convertible to size_t");

/// number of global parameter must be at least 2 (for the two local parameters)
static_assert(Acts::BoundParsDim > 1,
              "total number of global parameters must be >= 2");

/// several constants for the local parameters need to be defined
static_assert(Acts::eLOC_0 != Acts::eLOC_1,
              "local parameters must have different IDs");
static_assert(Acts::eLOC_R == Acts::eLOC_0 or Acts::eLOC_R == Acts::eLOC_1,
              "local radius must be a local parameter");
static_assert(Acts::eLOC_PHI == Acts::eLOC_0 or Acts::eLOC_PHI == Acts::eLOC_1,
              "local phi must be a local parameter");
static_assert(Acts::eLOC_RPHI == Acts::eLOC_0 or
                  Acts::eLOC_RPHI == Acts::eLOC_1,
              "local r x phi must be a local parameter");
static_assert(Acts::eLOC_Z == Acts::eLOC_0 or Acts::eLOC_Z == Acts::eLOC_1,
              "local z must be a local parameter");
static_assert(Acts::eLOC_X == Acts::eLOC_0 or Acts::eLOC_X == Acts::eLOC_1,
              "local x must be a local parameter");
static_assert(Acts::eLOC_Y == Acts::eLOC_0 or Acts::eLOC_Y == Acts::eLOC_1,
              "local y must be a local parameter");
static_assert(Acts::eLOC_D0 == Acts::eLOC_0 or Acts::eLOC_D0 == Acts::eLOC_1,
              "d0 must be a local parameter");
static_assert(Acts::eLOC_Z0 == Acts::eLOC_0 or Acts::eLOC_Z0 == Acts::eLOC_1,
              "z0 must be a local parameter");

/// check for par_type_t definition
static_assert(sizeof(Acts::par_type_t<Acts::eLOC_0>) > 0,
              "'par_type_t' is not defined");

// Matrix and vector types related to bound track parameters.

using BoundVector = ActsVector<BoundParametersScalar, eBoundParametersSize>;
using BoundRowVector =
    ActsRowVector<BoundParametersScalar, eBoundParametersSize>;
using BoundMatrix = ActsMatrix<BoundParametersScalar, eBoundParametersSize,
                               eBoundParametersSize>;
using BoundSymMatrix =
    ActsSymMatrix<BoundParametersScalar, eBoundParametersSize>;

// Matrix and vector types related to free track parameters.

using FreeVector = ActsVector<FreeParametersScalar, eFreeParametersSize>;
using FreeRowVector = ActsRowVector<FreeParametersScalar, eFreeParametersSize>;
using FreeMatrix =
    ActsMatrix<FreeParametersScalar, eFreeParametersSize, eFreeParametersSize>;
using FreeSymMatrix = ActsSymMatrix<FreeParametersScalar, eFreeParametersSize>;

// Matrix and vector types related to space points.

using SpacePointVector = ActsVector<SpacePointScalar, eSpacePointSize>;
using SpacePointRowVector = ActsRowVector<SpacePointScalar, eSpacePointSize>;
using SpacePointSymMatrix =
    ActsMatrix<SpacePointScalar, eSpacePointSize, eSpacePointSize>;
using SpacePointSymMatrix = ActsSymMatrix<SpacePointScalar, eSpacePointSize>;

// Mapping to bound track parameters.
//
// Assumes that matrices represent maps from another space into the space of
// bound track parameters. Thus, the bound parameters scalar type is sufficient
// to retain accuracy.

using FreeToBoundMatrix = ActsMatrix<BoundParametersScalar,
                                     eBoundParametersSize, eFreeParametersSize>;
using SpacePointToBoundMatrix =
    ActsMatrix<BoundParametersScalar, eBoundParametersSize, eSpacePointSize>;

// Mapping to free track parameters.
//
// Assumes that matrices represent maps from another space into the space of
// free track parameters. Thus, the free parameters scalar type is sufficient
// to retain accuracy.

using BoundToFreeMatrix =
    ActsMatrix<FreeParametersScalar, eFreeParametersSize, eBoundParametersSize>;
using SpacePointToFreeMatrix =
    ActsMatrix<FreeParametersScalar, eFreeParametersSize, eSpacePointSize>;

// Mapping to space points.
//
// Assumes that matrices represent maps from another space into the space point
// space. Thus, the space point scalar type is sufficient to retain accuracy.

using BoundToSpacePointMatrix =
    ActsMatrix<SpacePointScalar, eSpacePointSize, eBoundParametersSize>;
using FreeToSpacePointMatrix =
    ActsMatrix<SpacePointScalar, eSpacePointSize, eFreeParametersSize>;

}  // namespace Acts
