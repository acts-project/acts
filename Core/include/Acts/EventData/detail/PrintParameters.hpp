// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"

#include <cstdint>
#include <iosfwd>

namespace Acts {

class Surface;

namespace detail {

/// Print bound track parameters content to the output stream.
///
/// @param os The output stream
/// @param surface Bound parameters reference surface
/// @param params Bound parameters vector
/// @param cov Optional bound parameters covariance matrix
void printBoundParameters(std::ostream& os, const Surface& surface,
                          const BoundVector& params,
                          const BoundSymMatrix* cov = nullptr);

/// Print free track parameters content to the output stream.
///
/// @param os The output stream
/// @param params Free parameters vector
/// @param cov Optional free parameters covariance matrix
void printFreeParameters(std::ostream& os, const FreeVector& params,
                         const FreeMatrix* cov = nullptr);

/// Print bound measurement content to the output stream.
///
/// @param os The output stream
/// @param size Size of the measurement space
/// @param indices Which parameters are measured, must contain size elements
/// @param params Parameters vector
/// @param cov Optional Covariance matrix
void printMeasurement(std::ostream& os, BoundIndices size,
                      const uint8_t* indices,
                      const Eigen::Ref<const ActsVectorX<BoundScalar>>& params,
                      const Eigen::Ref<const ActsMatrixX<BoundScalar>>& cov);

/// Print free measurement content to the output stream.
///
/// @param os The output stream
/// @param size Size of the measurement space
/// @param indices Which parameters are measured, must contain size elements
/// @param params Parameters vector
/// @param cov Optional Covariance matrix
void printMeasurement(std::ostream& os, FreeIndices size,
                      const uint8_t* indices,
                      const Eigen::Ref<const ActsVectorX<FreeScalar>>& params,
                      const Eigen::Ref<const ActsMatrixX<FreeScalar>>& cov);

}  // namespace detail
}  // namespace Acts
