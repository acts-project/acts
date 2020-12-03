// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"

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

}  // namespace detail
}  // namespace Acts
