// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts::detail {

/// Outcome of one attempted Runge-Kutta step
///
/// Returned by the generated kernels instead of an @c Acts::Result, whose
/// variant is read back through the stack on the accepted path.
enum class Rk4Status {
  /// the error estimate exceeded the tolerance; nothing was written
  Rejected,
  /// the step was taken and the outputs written
  Accepted,
  /// a field lookup failed; the error is in the `fieldErr` output
  FieldError,
};

}  // namespace Acts::detail
