// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <tuple>
#include <utility>
#include <vector>

namespace Acts::detail {

// The container for derivative of local frame axis w.r.t. its
// rotation parameters. The first element is for x axis, second for y axis and
// last for z axis
using RotationToAxes =
    std::tuple<RotationMatrix3, RotationMatrix3, RotationMatrix3>;

/// @brief Evaluate the derivative of local frame axes vector w.r.t.
/// its rotation around global x/y/z axis
/// @Todo: add parameter for rotation axis order
///
/// @param rotation The rotation that help place the surface
///
/// @return Derivative of local frame x/y/z axis vector w.r.t. its
/// rotation angles (extrinsic Euler angles) around global x/y/z axis
RotationToAxes rotationToLocalAxesDerivative(const RotationMatrix3& rotation);

}  // namespace Acts::detail
